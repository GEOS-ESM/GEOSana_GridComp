#include <params.h>
      subroutine diagref(     
     &  lchnk,  pcols,                               
     &  pbot,   ubot,   vbot,   tbot,  qbot,  thbot,  zbot,     
     &  taux,   tauy,   shf,    lhf,   ts,    oro,             
     &  t10m,   q10m,   u10m,   v10m,  t2m,   q2m,    u2m,    v2m)

!----------------------------------------------------------------------- 
! 
! Purpose:  diagnostics of 10 meter and 2 meter fields
!
! Method:  Retrieve diagnostic mass and momentum fields from known 
!          fluxes and model lowest level fields based on similarity theory
! Assume:
!   1) Similarity of sub-grid tiles over land can be represented 
!      by mean variables
!   2) assume displacement is equal to zero
!   3) coefficients for similarity are the same as those used in ccm
! 
! Author: Jiundar Chern, Jan. 2001
! 
!-----------------------------------------------------------------------
      use precision

      implicit none

#include <pmgrid.h>

#include <comcon.h>

!
!-----------------------------------------------------------------------
! note: Constants cpvir,cpair,gravit,latvap,latice,rair, and zvir are 
!       defined in comcon.h
!  cpvir    = cpwv/cpair -1
!  cpwv     = Specific heat of water vapor
!  cpair    = Specific heat of dry air
!  gravit   = Gravitational acceleration
!  latvap   = Latent heat of vaporization
!  latice   = Latent heat of fusion
!  rair     = Gas constant for dry air
!  zvir     = rh2o/rair - 1
!  rh2o     = Gas constant for water vapor
!
!-----------------------------------------------------------------------
!
! Input arguments
!
      integer,  intent(in) :: lchnk            ! chunk indentifer
      integer,  intent(in) :: pcols            ! longitude dimension

      real(r8), intent(in) :: pbot(pcols)      ! Bottom level pressure
      real(r8), intent(in) :: ubot(pcols)      ! Bottom level u wind
      real(r8), intent(in) :: vbot(pcols)      ! Bottom level v wind
      real(r8), intent(in) :: tbot(pcols)      ! Bottom level temperature
      real(r8), intent(in) :: qbot(pcols)      ! Bottom level specific humidity
      real(r8), intent(in) :: thbot(pcols)     ! Bottom level potential temperature
      real(r8), intent(in) :: zbot(pcols)      ! Bottom level height above surface

      real(r8), intent(in) :: taux(pcols)      ! x surface stress (zonal)
      real(r8), intent(in) :: tauy(pcols)      ! y surface stress (meridional)
      real(r8), intent(in) :: shf(pcols)       ! sensible heat flux (W/m2)
      real(r8), intent(in) :: lhf(pcols)       ! latent heat flux (W/m2)
      real(r8), intent(in) :: ts(pcols)        ! surface temperature (K)
      real(r8), intent(in) :: oro(pcols)       ! land-sea mask 
!
! Output arguments
!
      real(r8), intent(out) :: t10m(pcols)     ! 10 meter temperature
      real(r8), intent(out) :: q10m(pcols)     ! 10 meter specify humidity
      real(r8), intent(out) :: u10m(pcols)     ! 10 meter U wind
      real(r8), intent(out) :: v10m(pcols)     ! 10 meter V wind

      real(r8), intent(out) :: t2m(pcols)      ! 2 meter temperature
      real(r8), intent(out) :: q2m(pcols)      ! 2 meter specify humidity
      real(r8), intent(out) :: u2m(pcols)      ! 2 meter U wind
      real(r8), intent(out) :: v2m(pcols)      ! 2 meter V wind
!
!---------------------------Local variable-----------------------------
!
      integer i                                ! Longitude indices
      real(r8), parameter :: z10   = 10.0      ! reference height at 10 meter
      real(r8), parameter :: z2    =  2.0      ! reference height at 2 meter
      real(r8), parameter :: xkar   = 0.4      ! Von Karman constant

      real(r8) th10m                           ! 10 meter potential temperature
      real(r8) th2m                            ! 10 meter potential temperature
    
      real(r8) rbot                            ! Density at bottom model level
      real(r8) cp                              ! Specific heat of moist air
      real(r8) ltheat                          ! Latent heat for given srf conditions
      real(r8) tau                             ! Surface stress
      real(r8) ustar                           ! u*
      real(r8) tstar                           ! t*
      real(r8) qstar                           ! q*
      real(r8) xsq                             ! Temporary variable
      real(r8) xqq                             ! Temporary variable
      real(r8) alz10                           ! ln(z10/zbot)
      real(r8) alz2                            ! ln(z2/zbot)
      real(r8) hol                             ! hgt at bottom model level / monin-obukhov length
      real(r8) hol10                           ! Ref hgt (10m) / monin-obukhov length
      real(r8) hol2                            ! Ref hgt (2m) / monin-obukhov length
      real(r8) psimh                           ! Stability function for momentum at bottom model level
      real(r8) psimh10                         ! Stability function for momentum at 10 meter
      real(r8) psimh2                          ! Stability function for momentum at 2 meter
      real(r8) psixh                           ! Stability function for heat & tracers at bottom model level
      real(r8) psixh10                         ! Stability function for heat & tracers at 10 meter
      real(r8) psixh2                          ! Stability function for heat & tracers at 2 meter
      real(r8) stable                          ! Stability factor
      real(r8) vmag                            ! wind magnitude at bottom model level
      real(r8) vmag10                          ! wind magnitude at 10 meter
      real(r8) vmag2                           ! wind magnitude at 2 meter
      real(r8) tmin                            ! min(tbot,ts)
      real(r8) tmax                            ! max(tbot,ts)
      logical  violate                         ! if true: violate the similarity
!
!--------------------------Statement functions--------------------------
!
      real(r8) psimhu                          ! Unstable part of psimh
      real(r8) psixhu                          ! Unstable part of psixh
      real(r8) xd                              ! Dummy argument
!
      psimhu(xd) = log((1.+xd*(2.+xd))*(1.+xd*xd)/8.) - 2.*atan(xd) + 1.571
      psixhu(xd) = 2. * log((1. + xd*xd)/2.)
!-----------------------------------------------------------------------

!
! Loop over Longitude
!
      do i = 1, pcols
!
!---------------------------------------------------------------
! Set up necessary variables
!---------------------------------------------------------------
!
         rbot    = pbot(i) / (rair*tbot(i))
         cp      = cpair*( 1. + cpvir*qbot(i) ) 
         if( tbot(i) > 273.15 ) then
           ltheat = latvap
         else
           ltheat = latvap + latice
         end if
!
!---------------------------------------------------------------
! Compute ustar, tstar, and qstar from taux, tauy, shf ,and lhf
!    tau     = sqrt( taux*taux+tauy*tauy )
!    ustar   = sqrt( tau/rbot )
!    tstar   = -shf/( cp*rbot*ustar )
!    qstar   = -lhf/( ltheat*rbot*ustar )
!---------------------------------------------------------------
!
         tau     = sqrt( taux(i) * taux(i) + tauy(i) * tauy(i) )
         ustar   = sqrt( tau / rbot )
         tstar   = - shf(i) / ( cp * rbot * ustar )
         qstar   = - lhf(i) / ( ltheat * rbot * ustar )
!
!---------------------------------------------------------------
! Compute stability and  stability function for momentum at model lowest level and 10 m
!    hol   = xkar *gravit *zbot * (tstar/thbot + qstar/(1./zvir+qbot(i))) / ustar**2
!    hol10 = hol*z10/zbot
!    hol2  = hol*z2/zbot
!---------------------------------------------------------------
!
         hol     = xkar *gravit *zbot(i) * ( tstar/thbot(i) + 
     &             qstar/(1./zvir+qbot(i)) ) / (ustar * ustar)
         hol     = sign( min(abs(hol),10.), hol )
         hol10   = hol * z10 / zbot(i)
         hol2    = hol * z2 / zbot(i)
         stable  = 0.5 + sign(0.5_r8 , hol)

         xsq     = max(sqrt(abs(1. - 16.*hol)) , 1.)
         xqq     = sqrt(xsq)
         psimh   = -5. * hol * stable + (1.-stable)*psimhu(xqq)
         psixh   = -5. * hol * stable + (1.-stable)*psixhu(xqq)

         xsq     = max(sqrt(abs(1. - 16.*hol10)) , 1.)
         xqq     = sqrt(xsq)
         psimh10 = -5. * hol10 * stable + (1.-stable)*psimhu(xqq)
         psixh10 = -5. * hol10 * stable + (1.-stable)*psixhu(xqq)

         xsq     = max(sqrt(abs(1. - 16.*hol2)) , 1.)
         xqq     = sqrt(xsq)
         psimh2  = -5. * hol2 * stable + (1.-stable)*psimhu(xqq)
         psixh2  = -5. * hol2 * stable + (1.-stable)*psixhu(xqq)
!
!---------------------------------------------------------------
! Compute wind at 10 and 2 meter
! vmag   = sqrt( ubot*ubot +vbot*vbot )
! vmag10 = vmag + ustar/xkar * ( alz10-psimh10+psimh )
! u10m   = vmag10 * ubot / vmag
! v10m   = vmag10 * vbot / vmag
! vmag2  = vmag + ustar/xkar * ( alz2-psimh2+psimh )
! u2m    = vmag2 * ubot / vmag
! v2m    = vmag2 * vbot / vmag
!---------------------------------------------------------------
!
         vmag    = max( 1.0, sqrt(ubot(i) * ubot(i) + vbot(i) * vbot(i)) )
         alz10   = log( z10 / zbot(i) ) 
         vmag10  = vmag + ustar * ( alz10 - psimh10 + psimh ) / xkar
         alz2    = log( z2 / zbot(i) ) 
         vmag2   = vmag + ustar * ( alz2 - psimh2 + psimh ) / xkar
!
         violate = .false.
 
! stability function violate the similarity
         if( vmag10 .lt. 0.0  .or. vmag2 .lt. 0.0 ) violate = .true.
! low wind and very stable condition stability function violate the similarity
         if( vmag .eq. 1.0  .and. hol .gt. 5.0 ) violate = .true.
         if( hol .eq. 10.0 ) violate = .true.
         if( violate ) then
           vmag10  = max( 0.0, vmag + ustar * alz10 / xkar)
           vmag2  = max( 0.0, vmag + ustar * alz2 / xkar )
         end if
           
         u10m(i) = vmag10 * ubot(i) / vmag
         v10m(i) = vmag10 * vbot(i) / vmag
         u2m(i)  = vmag2 * ubot(i) / vmag
         v2m(i)  = vmag2 * vbot(i) / vmag
!
!---------------------------------------------------------------
! Compute potential temperature and mixing ratio at 10 and 2 meter
! th10m  = thbot + tstar/xkar * ( alz10-psixh10+psixh )
! q10m   = qbot + qstar/xkar * ( alz10-psixh10+psixh )
! th2m   = thbot + tstar/xkar * ( alz2-psixh2+psixh )
! q2m    = qbot + qstar/xkar * ( alz2-psixh2+psixh )
!---------------------------------------------------------------
!
         if( violate ) then
           th10m     = thbot(i) + tstar/xkar *  alz10
           q10m(i)   = qbot(i) + qstar/xkar *  alz10
           th2m      = thbot(i) + tstar/xkar *  alz2
           q2m(i)    = qbot(i) + qstar/xkar *  alz2
         else
           th10m     = thbot(i) + tstar/xkar * ( alz10-psixh10+psixh )
           q10m(i)   = qbot(i) + qstar/xkar * ( alz10-psixh10+psixh )
           th2m      = thbot(i) + tstar/xkar * ( alz2-psixh2+psixh )
           q2m(i)    = qbot(i) + qstar/xkar * ( alz2-psixh2+psixh )
         end if
!
!---------------------------------------------------------------
! Compute temperature at 10 and 2 meter
! dT/T   = dth/th  + rd/cp * dP/P
! t10m   = th10m * tbot / thbot + g * (zbot-z10) / cp / (1.0+zvir* qbot)
! t2m    = th2m * tbot / thbot + g * (zbot-z2) / cp / (1.0+zvir* qbot)
!---------------------------------------------------------------
!
         t10m(i)   = th10m * tbot(i) / thbot(i) 
     &    + gravit * ( zbot(i) - z10 ) / cpair / (1.0 + zvir * qbot(i))
         t2m(i)    = th2m * tbot(i) / thbot(i)   
     &    + gravit * ( zbot(i) - z2 ) / cpair / (1.0 + zvir * qbot(i))

!---------------------------------------------------------------
! monotonic constrains
!---------------------------------------------------------------
         tmin=min(ts(i),tbot(i))-1.0
         tmax=max(ts(i),thbot(i))+1.0

         if( t2m(i) .lt. tmin .or. t2m(i) .gt. tmax ) then
           t10m(i)  = tbot(i) + (ts(i)-tbot(i))* (1.0-10.0/zbot(i))
           t2m(i)  = tbot(i) + (ts(i)-tbot(i))* (1.0-2.0/zbot(i))

           if( hol .gt. 0.0 ) then
             q10m(i)   = qbot(i) + qstar/xkar * alz10
             q2m(i)    = qbot(i) + qstar/xkar * alz2
             vmag10  = max( 0.0, vmag + ustar * alz10 / xkar)
             vmag2  = max( 0.0, vmag + ustar * alz2 / xkar )
             u10m(i) = vmag10 * ubot(i) / vmag
             v10m(i) = vmag10 * vbot(i) / vmag
             u2m(i)  = vmag2 * ubot(i) / vmag
             v2m(i)  = vmag2 * vbot(i) / vmag
           end if
         end if
      end do

      return
      end subroutine diagref
