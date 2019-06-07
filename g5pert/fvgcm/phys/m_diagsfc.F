!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_diagsfc --- diagnostics of surface fluxes and 2-meter and 10-meter fields
!
! !INTERFACE:
!
      MODULE  m_diagsfc

! !USES:

!      Implicit NONE

!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  diagsfc           ! diagnostics of surface fluxes and 2-meter and 10-meter 
                                ! fields from analysis and diag datasets
!
! !DESCRIPTION: This module computes the diagnostics of surface fluxes and 2-meter 
!               and 10-meter fields from analysis and diag datasets
!
! !REVISION HISTORY:
!
! 10Jan2002  j Chern  Initial specs and prologues.
!
!EOP
!-------------------------------------------------------------------------

      CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  diagsfc -- diagnostics of surface fluxes and 2-meter and 10-meter fields
!
! !INTERFACE:
!
      subroutine diagsfc( im,      jm,                           
     &                   pbot,    tbot,    qbot,   ubot,   vbot,
     &                   zbot,    lwi,     ts,     ps,     bma,
     &                   z0h,     z0m,     zpd,               
     &                   shf,     lhf,     taux,   tauy,     
     &                   t2m,     q2m,     u2m,    v2m,     
     &                   t10m,    q10m,    u10m,   v10m    )
!
! !USES:
!
      use precision
      implicit none 
#include <comcon.h>
!
! !INPUT PARAMETERS:
!
      integer,  intent(in)  :: im               !  zonal dimension
      integer,  intent(in)  :: jm               !  meridional dimension

      real(r8), intent(in)  ::  pbot(im,jm)     ! Bottom level pressure (Pa)
      real(r8), intent(in)  ::  tbot(im,jm)     ! Bottom level temperature (K)
      real(r8), intent(in)  ::  qbot(im,jm)     ! Bottom level specific humidity (kg/kg)
      real(r8), intent(in)  ::  ubot(im,jm)     ! Bottom level U-wind (m/s)
      real(r8), intent(in)  ::  vbot(im,jm)     ! Bottom level V-wind (m/s)
      real(r8), intent(in)  ::  zbot(im,jm)     ! Bottom level height (m)
      real(r8), intent(in)  ::  lwi(im,jm)      ! Land-water-ice mask:
                                                !   lwi = 0  over ocean
                                                !   lwi = 1  over land
                                                !   lwi = 2  over sea ice
                                                ! NOTE: same as ORO in FVCCM.
      real(r8), intent(in)  ::  ts(im,jm)       ! Surface temperature (K)
      real(r8), intent(in)  ::  ps(im,jm)       ! Surface pressure (Pa)
      real(r8), intent(in)  ::  bma(im,jm)      ! Bulk moisture availability (fraction)
      real(r8), intent(in)  ::  z0h(im,jm)      ! Roughness length, sensible heat (m)
      real(r8), intent(in)  ::  z0m(im,jm)      ! Roughness length, momentum (m)
      real(r8), intent(in)  ::  zpd(im,jm)      ! Displacement height (m)

!
! !OUTPUT PARAMETERS:
!
      real(r8), intent(out)  :: shf(im,jm)      ! sensible heat flux (W/m2)
      real(r8), intent(out)  :: lhf(im,jm)      ! latent heat flux (W/m2)
      real(r8), intent(out)  :: taux(im,jm)     ! X surface stress (N/m2)
      real(r8), intent(out)  :: tauy(im,jm)     ! Y surface stress (N/m2)
      real(r8), intent(out)  :: t2m(im,jm)      ! 2-meter temperature (K)
      real(r8), intent(out)  :: q2m(im,jm)      ! 2-meter specific humidity (kg/kg)
      real(r8), intent(out)  :: u2m(im,jm)      ! 2-meter u-wind (m/s)
      real(r8), intent(out)  :: v2m(im,jm)      ! 2-meter v-wind (m/s)
      real(r8), intent(out)  :: t10m(im,jm)     ! 10-meter temperature (K)
      real(r8), intent(out)  :: q10m(im,jm)     ! 10-meter specific humidity (kg/kg)
      real(r8), intent(out)  :: u10m(im,jm)     ! 10-meter u-wind (m/s)
      real(r8), intent(out)  :: v10m(im,jm)     ! 10-meter v-wind (m/s)

! !DESCRIPTION: Given the near surface analysis fields and diag fields
!               to compute surface fluxes and 2-meter and 10-meter fields
!               based on similarity theory.
!
! !REVISION HISTORY:
!
! 10Jan2001  Chern     initial code
!
!EOP
!-------------------------------------------------------------------------

!--------------------------Local Variables------------------------------
!
      integer i,j               ! column indices
      integer npts              ! Number of land/ocean/ice points
      integer indx(im)          ! column index array (land/ocean/ice)
      real(r8) oro(im)          ! land-sea-ice mask
      real(r8) ltheat(im)       ! Latent heat for given sfc conditions
      real(r8) thbot(im,jm)     ! Latent heat for given sfc conditions
      real(r8) psomc(im,jm)     ! (ps/pbot)**cappa

!--------------------------Statement Functions--------------------------
      logical ocean
      logical land  
      logical seaice
      ocean(i)  = nint(oro(i))==0
      land(i)   = nint(oro(i))==1
      seaice(i) = nint(oro(i))==2

!-----------------------------------------------------------------------
! This progrm is called offline, therefore, we need to set all of constants
! in comcon.h 
!-------------------------------------------------------------------------
      cpair  = 1004.64
      epsilo = 0.622
      gravit = 9.80616
      latvap = 2.5104e06
      latice = 3.336e5

      omega  = 7.292e-5
      rair   = 287.04
      rearth = 6.37122e+6
      rhoh2o = 1.e3

      stebol = 5.67e-8
      tmelt  = 273.16

      cappa  = rair/cpair
      ez     = omega/sqrt(0.375)
      ra     = 1./rearth
      rga    = 1./gravit

      rh2o   = 4.61e2
      zvir   = rh2o/rair - 1.
      cpwv   = 1.81e3
      cpvir  = cpwv/cpair - 1.
!
      do j = 1, jm
        do i = 1, im
          psomc(i,j) = (ps(i,j)/pbot(i,j))**cappa
          thbot(i,j) = tbot(i,j) * psomc(i,j)
        end do
      end do

      do j = 1, jm

        do i = 1, im
          oro(i) = lwi(i,j)
          if( ts(i,j) > 273.15 ) then
            ltheat(i) = latvap
          else
            ltheat(i) = latvap + latice
          end if
          if( ocean(i) ) ltheat(i) = latvap
          if( seaice(i) ) ltheat(i) = latvap + latice
        end do
!
! Set up index array of land surfaces
!
        npts = 0
        do i = 1, im
          if( land(i) ) then
            npts       = npts + 1
            indx(npts) = i
          end if
        end do
        if( npts > 0) then
!
! Compute surface fluxes, 2-meter and 10-meter diagnostics over land
!
          call diaglnd(im,
     &         indx      ,npts      ,pbot(1,j)  ,ubot(1,j) ,vbot(1,j) ,
     &         tbot(1,j) ,qbot(1,j) ,thbot(1,j) ,zbot(1,j) ,ts(1,j)   ,
     &         ltheat    ,z0m(1,j)  ,z0h(1,j)   ,zpd(1,j)  ,bma(1,j)  , 
     &         shf(1,j)  ,lhf(1,j)  ,taux(1,j)  ,tauy(1,j) ,            
     &         t2m(1,j)  ,q2m(1,j)  ,u2m(1,j)   ,v2m(1,j)  ,            
     &         t10m(1,j) ,q10m(1,j) ,u10m(1,j)  ,v10m(1,j) )           
        end if
!
! Set up index array of ocean surfaces
!
        npts = 0
        do i = 1, im
          if( ocean(i) ) then
            npts       = npts + 1
            indx(npts) = i
          end if
        end do
        if( npts > 0) then
!
! Compute surface fluxes, 2-meter and 10-meter diagnostics over ocean
!
          call diagoce(im,
     &         indx      ,npts      ,pbot(1,j)  ,ubot(1,j) ,vbot(1,j) ,
     &         tbot(1,j) ,qbot(1,j) ,thbot(1,j) ,zbot(1,j) ,ts(1,j)   ,
     &         ltheat    ,shf(1,j)  ,lhf(1,j)   ,taux(1,j) ,tauy(1,j) ,
     &         t2m(1,j)  ,q2m(1,j)  ,u2m(1,j)   ,v2m(1,j)  ,          
     &         t10m(1,j) ,q10m(1,j) ,u10m(1,j)  ,v10m(1,j) )           
        end if
!
! Set up index array of seaice surfaces
!
        npts = 0
        do i = 1, im
          if( seaice(i) ) then
            npts       = npts + 1
            indx(npts) = i
          end if
        end do

        if( npts > 0) then
!
! Compute surface fluxes, 2-meter and 10-meter diagnostics over seaice
!
          call diagice(im,
     &         indx      ,npts      ,pbot(1,j)  ,ubot(1,j) ,vbot(1,j) ,
     &         tbot(1,j) ,qbot(1,j) ,thbot(1,j) ,zbot(1,j) ,ts(1,j)   ,
     &         ltheat    ,shf(1,j)  ,lhf(1,j)   ,taux(1,j) ,tauy(1,j) ,
     &         t2m(1,j)  ,q2m(1,j)  ,u2m(1,j)   ,v2m(1,j)  ,          
     &         t10m(1,j) ,q10m(1,j) ,u10m(1,j)  ,v10m(1,j) )           
        end if

      end do   ! end j loop

      return

      end subroutine diagsfc

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  diaglnd -- diagnostics of surface fluxes and 2-meter and 10-meter
!                        fields over land
!
! !INTERFACE:
!

      subroutine diaglnd(pcols    ,
     &           indx    ,npts    ,pmidm1  ,ubot    ,vbot    , 
     &           tbot    ,qbot    ,thbot   ,zbot    ,ts      ,
     &           ltheat  ,z0m     ,z0h     ,zpd     ,bma     ,
     &           shf     ,lhf     ,taux    ,tauy    ,        
     &           t2m     ,q2m     ,u2m     ,v2m     ,       
     &           t10m    ,q10m    ,u10m    ,v10m    )

!
! !USES:
!
      use precision
      implicit none
#include <comcon.h>
#include <parpbl.h>
!
! !INPUT PARAMETERS:
!
      integer , intent(in) :: pcols          ! number of column
      integer , intent(in) :: indx(pcols)    ! column index array (land)
      integer , intent(in) :: npts           ! Number of land points
      real(r8), intent(in) :: pmidm1(pcols)  ! Bottom level pressure
      real(r8), intent(in) :: ubot(pcols)    ! Bottom level u wind
      real(r8), intent(in) :: vbot(pcols)    ! Bottom level v wind
      real(r8), intent(in) :: tbot(pcols)    ! Bottom level temperature
      real(r8), intent(in) :: qbot(pcols)    ! Bottom level specific humidity
      real(r8), intent(in) :: thbot(pcols)   ! Bottom level potential temperature
      real(r8), intent(in) :: zbot(pcols)    ! Bottom level height above surface
      real(r8), intent(in) :: ts(pcols)      ! Surface temperature
      real(r8), intent(in) :: ltheat(pcols)  ! Latent heat for given srf conditions
      real(r8), intent(in) :: z0m(pcols)     ! Roughness length, momentum (m)
      real(r8), intent(in) :: z0h(pcols)     ! Roughness length, sensible heat (m)
      real(r8), intent(in) :: zpd(pcols)     ! Displacement height (m) 
      real(r8), intent(in) :: bma(pcols)     ! Bulk moisture availability (fraction)
    
!
! !OUTPUT PARAMETERS:
!
      real(r8), intent(out):: shf(pcols)     ! Initial sensible heat flux (W/m2)
      real(r8), intent(out):: lhf(pcols)     ! Initial latent heat flux (W/m2)
      real(r8), intent(out):: taux(pcols)    ! X surface stress (N/m2)
      real(r8), intent(out):: tauy(pcols)    ! Y surface stress (N/m2)
      real(r8), intent(out):: t2m(pcols)     ! 2-meter temperature (K)
      real(r8), intent(out):: q2m(pcols)     ! 2-meter specific humidity (kg/kg)
      real(r8), intent(out):: u2m(pcols)     ! 2-meter u-wind (m/s)
      real(r8), intent(out):: v2m(pcols)     ! 2-meter v-wind (m/s)
      real(r8), intent(out):: t10m(pcols)    ! 10-meter temperature (K)
      real(r8), intent(out):: q10m(pcols)    ! 10-meter specific humidity (kg/kg)
      real(r8), intent(out):: u10m(pcols)    ! 10-meter u-wind (m/s)
      real(r8), intent(out):: v10m(pcols)    ! 10-meter v-wind (m/s)
! !DESCRIPTION: Given the near surface analysis and diag fields to compute surface fluxes
!               and 2-meter and 10-meter fields over land based on similarity theory.
!
! Method: similarity theory, coefficients same as CCM
! Assume:
!   1) ignore subgrid (tile) heterogeneity, assume similarity theory works for bulk grid
!      mean.
!   2) Assume subgrid soil and vegetation resistances can be presented by a bulk
!      moisture availability
!   3) wind speeds should all be above a minimum speed (umin)
!
!-----------------------------------------------------------------------
! !REVISION HISTORY:
!
! 10Jan2001  Chern     initial code
!
!EOP

!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
      real(r8), parameter :: z2m  = 2.0 !  2-m height
      real(r8), parameter :: z10m = 10.0 !  10-m height
      real(r8), parameter :: crit = 0.01 !  criterion for converge
      integer i,ii            ! column indices
      integer nn              ! iteration #
      logical conv            ! iteration converge if .true.
      real(r8) ssq            ! Surface saturation specific humidity
      real(r8) qsfc           ! Surface specific humidity
      real(r8) ustar          ! ustar
      real(r8) tstar          ! tstar
      real(r8) qstar          ! qstar
      real(r8) vmag           ! Bottom level wind magnitude  (m/s)
      real(r8) vmag2          ! 2-meter wind magnitude (m/s)
      real(r8) vmag10         ! 10-meter wind magnitude (m/s)
      real(r8) thvbot         ! Bottom lev virtual potential temp
      real(r8) th2m           ! 2-meter potential temp (K)
      real(r8) th10m          ! 10-meter potential temp (K)
      real(r8) delt           ! potential T difference (K)
      real(r8) delq           ! specific humidity difference (kg/kg)
      real(r8) rdn            ! sqrt of neutral exchange coeff (momentum)
      real(r8) rhn            ! sqrt of neutral exchange coeff (heat)
      real(r8) ren            ! sqrt of neutral exchange coeff (tracers)          
      real(r8) rd             ! sqrt of exchange coeff (momentum)
      real(r8) rh             ! sqrt of exchange coeff (heat)
      real(r8) re             ! sqrt of exchange coeff (tracers)
      real(r8) rd0            ! old sqrt of exchange coeff (momentum)
      real(r8) rh0            ! old sqrt of exchange coeff (heat)
      real(r8) hol            ! (zbot -zpd) / monin-obukhov length
      real(r8) dz             ! zbot -zpd
      real(r8) z0mol          ! z0m / monin-obukhov length
      real(r8) z0hol          ! z0h / monin-obukhov length
      real(r8) z2ol           ! z2m / monin-obukhov length
      real(r8) z10ol          ! z10m / monin-obukhov length
      real(r8) xsq            ! Temporary variable
      real(r8) xqq            ! Temporary variable
      real(r8) alzm           ! ln( (zbot-zpd)/z0m )
      real(r8) alzh           ! ln( (zbot-zpd)/z0h )
      real(r8) alzm2          ! ln( (z2m)/z0m )
      real(r8) alzh2          ! ln( (z2m)/z0h )
      real(r8) alzm10         ! ln( (z10m)/z0m )
      real(r8) alzh10         ! ln( (z10m)/z0h )
      real(r8) cp             ! Specific heat of moist air
      real(r8) tau            ! Reference height stress
      real(r8) psimh          ! Stability function at bottom level (momentum)
      real(r8) psixh          ! Stability function at bottom level (heat & tracers) 
      real(r8) psimhz0        ! Stability function at z0m + zpd (momentum)
      real(r8) psixhz0        ! Stability function at z0h + zpd (heat & tracers) 
      real(r8) psimh2         ! Stability function at 2 meter (momentum)
      real(r8) psixh2         ! Stability function at 2 meter (heat & tracers) 
      real(r8) psimh10        ! Stability function at 10 meter (momentum)
      real(r8) psixh10        ! Stability function at 10 meter (heat & tracers) 
      real(r8) stable         ! Stability factor
      real(r8) rbot(pcols)    ! Density at bottom model level
      real(r8) del            ! working space
      real(r8) z0h1           ! z0h1=0.1*z0h
      real(r8) bma1           ! 
!-----------------------------------------------------------------------

!--------------------------Statement functions--------------------------
      real(r8) psimhu         ! Unstable part of psimh
      real(r8) psixhu         ! Unstable part of psixh
      real(r8) qsat           ! Saturation specific humidty of air
      real(r8) xd             ! Dummy argument
      real(r8) Tk             ! Temperature (K)
!
      qsat(Tk)   = 640380. / exp(5107.4/Tk)
      psimhu(xd) = log((1.+xd*(2.+xd))*(1.+xd*xd)/8.) - 2.*atan(xd) 
     &           + 1.571
      psixhu(xd) = 2. * log((1. + xd*xd)/2.)
!-----------------------------------------------------------------------
!
! Loop over land points
!
      do ii=1,npts
        i = indx(ii)
!
!---------------------------------------------------------------
! Set up necessary variables
!---------------------------------------------------------------
!
        conv   = .false.
        rbot(i)= pmidm1(i) / (rair*tbot(i))
        vmag   = max(umin, sqrt(ubot(i)**2+vbot(i)**2))
        thvbot = thbot(i) * (1.0 + zvir*qbot(i))
        ssq    = qsat(ts(i)) / rbot(i)
        bma1   = min( bma(i), qbot(i)/ssq )
        ssq    = min(0.030,ssq)
        qsfc   = ( 1.0 - bma1 ) * qbot(i) + bma1 * ssq
        delt   = thbot(i) - ts(i)
        delq   = qbot(i) - qsfc
        cp     = cpair*(1. + cpvir*qsfc) 
        dz     = zbot(i) - zpd(i)
        z0h1   = 0.1 * z0m(i)
!
!---------------------------------------------------------------
! First iteration to converge on Z/L and hence the fluxes
!---------------------------------------------------------------
!
! Initial guess of neutral exchange coefficients 
!
        alzm  = log( dz / z0m(i) )
        alzh  = log( dz / z0h1 )
        rdn   = xkar / alzm
        rhn   = xkar / alzh
        ren   = rhn
!
! Initial guess of ustar,tstar and qstar
!
        ustar = rdn * vmag
        tstar = rhn * delt
        qstar = ren * delq
        rd0   = rdn
        rh0   = rhn
!  iteration start
        do nn = 1, 10
!
! Compute stability and evaluate all stability functions
! Stable if (hol > 0 )
!
          hol      = xkar * gravit * dz * ( tstar / thvbot +      
     &               qstar / (1. / zvir + qbot(i)) ) / (ustar * ustar)
          hol      = sign( min(abs(hol),10._r8), hol )
          stable   = 0.5 + sign(0.5_r8 , hol)
          xsq      = max(sqrt(abs(1. - 16.*hol)) , 1._r8)
          xqq      = sqrt(xsq)
          psimh    = -5. * hol * stable + (1.-stable)*psimhu(xqq)
          psixh    = -5. * hol * stable + (1.-stable)*psixhu(xqq)
  
          z0mol    = z0m(i) / dz * hol
          z0mol    = sign( min(abs(z0mol),10._r8), z0mol )
          xsq      = max(sqrt(abs(1. - 16.*z0mol)) , 1._r8)
          xqq      = sqrt(xsq)
          psimhz0  = -5. * z0mol * stable + (1.-stable)*psimhu(xqq)
          z0hol    = z0h1 / dz * hol
          z0hol    = sign( min(abs(z0hol),10._r8), z0hol )
          xsq      = max(sqrt(abs(1. - 16.*z0hol)) , 1._r8)
          xqq      = sqrt(xsq)
          psixhz0  = -5. * z0hol * stable + (1.-stable)*psixhu(xqq)
          if ( conv ) exit
!
! Compute exchange coefficients with stability
!
          rd       = xkar / (alzm - psimh + psimhz0) 
          rh       = xkar / (alzh - psixh + psixhz0) 
          re       = rh 
!
! Update ustar, tstar, qstar 
!
          ustar    = rd * vmag 
          tstar    = rh * delt 
          qstar    = re * delq 
!
          del      = ( rd - rd0 ) / rd + ( rh - rh0 ) / rh
          if ( nn > 2 .and. del < crit ) then
            conv = .true.
!!!            print*,'took nn iteration to converge', nn
          else
!!!            print*,'nn,del,crit=',nn,del,crit
!!!            print*,'rd0,rh0,rd,rh=',rd0,rh0,rd,rh
            rd0     = rd
            rh0     = rh
          end if
        end do
!
!---------------------------------------------------------------
!    if not converge, assuming netural profile
!---------------------------------------------------------------
!
        alzm2       = log( z2m / z0m(i) )
        alzh2       = log( z2m / z0h1 )
        alzm10      = log( z10m / z0m(i) )
        alzh10      = log( z10m / z0h1 )
        if( .not. conv ) then
          ustar     = rdn * vmag
          tstar     = rhn  * delt
          qstar     = ren  * delq
          psimh2    = psimhz0
          psixh2    = psixhz0
          psimh10   = psimhz0
          psixh10   = psixhz0
        else
          z2ol      = z2m / dz * hol
          z2ol      = sign( min(abs(z2ol),10._r8), z2ol )
          xsq       = max(sqrt(abs(1. - 16.*z2ol)) , 1._r8)
          xqq       = sqrt(xsq)
          psimh2    = -5. * z2ol * stable + (1.-stable)*psimhu(xqq)
          psixh2    = -5. * z2ol * stable + (1.-stable)*psixhu(xqq)
          z10ol     = z10m / dz * hol
          z10ol     = sign( min(abs(z10ol),10._r8), z10ol )
          xsq       = max(sqrt(abs(1. - 16.*z10ol)) , 1._r8)
          xqq       = sqrt(xsq)
          psimh10   = -5. * z10ol * stable + (1.-stable)*psimhu(xqq)
          psixh10   = -5. * z10ol * stable + (1.-stable)*psixhu(xqq)
        end if
!
!---------------------------------------------------------------
! Compute surface stress components
!---------------------------------------------------------------
!
        tau     =  rbot(i) * ustar * ustar 
        taux(i) = -tau * ubot(i) / vmag 
        tauy(i) = -tau * vbot(i) / vmag 
!
!---------------------------------------------------------------
! Compute sensible and latent heat fluxes
! (Define positive latent and sensible heat as upwards into the atm)
!---------------------------------------------------------------
!
        shf(i)  = -cp * tau * tstar / ustar 
        lhf(i)  = -ltheat(i) * tau * qstar / ustar
!DEBUG
#ifdef DEBUG
        if( shf(i) > 1000.0  .or. lhf(i) > 1000.0) then
         print*,'pcols,npts,pmidm1(i)=',pcols,npts,pmidm1(i)
         print*,'ubot,vbot,tbot=',ubot(i),vbot(i),tbot(i)
         print*,'qbot,thbot,zbot=',qbot(i),thbot(i),zbot(i)
         print*,'ts,ltheat,rbot=',ts(i),ltheat,rbot(i)
         print*,'z0m,z0h,bma1=',z0m(i),z0h1,bma1
         print*,'shf(i),lhf(i),tau=',shf(i),lhf(i),tau
         print*,'tstar,qstar,ustar=',tstar,qstar,ustar
         print*,'vmag,thvbot,ssq=',vmag,thvbot,ssq
         print*,'qsfc,delt,delq=',qsfc,delt,delq
         print*,'conv,dz,zpd,=',conv,dz,zpd(i)
         print*,'rd,rh,re=',rd,rh,re
         print*,'psimh10,psixh10,hol=',psimh10,psixh10,hol
       end if
#endif
!
!---------------------------------------------------------------
! Compute  2-meter and 10 meter fileds
!---------------------------------------------------------------
!
        vmag2   = ustar * (alzm2 - psimh2 + psimhz0) / xkar
        vmag10  = ustar * (alzm10 - psimh10 + psimhz0) / xkar
        u2m(i)  = vmag2 * ubot(i) / vmag
        v2m(i)  = vmag2 * vbot(i) / vmag
        u10m(i) = vmag10 * ubot(i) / vmag
        v10m(i) = vmag10 * vbot(i) / vmag
        th2m    = ts(i) + tstar * ( alzh2 - psixh2 + psixhz0 ) / xkar
        th10m   = ts(i) + tstar * ( alzh10 - psixh10 + psixhz0 ) / xkar
        q2m(i)  = qsfc + qstar * ( alzh2 - psixh2 + psixhz0 ) / xkar
        q10m(i) = qsfc + qstar * ( alzh10 - psixh10 + psixhz0 ) / xkar
!---------------------------------------------------------------
! Compute temperature at 10 and 2 meter
! dT/T   = dth/th  + rd/cp * dP/P
! t10m   = th10m * tbot / thbot + g * (zbot-z10) / cp / (1.0+zvir* qbot)
! t2m    = th2m * tbot / thbot + g * (zbot-z2) / cp / (1.0+zvir* qbot)
!---------------------------------------------------------------
!
        t2m(i)  = th2m * tbot(i) / thbot(i) 
     &          + gravit * ( zbot(i) - z2m ) / cpair 
     &          / (1.0 + zvir * qbot(i))
        t10m(i) = th10m * tbot(i) / thbot(i) 
     &          + gravit * ( zbot(i) - z10m ) / cpair 
     &          / (1.0 + zvir * qbot(i))
      end do
!
      return
      end subroutine diaglnd

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  diagoce -- diagnostics of surface fluxes and 2-meter and 10-meter 
!                        fields over ocean
!
! !INTERFACE:
!

      subroutine diagoce(pcols    ,
     &           indx    ,npts    ,pmidm1  ,ubot    ,vbot    ,
     &           tbot    ,qbot    ,thbot   ,zbot    ,ts      ,
     &           ltheat  ,shf     ,lhf     ,taux    ,tauy    ,
     &           t2m     ,q2m     ,u2m     ,v2m     ,        
     &           t10m    ,q10m    ,u10m    ,v10m    )

!
! !USES:
!
      use precision
      implicit none
#include <comcon.h>
#include <parpbl.h>
!
! !INPUT PARAMETERS:
!
      integer , intent(in) :: pcols          ! number of column 
      integer , intent(in) :: indx(pcols)    ! column index array (ocean)
      integer , intent(in) :: npts           ! Number of ocean points
      real(r8), intent(in) :: pmidm1(pcols)  ! Bottom level pressure
      real(r8), intent(in) :: ubot(pcols)    ! Bottom level u wind
      real(r8), intent(in) :: vbot(pcols)    ! Bottom level v wind
      real(r8), intent(in) :: tbot(pcols)    ! Bottom level temperature
      real(r8), intent(in) :: qbot(pcols)    ! Bottom level specific humidity
      real(r8), intent(in) :: thbot(pcols)   ! Bottom level potential temperature
      real(r8), intent(in) :: zbot(pcols)    ! Bottom level height above surface
      real(r8), intent(in) :: ts(pcols)      ! Surface temperature
      real(r8), intent(in) :: ltheat(pcols)  ! Latent heat for given srf conditions

!
! !OUTPUT PARAMETERS:
!
      real(r8), intent(out) :: shf(pcols)    ! sensible heat flux (W/m2)
      real(r8), intent(out) :: lhf(pcols)    ! latent heat flux (W/m2)
      real(r8), intent(out) :: taux(pcols)   ! X surface stress (N/m2)
      real(r8), intent(out) :: tauy(pcols)   ! Y surface stress (N/m2)
      real(r8), intent(out) :: t2m(pcols)    ! 2-meter temperature (K)
      real(r8), intent(out) :: q2m(pcols)    ! 2-meter specific humidity (kg/kg)
      real(r8), intent(out) :: u2m(pcols)    ! 2-meter u-wind (m/s)
      real(r8), intent(out) :: v2m(pcols)    ! 2-meter v-wind (m/s)
      real(r8), intent(out) :: t10m(pcols)   ! 10-meter temperature (K)
      real(r8), intent(out) :: q10m(pcols)   ! 10-meter specific humidity (kg/kg)
      real(r8), intent(out) :: u10m(pcols)   ! 10-meter u-wind (m/s)
      real(r8), intent(out) :: v10m(pcols)   ! 10-meter v-wind (m/s)

! !DESCRIPTION: Given the near surface analysis fields to compute surface fluxes 
!               and 2-meter and 10-meter fields over ocean based on similarity theory.
!
! Method: Use the same method as CCM
! Assume:
!   1) Neutral 10m drag coeff: 
!         cdn = .0027/U10N + .000142 + .0000764 U10N
!   2) Neutral 10m stanton number: 
!         ctn = .0327 sqrt(cdn), unstable
!         ctn = .0180 sqrt(cdn), stable
!   3) Neutral 10m dalton number:  
!         cen = .0346 sqrt(cdn)
!   4) The saturation humidity of air at T(K): 
!         qsat(T)  (kg/m^3)
! Note:
!   1) here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!   2) wind speeds should all be above a minimum speed (umin)
! 
!-----------------------------------------------------------------------
! !REVISION HISTORY:
!
! 10Jan2001  Chern     initial code
!
!EOP

!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
      real(r8), parameter :: z2m  = 2.0 !  2-m height
      real(r8), parameter :: z10m = 10.0 !  10-m height
      integer i,ii            ! column indices
      real(r8) ssq            ! Surface saturation specific humidity
      real(r8) ustar          ! ustar
      real(r8) tstar          ! tstar
      real(r8) qstar          ! qstar
      real(r8) u10n           ! neutral 10 m wind speed over ocean
      real(r8) vmag           ! Bottom level wind magnitude
      real(r8) vmag2          ! 2-meter wind magnitude
      real(r8) vmag10         ! 10-meter wind magnitude
      real(r8) thvbot         ! Bottom lev virtual potential temp
      real(r8) th2m           ! 2-meter potential temp
      real(r8) th10m          ! 10-meter potential temp
      real(r8) delt           ! potential T difference (K)
      real(r8) delq           ! specific humidity difference (kg/kg)
      real(r8) rdn            ! sqrt of neutral exchange coeff (momentum)
      real(r8) rhn            ! sqrt of neutral exchange coeff (heat)
      real(r8) ren            ! sqrt of neutral exchange coeff (tracers)          
      real(r8) rd             ! sqrt of exchange coeff (momentum)
      real(r8) rh             ! sqrt of exchange coeff (heat)
      real(r8) re             ! sqrt of exchange coeff (tracers)
      real(r8) hol            ! Ref hgt (10m) / monin-obukhov length
      real(r8) hol2           ! 2 meter / monin-obukhov length
      real(r8) hol10          ! 10 meter / monin-obukhov length
      real(r8) xsq            ! Temporary variable
      real(r8) xqq            ! Temporary variable
      real(r8) alz            ! ln(zbot/z10)
      real(r8) alz10          ! ln(z10m/zbot)
      real(r8) alz2           ! ln(z2m/zbot)
      real(r8) cp             ! Specific heat of moist air
      real(r8) tau            ! Reference height stress
      real(r8) psimh          ! Stability funct at bottom level(momentum)
      real(r8) psimh2         ! Stability funct at 2 meter (momentum)
      real(r8) psimh10        ! Stability funct at 10 meter (momentum)
      real(r8) psixh          ! Stability funct at bottom level (heat & tracers) 
      real(r8) psixh2         ! Stability funct at 2 meter (heat & tracers) 
      real(r8) psixh10        ! Stability funct at 10 meter (heat & tracers) 
      real(r8) stable         ! Stability factor
      real(r8) rbot(pcols)    ! Density at bottom model level
!-----------------------------------------------------------------------
    
!--------------------------Statement functions--------------------------
      real(r8) psimhu         ! Unstable part of psimh
      real(r8) psixhu         ! Unstable part of psixh
      real(r8) qsat           ! Saturation specific humidty of air
      real(r8) cdn            ! Neutral drag coeff at bottom model level
      real(r8) xd             ! Dummy argument
      real(r8) Tk             ! Temperature (K)
      real(r8) Umps           ! Wind velocity (m/sec)
!
      qsat(Tk)   = 640380. / exp(5107.4/Tk)
      cdn(Umps)  = 0.0027 / Umps + .000142 + .0000764 * Umps
      psimhu(xd) = log((1.+xd*(2.+xd))*(1.+xd*xd)/8.) - 2.*atan(xd)
     &           + 1.571
      psixhu(xd) = 2. * log((1. + xd*xd)/2.)
!-----------------------------------------------------------------------
!
! Loop over ocean points
!
      do ii=1,npts
        i = indx(ii)
!
!---------------------------------------------------------------
! Set up necessary variables
!---------------------------------------------------------------
!
        rbot(i)= pmidm1(i) / (rair*tbot(i))
        vmag   = max(umin, sqrt(ubot(i)**2+vbot(i)**2))
        thvbot = thbot(i) * (1.0 + zvir*qbot(i))
        ssq    = 0.98 * qsat(ts(i)) / rbot(i)
        delt   = thbot(i) - ts(i)
        delq   = qbot(i) - ssq  
        alz    = log(zbot(i)/zref) 
        cp     = cpair*(1. + cpvir*ssq) 
!
!---------------------------------------------------------------
! First iteration to converge on Z/L and hence the fluxes
!---------------------------------------------------------------
!
! Initial guess for roots of neutral exchange coefficients, 
! assume z/L=0. and u10n is approximated by vmag.
! Stable if (thbot > ts ).
!
        stable = 0.5 + sign(0.5_r8 , delt)
        rdn  = sqrt(cdn(vmag))
        rhn  = (1.-stable) * 0.0327 + stable * 0.018 
        ren  = 0.0346 
!
! Initial guess of ustar,tstar and qstar
!
        ustar = rdn*vmag
        tstar = rhn*delt
        qstar = ren*delq
!
! Compute stability and evaluate all stability functions
! Stable if (thbot > ts or hol > 0 )
!
        hol = xkar * gravit * zbot(i) * ( tstar / thvbot + 
     &        qstar / (1./zvir+qbot(i)) ) / ustar**2
        hol = sign( min(abs(hol),10._r8), hol )
        stable = 0.5 + sign(0.5_r8 , hol)
        xsq   = max(sqrt(abs(1. - 16.*hol)) , 1._r8)
        xqq   = sqrt(xsq)
        psimh = -5. * hol * stable + (1.-stable)*psimhu(xqq)
        psixh = -5. * hol * stable + (1.-stable)*psixhu(xqq)
!
! Shift 10m neutral wind speed using old rdn coefficient
!
        rd   = rdn / (1.+rdn/xkar*(alz-psimh))
        u10n = vmag * rd / rdn
!
! Update the neutral transfer coefficients at 10m and neutral stability
!
        rdn = sqrt(cdn(u10n))
        ren = 0.0346
        rhn = (1.-stable) * 0.0327 + stable * 0.018 
!
! Shift all coeffs to measurement height and stability
!
        rd = rdn / (1.+rdn/xkar*(alz-psimh)) 
        rh = rhn / (1.+rhn/xkar*(alz-psixh)) 
        re = ren / (1.+ren/xkar*(alz-psixh))
!
! Update ustar, tstar, qstar using updated, shifted coeffs 
!
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 
!
!---------------------------------------------------------------
! Second iteration to converge on Z/L and hence the fluxes
!---------------------------------------------------------------
!
! Recompute stability & evaluate all stability functions  
! Stable if (thbot > ts or hol > 0 )
! 
        hol = xkar * gravit * zbot(i) * ( tstar / thvbot 
     &      + qstar / (1./zvir+qbot(i)) ) / ustar**2
        hol = sign( min(abs(hol),10._r8), hol )
        stable = 0.5 + sign(0.5_r8 , hol)
        xsq   = max(sqrt(abs(1. - 16.*hol)) , 1._r8)
        xqq   = sqrt(xsq)
        psimh = -5. * hol * stable + (1.-stable)*psimhu(xqq)
        psixh = -5. * hol * stable + (1.-stable)*psixhu(xqq)
!
! Shift 10m neutral wind speed using old rdn coefficient
!
        rd   = rdn / (1.+rdn/xkar*(alz-psimh))
        u10n = vmag * rd / rdn
!
! Update the neutral transfer coefficients at 10m and neutral stability
!
        rdn = sqrt(cdn(u10n))
        ren = 0.0346
        rhn = (1.-stable) * 0.0327 + stable * 0.018 
!
! Shift all coeffs to measurement height and stability

        rd = rdn / (1.+rdn/xkar*(alz-psimh)) 
        rh = rhn / (1.+rhn/xkar*(alz-psixh)) 
        re = ren / (1.+ren/xkar*(alz-psixh))
!
!---------------------------------------------------------------
! Compute the fluxes
!---------------------------------------------------------------
!
! Update ustar, tstar, qstar using updated, shifted coeffs 
!
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 
!
! Compute surface stress components
!
        tau     =  rbot(i) * ustar * ustar 
        taux(i) = -tau * ubot(i) / vmag 
        tauy(i) = -tau * vbot(i) / vmag 
!
! Compute heat flux components at current surface temperature
! (Define positive latent and sensible heat as upwards into the atm)
!
        shf(i) = -cp * tau * tstar / ustar 
        lhf(i) = -ltheat(i) * tau * qstar / ustar
!
!---------------------------------------------------------------
! Compute 2-meter and 10 meter fields
!---------------------------------------------------------------
!
        hol     = xkar *gravit *zbot(i) * ( tstar/thbot(i) + 
     &            qstar/(1./zvir+qbot(i)) ) / (ustar * ustar)
        hol     = sign( min(abs(hol),10.), hol )
        hol10   = hol * z10m / zbot(i)
        hol2    = hol * z2m / zbot(i)
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
        alz2    = log( z2m / zbot(i) )
        alz10   = log( z10m / zbot(i) )
        vmag2   = vmag + ustar * ( alz2 - psimh2 + psimh ) / xkar
        vmag10  = vmag + ustar * ( alz10 - psimh10 + psimh ) / xkar
        u10m(i) = vmag10 * ubot(i) / vmag
        v10m(i) = vmag10 * vbot(i) / vmag
        u2m(i)  = vmag2 * ubot(i) / vmag
        v2m(i)  = vmag2 * vbot(i) / vmag
!---------------------------------------------------------------
! Compute potential temperature and mixing ratio at 10 and 2 meter
! th10m  = thbot + tstar/xkar * ( alz10-psixh10+psixh )
! q10m   = qbot + qstar/xkar * ( alz10-psixh10+psixh )
! th2m   = thbot + tstar/xkar * ( alz2-psixh2+psixh )
! q2m    = qbot + qstar/xkar * ( alz2-psixh2+psixh )
!---------------------------------------------------------------
!
        th2m    = thbot(i) + tstar/xkar * ( alz2-psixh2+psixh )
        q2m(i)  = qbot(i) + qstar/xkar * ( alz2-psixh2+psixh )
        th10m   = thbot(i) + tstar/xkar * ( alz10-psixh10+psixh )
        q10m(i) = qbot(i) + qstar/xkar * ( alz10-psixh10+psixh )
!
!---------------------------------------------------------------
! Compute temperature at 10 and 2 meter
! dT/T   = dth/th  + rd/cp * dP/P
! t10m   = th10m * tbot / thbot + g * (zbot-z10) / cp / (1.0+zvir* qbot)
! t2m    = th2m * tbot / thbot + g * (zbot-z2) / cp / (1.0+zvir* qbot)
!---------------------------------------------------------------
!
        t2m(i)  = th2m * tbot(i) / thbot(i)  + gravit 
     &          * ( zbot(i) - z2m ) / cpair / (1.0 + zvir * qbot(i))
        t10m(i) = th10m * tbot(i) / thbot(i) + gravit 
     &          * ( zbot(i) - z10m ) / cpair / (1.0 + zvir * qbot(i))

      end do
!
      return
      end subroutine diagoce

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  diagice -- diagnostics of surface fluxes and 2-meter and 10-meter
!                        fields over sea ice
!
! !INTERFACE:
!

      subroutine diagice(pcols    ,
     &           indx    ,npts    ,pmidm1  ,ubot    ,vbot    ,
     &           tbot    ,qbot    ,thbot   ,zbot    ,ts      ,
     &           ltheat  ,shf     ,lhf     ,taux    ,tauy    ,
     &           t2m     ,q2m     ,u2m     ,v2m     ,        
     &           t10m    ,q10m    ,u10m    ,v10m    )

! !USES:
!
      use precision
      implicit none
#include <comcon.h>
#include <parpbl.h>
!
! !INPUT PARAMETERS:
!
      integer , intent(in) :: pcols          ! number of column 
      integer , intent(in) :: indx(pcols)    ! column index array (seaice)
      integer , intent(in) :: npts           ! Number of seaice points
      real(r8), intent(in) :: pmidm1(pcols)  ! Bottom level pressure
      real(r8), intent(in) :: ubot(pcols)    ! Bottom level u wind
      real(r8), intent(in) :: vbot(pcols)    ! Bottom level v wind
      real(r8), intent(in) :: tbot(pcols)    ! Bottom level temperature
      real(r8), intent(in) :: qbot(pcols)    ! Bottom level specific humidity
      real(r8), intent(in) :: thbot(pcols)   ! Bottom level potential temperature
      real(r8), intent(in) :: zbot(pcols)    ! Bottom level height above surface
      real(r8), intent(in) :: ts(pcols)      ! Surface temperature
      real(r8), intent(in) :: ltheat(pcols)  ! Latent heat for given srf conditions
!
! !OUTPUT PARAMETERS:
!
      real(r8), intent(out)  :: shf(pcols)   ! sensible heat flux (W/m2)
      real(r8), intent(out)  :: lhf(pcols)   ! latent heat flux (W/m2)
      real(r8), intent(out)  :: taux(pcols)  ! X surface stress (N/m2)
      real(r8), intent(out)  :: tauy(pcols)  ! Y surface stress (N/m2)
      real(r8), intent(out)  :: t2m(pcols)   ! 2-meter temperature (K)
      real(r8), intent(out)  :: q2m(pcols)   ! 2-meter specific humidity (kg/kg)
      real(r8), intent(out)  :: u2m(pcols)   ! 2-meter u-wind (m/s)
      real(r8), intent(out)  :: v2m(pcols)   ! 2-meter v-wind (m/s)
      real(r8), intent(out)  :: t10m(pcols)  ! 10-meter temperature (K)
      real(r8), intent(out)  :: q10m(pcols)  ! 10-meter specific humidity (kg/kg)
      real(r8), intent(out)  :: u10m(pcols)  ! 10-meter u-wind (m/s)
      real(r8), intent(out)  :: v10m(pcols)  ! 10-meter v-wind (m/s)

! !DESCRIPTION: Given the near surface analysis fields to compute surface fluxes
!               and 2-meter and 10-meter fields over ocean based on similarity theory.
!
! Method: Use the same method as CCM
!-----------------------------------------------------------------------
! !REVISION HISTORY:
!
! 10Jan2001  Chern     initial code
!
!EOP

!---------------------------Local variables-----------------------------
      real(r8), parameter :: z2m  = 2.0 !  2-m height
      real(r8), parameter :: z10m = 10.0 !  10-m height
      integer i,ii            ! Column indices
      real(r8) ssq            ! Surface saturation specific humidity
      real(r8) ustar          ! ustar
      real(r8) tstar          ! tstar
      real(r8) qstar          ! qstar
      real(r8) vmag           ! Surface wind magnitude
      real(r8) vmag2          ! 2-meter wind magnitude
      real(r8) vmag10         ! 10-meter wind magnitude
      real(r8) thvbot         ! Bottom level virtual potential temperature
      real(r8) th2m           ! 2-meter potential temp
      real(r8) th10m          ! 10-meter potential temp
      real(r8) delt           ! potential T difference (k)
      real(r8) delq           ! specific humidity difference (kg/kg)
      real(r8) rdn            ! sqrt of neutral exchange coeff (momentum)
      real(r8) rhn            ! sqrt of neutral exchg coeff (heat & tracers)
      real(r8) ren            ! sqrt of neutral exchange coeff (water)
      real(r8) rd             ! sqrt of exchange coefficient (momentum)
      real(r8) rh             ! sqrt of exchange coefficient(heat & tracers)
      real(r8) re             ! sqrt of exchange coefficient (water) 
      real(r8) hol            ! bottom level / monin-obukhov length
      real(r8) hol2           ! 2 meter / monin-obukhov length
      real(r8) hol10          ! 10 meter / monin-obukhov length
      real(r8) xsq            ! Temporary variable
      real(r8) xqq            ! Temporary variable
      real(r8) alz            ! ln(zbot/10)
      real(r8) alz10          ! ln(z10m/zbot)
      real(r8) alz2           ! ln(z2m/zbot)
      real(r8) cp             ! Specific heat of moist air
      real(r8) tau            ! Reference height stress
      real(r8) psimh          ! Stability function at ref lev (momentum)
      real(r8) psimh2         ! Stability funct at 2 meter (momentum)
      real(r8) psimh10        ! Stability funct at 10 meter (momentum)
      real(r8) psixh          ! Stability function at ref lev (heat & tracers) 
      real(r8) psixh2         ! Stability funct at 2 meter (heat & tracers)
      real(r8) psixh10        ! Stability funct at 10 meter (heat & tracers)
      real(r8) stable         ! Stability factor
      real(r8) rbot(pcols)    ! Density at bottom model level
!-----------------------------------------------------------------------

!------------------------------Functions--------------------------------
      real(r8) psimhu         ! Unstable part of psimh
      real(r8) psixhu         ! Unstable part of psixh
      real(r8) qsat           ! Saturation humidty of air
      real(r8) xd             ! Dummy argument
      real(r8) Tk             ! Temperature (K)
!
      qsat(Tk)    = 640380. / exp(5107.4/Tk)
      psimhu(xd)  = log((1.+xd*(2.+xd))*(1.+xd*xd)/8.) - 2.*atan(xd) 
     &            + 1.571
      psixhu(xd)  = 2. * log((1. + xd*xd)/2.)
!-----------------------------------------------------------------------
!
! Loop over ice points
!
      do ii=1,npts
        i = indx(ii)
!
!-----------------------------------------------------------------------
! Determine some necessary variables
!-----------------------------------------------------------------------
!
        rbot(i)= pmidm1(i)/(rair*tbot(i))
        vmag   = max(umin, sqrt(ubot(i)**2+vbot(i)**2))
        thvbot = thbot(i) * (1.0 + zvir*qbot(i))
        ssq    = qsat(ts(i)) / rbot(i)
        delt   = thbot(i) - ts(i)
        delq   = qbot(i) - ssq 
        alz    = log(zbot(i)/zref) 
        cp     = cpair*(1. + cpvir*ssq) 
!
!---------------------------------------------------------------
! First iteration to converge on Z/L and hence the fluxes
!---------------------------------------------------------------
!
! Determine roots of neutral exchange coefficients
!
        rdn = xkar/log(zref/zzsice)
        rhn = rdn
        ren = rdn
!
! Determine initial guess of ustar,tstar and qstar
!
        ustar = rdn*vmag
        tstar = rhn*delt
        qstar = ren*delq
!
! Compute stability and evaluate all stability functions
! Stable if (thbot > ts or hol > 0 )
!
        hol = xkar * gravit * zbot(i) * ( tstar / thvbot 
     &      + qstar / (1./zvir+qbot(i)) ) / ustar**2
        hol = sign( min(abs(hol),10._r8), hol )
        stable = 0.5 + sign(0.5_r8 , hol)
        xsq   = max(sqrt(abs(1. - 16.*hol)) , 1._r8)
        xqq   = sqrt(xsq)
        psimh = -5. * hol * stable + (1.-stable)*psimhu(xqq)
        psixh = -5. * hol * stable + (1.-stable)*psixhu(xqq)
!
! Shift all coeffs to measurement height and stability
!
        rd = rdn / (1.+rdn/xkar*(alz-psimh)) 
        rh = rhn / (1.+rhn/xkar*(alz-psixh)) 
        re = ren / (1.+ren/xkar*(alz-psixh))
!
! Update ustar, tstar, qstar using updated, shifted coeffs 
!
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 
!
!---------------------------------------------------------------
! Second iteration to converge on Z/L and hence the fluxes
!---------------------------------------------------------------
!
! Recompute stability & evaluate all stability functions  
! Stable if (thbot > ts or hol > 0 )
! 
        hol = xkar * gravit * zbot(i) * (tstar/thvbot 
     &      + qstar/(1./zvir+qbot(i))) / ustar**2
        hol = sign( min(abs(hol),10._r8), hol )
        stable = 0.5 + sign(0.5_r8 , hol)
        xsq   = max(sqrt(abs(1. - 16.*hol)) , 1._r8)
        xqq   = sqrt(xsq)
        psimh = -5. * hol * stable + (1.-stable)*psimhu(xqq)
        psixh = -5. * hol * stable + (1.-stable)*psixhu(xqq)
!
! Shift all coeffs to measurement height and stability
!
        rd = rdn / (1.+rdn/xkar*(alz-psimh)) 
        rh = rhn / (1.+rhn/xkar*(alz-psixh)) 
        re = ren / (1.+ren/xkar*(alz-psixh)) 
!
!---------------------------------------------------------------
! Compute the fluxes
!---------------------------------------------------------------
!
! Update ustar, tstar, qstar using updated, shifted coeffs 
!
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 
!
! Compute surface stress components
!
        tau     =  rbot(i) * ustar * ustar 
        taux(i) = -tau * ubot(i) / vmag 
        tauy(i) = -tau * vbot(i) / vmag 
!
! Compute heat flux components at current surface temperature
! (Define positive latent and sensible heat as upwards into the atm)
!
        shf(i) = -cp * tau * tstar / ustar 
        lhf(i) = -ltheat(i) * tau * qstar / ustar
!
!---------------------------------------------------------------
! Compute 2-meter and 10 meter fields
!---------------------------------------------------------------
!
        hol     = xkar *gravit *zbot(i) * ( tstar/thbot(i) + 
     &            qstar/(1./zvir+qbot(i)) ) / (ustar * ustar)
        hol     = sign( min(abs(hol),10.), hol )
        hol10   = hol * z10m / zbot(i)
        hol2    = hol * z2m / zbot(i)
        stable  = 0.5 + sign(0.5_r8 , hol)

        xsq     = max(sqrt(abs(1. - 16.*hol)) , 1._r8)
        xqq     = sqrt(xsq)
        psimh   = -5. * hol * stable + (1.-stable)*psimhu(xqq)
        psixh   = -5. * hol * stable + (1.-stable)*psixhu(xqq)

        xsq     = max(sqrt(abs(1. - 16.*hol10)) , 1._r8)
        xqq     = sqrt(xsq)
        psimh10 = -5. * hol10 * stable + (1.-stable)*psimhu(xqq)
        psixh10 = -5. * hol10 * stable + (1.-stable)*psixhu(xqq)

        xsq     = max(sqrt(abs(1. - 16.*hol2)) , 1._r8)
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
        alz2    = log( z2m / zbot(i) )
        alz10   = log( z10m / zbot(i) )
        vmag2   = vmag + ustar * ( alz2 - psimh2 + psimh ) / xkar
        vmag10  = vmag + ustar * ( alz10 - psimh10 + psimh ) / xkar
        u10m(i) = vmag10 * ubot(i) / vmag
        v10m(i) = vmag10 * vbot(i) / vmag
        u2m(i)  = vmag2 * ubot(i) / vmag
        v2m(i)  = vmag2 * vbot(i) / vmag
!---------------------------------------------------------------
! Compute potential temperature and mixing ratio at 10 and 2 meter
! th10m  = thbot + tstar/xkar * ( alz10-psixh10+psixh )
! q10m   = qbot + qstar/xkar * ( alz10-psixh10+psixh )
! th2m   = thbot + tstar/xkar * ( alz2-psixh2+psixh )
! q2m    = qbot + qstar/xkar * ( alz2-psixh2+psixh )
!---------------------------------------------------------------
!
        th2m    = thbot(i) + tstar/xkar * ( alz2-psixh2+psixh )
        q2m(i)  = qbot(i) + qstar/xkar * ( alz2-psixh2+psixh )
        th10m   = thbot(i) + tstar/xkar * ( alz10-psixh10+psixh )
        q10m(i) = qbot(i) + qstar/xkar * ( alz10-psixh10+psixh )
!
!---------------------------------------------------------------
! Compute temperature at 10 and 2 meter
! dT/T   = dth/th  + rd/cp * dP/P
! t10m   = th10m * tbot / thbot + g * (zbot-z10) / cp / (1.0+zvir* qbot)
! t2m    = th2m * tbot / thbot + g * (zbot-z2) / cp / (1.0+zvir* qbot)
!---------------------------------------------------------------
!
        t2m(i)  = th2m * tbot(i) / thbot(i)  + gravit * 
     &          ( zbot(i) - z2m ) / cpair / (1.0 + zvir * qbot(i))
        t10m(i) = th10m * tbot(i) / thbot(i) + gravit * 
     &          ( zbot(i) - z10m ) / cpair / (1.0 + zvir * qbot(i))
      end do
!
      return
      end subroutine diagice

      end module m_diagsfc

