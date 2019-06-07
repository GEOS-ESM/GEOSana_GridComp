#include <params.h>
      subroutine tphysaa(lat ,nstep   ,tdt     ,
     &               gravit  ,rair    ,zvir    ,
     &               cpair   ,latvap  ,rga     ,tmelt   ,
     &               phis    ,ts      ,pblht   ,tpert   ,qpert   ,
     &               pintm1  ,pilnm1  ,pdelm1  ,pmidm1  ,pmlnm1  ,
     &               tm1     ,qm1     ,tvm     ,
     &               um1     ,vm1     ,psomc   ,
     &               taux    ,tauy    ,shflx   ,cflx    ,
     &               tstar   ,qstar            )
!-----------------------------------------------------------------------
!
! Purpose:
! This program is used to correct the excessive Precipitation after FVDAS.
!
! Method:
! Call moist physics and PBL subroutines to adjust q and t
!
! Author: Jiundar Chern  April 2001
!
!-----------------------------------------------------------------------
!
       implicit none

!--------- Common Blocks -----------------------------------------------
#include <pmgrid.h>       ! use only: plond, plev, plevp, pcnst, pnats
!
! Input Arguments for moist physics
!
      integer, intent(in) :: lat                   ! Latitude index (S->N)
      integer, intent(in) :: nstep                 ! Current time step
      real,    intent(in) :: tdt                   ! 2 delta t (model time increment)
      real,    intent(in) :: gravit                ! Gravitational acceleration
      real,    intent(in) :: rair                  ! Gas constant for dry air
      real,    intent(in) :: zvir                  ! rh2o/rair - 1
      real,    intent(in) :: cpair                 ! Specific heat of dry air
      real,    intent(in) :: latvap                ! Latent heat of vaporization
      real,    intent(in) :: rga                   ! Reciprocal of gravitational acceleration
      real,    intent(in) :: tmelt                 ! Freezing point of water
      real,    intent(in) :: phis(plond)           ! Surface geopotential
      real,    intent(in) :: ts(plond)             ! Surface temperature
      real,    intent(in) :: pblht(plond)          ! Planetary boundary layer height
      real,    intent(in) :: pintm1(plond,plevp)   ! Interface pressures
      real,    intent(in) :: pilnm1(plond,plevp)   ! log(pintm1)
      real,    intent(in) :: pdelm1(plond,plev)    ! Layer pressure depth
      real,    intent(in) :: pmidm1(plond,plev)    ! Midpoint pressures
      real,    intent(in) :: pmlnm1(plond,plev)    ! log(pmidm1)
      real,    intent(in) :: tm1(plond,plev)       ! Temperature
      real,    intent(in) :: qm1(plond,plev,pcnst+pnats) ! Specific humidity
      real,    intent(in) :: tvm(plond,plev)       ! Virtual temperature
!
! Input Arguments for PBL
!
      real,    intent(in) :: taux(plond)            ! X surface stress (zonal)
      real,    intent(in) :: tauy(plond)            ! Y surface stress (zonal)
      real,    intent(in) :: shflx(plond)           ! Surface sensible heat flux (J/m2/s)
      real,    intent(in) :: cflx(plond,pcnst)      ! Surface constituent flux (kg/m^2/s)
      real,    intent(in) :: um1(plond,plev)        ! U horizontal wind component (m/s)
      real,    intent(in) :: vm1(plond,plev)        ! V horizontal wind component (m/s)
      real,    intent(in) :: psomc(plond,plev)      ! pk(i,plevp) / pkz(i,k)
!
! Output/Input Arguments
!
      real,    intent(inout) :: tpert(plond)          ! Thermal temperature excess
      real,    intent(inout) :: qpert(plond,pcnst)    ! Thermal humidity & constituent excess

!
! Output arguments
!
      real,    intent(out) :: tstar(plond,plev)    ! adjusted temperature
      real,    intent(out) :: qstar(plond,plev,pcnst+pnats)      ! adjusted q fields
!
! Local Variables
!

      real zm(plond,plev)           ! Height of midpoints (above surface)
      real zi(plond,plevp)          ! Height above sfc in interface

      real precl(plond)             ! Large-scale precipitation rate
      real precc(plond)             ! Convective-scale preciptn rate
      real cnt(plond)               ! Top level of convective activity
      real cnb(plond)               ! Bottom level of convective activity
      real cmfmc(plond,plev)        ! Moist convection cloud mass flux Total
      real cmfdt(plond,plev)        ! dT/dt due to moist convection
      real cmfdq(plond,plev)        ! dq/dt due to moist convection
      real cmfdqr(plond,plev)       ! dq/dt due to moist convective rainout 
      real qc(plond,plev)           ! dq/dt due to rainout terms
      real zmdt(plond,plev)         ! Zhang convective temperature tendency
      real zmdq(plond,plev)         ! Zhang convective moisture tendency

      real precc2(plond)            ! Convective-scale preciptn rate
      real cnt2(plond)              ! Top level of convective activity
      real cnb2(plond)              ! Bottom level of convective activity
      real cmfmc2(plond,plev)       ! Moist convection cloud mass flux Hack
      real cmfdqr2(plond,plev)      ! dq/dt due to moist convective rainout 
      real cmfsl2(plond,plev)       ! Moist convection lw stat energy flux
      real cmflq2(plond,plev)       ! Moist convection total water flux
      real qc2(plond,plev)          ! dq/dt due to rainout terms
      real tpert2(plond)            ! Perturbation T
      real qpert2(plond,pcnst)      ! Perturbation q
      real dlf(plond,plev)          ! Detraining cld H20 from convection
      real pflx(plond,plev)         ! Conv rain flux thru out btm of lev
C JDC ADDED
      real htlcl(plond)         ! height above surface at LCL level
      real htmmse(plond)        ! height above surface at maximum 
                                ! moist static energy level
      real cape(plond)          ! maximum CAPE at a column
      real cape2d(plond,plev)   ! Convective Available Potential Energy (CAPE)
      real evp(plond,plev)      ! Evaporation for downdraft - Z&M scheme
      real cme(plond,plev)      ! Condensation - evaporation - Z&M scheme
      real mup(plond,plev)      ! Updraft mass flux - Z&M scheme
      real mdn(plond,plev)      ! Downdraft mass flux - Z&M scheme
      real qlup(plond,plev)     ! Updraft liquid water - Z&M scheme
      real dup(plond,plev)      ! Updraft detrainment - Z&M scheme
      real eup(plond,plev)      ! Updraft entrainment - Z&M scheme
      real edn(plond,plev)      ! Downdraft entrainment - Z&M scheme
      real zmeps(plond,plev)    ! fractional entrainment rate - Z&M scheme
      real mbup(plond)          ! Cloud base mass flux - Z&M scheme
      real cmfetr(plond,plev)   ! entrainment - Hack scheme
      real cmfdtr(plond,plev)   ! detrainment - Hack scheme
      real evapl(plond,plev)    ! Large scale rain evaporation
      real relhum(plond,plev)   ! relative humidity after cloud physics
      real rhclr(plond,plev)    ! relative humidity in clear region
      real convcld(plond,plev)  ! convective cloud amount
      real stratcld(plond,plev) ! stratiform cloud amount
      real cldtmp(plond,plev)   ! Temporary array
      real dtpblcg(plond,plev)    ! T tendency due to counter gradient term
      real dqpblcg(plond,plev)    ! q tendency due to counter gradient term

      real rpdel(plond,plev)        ! 1./(pint(k+1)-pint(k))
      real fac
      real rtdt
      real rrg
      integer  i,k,m                ! Lon, lev, constituent indices
!
!--- Rain re-evaporation section -----------------
      integer  k1
      real sumflx(plond)          ! flux integral
      real est(plond,plev)        ! Saturation vapor pressure
      real qsat(plond,plev)
      real pflxtmp                ! Conv rain flux thru out btm of lev
      real ke
      real rlovrcp
      real dpovrg
      real envevap
      real tcrit

!
! variables for PBL      
!
      real told(plond,plev)         ! temperature before vertical diffusion
      real qold(plond,plev,pcnst)   ! Specific humidity before diffusion
      real tvm1(plond,plev)         ! virtual temperature before diffusion
      real rpdeli(plond,plev)       ! 1.0 / pdeli (thickness between midpoints)
      real thm(plond,plev)          ! Potential temperature
      real ustar(plond)             ! surface friction velocity
      real kvh(plond,plevp)         ! diffusion coefficient for heat
      real kvm(plond,plevp)         ! diffusion coefficient for momentum
      real cgs(plond,plevp)         ! counter-gradient star (cg/flux)
      real up1(plond,plev)          ! u-wind after vertical diffusion
      real vp1(plond,plev)          ! v-wind after vertical diffusion
      real thp(plond,plev)          ! potential temperature after vdiff

      logical deep_conv 
      data    deep_conv / .true. /

      rtdt = 1./tdt

      do m = 1, pcnst + pnats
        do k = 1, plev
          do i = 1, plon
            qstar(i,k,m) = max( 0., qm1(i,k,m) )
          end do
        end do
      end do

      do k = 1, plev
        do i = 1, plon
          tstar(i,k)     = tm1(i,k)
        end do
      end do
!
! Zero out precip and convective fields before accumulating terms
!
      call resetr(precl  ,plond      ,0.)
      call resetr(precc  ,plond      ,0.)
      call resetr(qc     ,plond*plev ,0.)
      call resetr(cmfdt  ,plond*plev ,0.)
      call resetr(cmfdq  ,plond*plev ,0.)
      call resetr(zmdt   ,plond*plev ,0.)
      call resetr(zmdq   ,plond*plev ,0.)
      call resetr(cmfdqr ,plond*plev ,0.)
      call resetr(cmfmc  ,plond*plev ,0.)
      call resetr(cldtmp ,plond*plev ,0.)
      call resetr(cme    ,plond*plev ,0.)
      call resetr(dlf    ,plond*plev ,0.)
      call resetr(pflx   ,plond*plev ,0.)
      call resetr(evp    ,plond*plev ,0.)
      call resetr(mdn    ,plond*plev ,0.)
      call resetr(mup    ,plond*plev ,0.)
      call resetr(qlup   ,plond*plev ,0.)
      call resetr(dup    ,plond*plev ,0.)
      call resetr(eup    ,plond*plev ,0.)
      call resetr(edn    ,plond*plev ,0.)
      call resetr(cape2d ,plond*plev ,0.)
      call resetr(mbup   ,plond      ,0.)
      call resetr(cape   ,plond      ,0.)

!
! Precompute frequently used variables
!
      do k = 1, plev
         do i = 1, plon
            rpdel(i,k) = 1. / pdelm1(i,k)
         end do
      end do
!
! Calculate geopotential height for cmfmca (requires virtual temperature)
!
      call zint(pilnm1  ,pmlnm1  ,rair    ,gravit  ,tvm      ,
     &          zm      ,zi      )      

      rrg    = rair / gravit
      do k = 1, plev
         do i = 1, plon
! SJL: finite-volume algorithm
            zm(i,k) = zi(i,k+1) + rrg * tvm(i,k)
     &            * ( 1. - pintm1(i,k) * (pilnm1(i,k+1) - pilnm1(i,k))
     &            * rpdel(i,k) )
         end do
      end do
!
! Since the PBL doesn't pass constituent perturbations, they 
! are zeroed here for input to the moist convection routine
!
      if( pcnst .ge. 2) then
        do m = 2, pcnst
          do i = 1, plon
            qpert(i,m) = 0.0
          end do
        end do
      endif
!
!=======================================================================
!                        Moist physics
!=======================================================================
!
! Begin with Zhang-McFarlane (1996) convection parameterization
!
      if(deep_conv) then
        call conv_ccm(tstar   ,qstar    ,precc   ,cnt     ,cnb     ,
     $               pblht   ,zm       ,phis    ,zi      ,zmdq    ,
     $               zmdt    ,pmidm1   ,pintm1  ,pdelm1  ,ts      ,
     $               .5*tdt  ,cmfmc    ,cme     ,nstep   ,lat     ,
     $               tpert   ,qpert    ,dlf     ,pflx    ,
     $               htlcl    ,htmmse  ,cape    ,cape2d  ,
     $               evp     ,mup      ,mdn     ,qlup    ,dup     ,
     $               eup     ,edn      ,mbup    ,zmeps)

!-- This code segment is added by SJL ---------------------------------
! SJL -- added re-evaporation of convective rain
! Algorithm based on Hack's ccm3.11.39 branch physics
! Evaporate some of the precip directly into the environment (Sundqvist)
! pflx = kg/m^2/s
! Setup some constants

      ke   = 2.4e-6
      rlovrcp = latvap / cpair
      tcrit = tmelt - 1.
      k1 = plev/6                 ! SJL

      call aqsat (tstar,  pmidm1,  est    ,qsat    ,plond,
     &            plon,   plev,    k1,     plev    )

      do i=1,plon
         sumflx(i) = 0.
      enddo

      do k=k1,plev
        do i=1,plon
        if ( tstar(i,k) > tcrit ) then
          if( precc(i) > 0.001*sumflx(i) ) then
!
! Re-evaporation of convective rain
!
          dpovrg  =  pdelm1(i,k)*rga
! Rain flux at the lower edge of the layer
          pflxtmp = max(pflx(i,k-1)-sumflx(i), 0.)
! Evaporation rate:
          envevap = max( ke*(1.-qstar(i,k,1)/qsat(i,k))*sqrt(pflxtmp), 0.)
! Evaporation is not to supersaturate the layer
! The tuning factor 0.5 is needed because qsat will become lower after evap
! In addition, evaporation of subgrid convective rain should not super-saturate
! the whole grid.
          envevap = max(min(envevap,0.5*(qsat(i,k)-qstar(i,k,1))*rtdt),0.)
! Should not evaporate more than it is available
          envevap = min(pflxtmp/dpovrg, envevap)
! The factor 0.999 is to prevent rounding-level negative precc
          envevap = min(envevap,0.999*(1000.*precc(i)-sumflx(i))/dpovrg)
          envevap = min(envevap,max(0.,0.999*(pflx(i,k-1)-sumflx(i))/dpovrg))
          pflx(i,k-1) = pflx(i,k-1) - envevap*dpovrg - sumflx(i)
          sumflx(i) = sumflx(i) + envevap*dpovrg
          qstar(i,k,1)  = qstar(i,k,1) + envevap*tdt
          tstar(i,k)  = tstar(i,k) - envevap*tdt*rlovrcp
          zmdt(i,k)   = zmdt(i,k) - envevap*rlovrcp
          zmdq(i,k)   = zmdq(i,k) + envevap
          cme(i,k)    = cme(i,k) - envevap
          cldtmp(i,k) = envevap
          endif
        endif
        end do
      end do
!
! adjust precc by the amount of precip evaporated
!
      do i = 1, plon
         precc(i) = precc(i) - sumflx(i)* 0.001
         pflx(i,plev) = pflx(i,plev) - sumflx(i)
      end do

! End re-evaporation of convective rain section
!---------------------------------------------------------------------

      call fillz(plon, plev, 1, qstar, pdelm1)

!
! Convert mass flux from reported mb/s to kg/m^2/s
!
      fac = 100. / gravit
      do k = 1, plev
        do i = 1, plon
          cmfmc(i,k) = cmfmc(i,k) * fac
        end do
      end do
 
c Add production of rain by conv_ccm to qc.  PJR says the vertical loop 
c indices are inconsistent.
 
      do k = 2, plev
        do i = 1, plon
          cmfdqr(i,k) = (pflx(i,k)-pflx(i,k-1))*gravit*rpdel(i,k)
        end do
      end do
! ------- End Zhang convection scheme -------------------------

      endif          ! end deep convection check
 
      do i = 1, plon
        tpert2(i) = 0.
      end do
      do m = 1, pcnst
        do i = 1, plon
           qpert2(i,m) = qpert(i,m)
        end do
      end do
!
! Call Hack (1994) convection scheme to deal with shallow/mid-level convection
! Begin by zeroing local copies of mass flux, energy fluxes, etc.
!
      call resetr(cmfmc2  ,plond*plev ,0.)
      call resetr(cmfdqr2 ,plond*plev ,0.)
      call resetr(cmfsl2  ,plond*plev ,0.)
      call resetr(cmflq2  ,plond*plev ,0.)
      call resetr(qc2     ,plond*plev ,0.)
      call resetr(cmfetr  ,plond*plev ,0.)
      call resetr(cmfdtr  ,plond*plev ,0.)
      call resetr(precc2  ,plond      ,0.)
 
      call cmfmca(lat     ,nstep   ,tdt     ,pmidm1  ,pdelm1  ,
     $            rpdel   ,zm      ,tpert2  ,qpert2  ,phis    ,
     $            pblht   ,tstar   ,qstar   ,cmfdt   ,cmfdq   ,
     $            cmfmc2  ,cmfdqr2 ,cmfsl2  ,cmflq2  ,precc2  ,
     $            qc2     ,cnt2    ,cnb2    ,cmfetr  ,cmfdtr  )

!
! Merge shallow/mid-level output with prior results from Zhang-McFarlane
!
 
      if( deep_conv ) then
        do i = 1, plon
          precc(i) = precc(i) + precc2(i)
          if (cnt2(i) .lt. cnt(i)) cnt(i) = cnt2(i)
          if (cnb2(i) .gt. cnb(i)) cnb(i) = cnb2(i)
        end do
      else
        do i = 1, plon
          precc(i) = precc2(i)
          cnt(i)   = cnt2(i)
          cnb(i)   = cnb2(i)
        end do
      endif
!
      do k = 1, plev
        do i = 1, plon
          cmfmc(i,k)  = cmfmc(i,k)  + cmfmc2(i,k)
          cmfdqr(i,k) = cmfdqr(i,k) + cmfdqr2(i,k)
        end do 
      end do
!
! Large-scale (stable) condensation procedure
! Evaporation of stratiform precipitation is included in this call
!
      call resetr(qc     ,plond*plev ,0.)
      call resetr(evapl  ,plond*plev ,0.)

      call cond(tdt    ,pmidm1  ,pdelm1  ,tstar   ,qstar   ,
     $          qc     ,precl   ,evapl)
!
      do k = 1, plev
        do i = 1, plon
          qstar(i,k,1) = max(1.0E-12, qstar(i,k,1))
        end do
      end do
!
!=======================================================================
!                     vertical diffusion
!=======================================================================
!
      do k = 1, plev-1
        do i = 1, plon
          rpdeli(i,k) = 1. / ( pmidm1(i,k+1) - pmidm1(i,k) )
        end do
      end do
 
      do m = 1, pcnst
        do k = 1, plev
           do i = 1, plon
             qold(i,k,m) = qstar(i,k,m)
           end do
         end do
      end do

      do k = 1, plev
        do i = 1, plon
          told(i,k) = tstar(i,k)
          tvm1(i,k) = told(i,k) * ( 1.0 + zvir * qold(i,k,1) )
          thm(i,k)  = told(i,k) * psomc(i,k)
        end do
      end do

!
! Compute geopotential height for PBL (requires virtual temperature)
!
      call zint(pilnm1  ,pmlnm1  ,rair    ,gravit  ,tvm1     ,
     &          zm      ,zi      )      

      do k = 1, plev
         do i = 1, plon
! SJL: finite-volume algorithm
            zm(i,k) = zi(i,k+1) + rrg * tvm1(i,k)
     &            * ( 1. - pintm1(i,k) * (pilnm1(i,k+1) - pilnm1(i,k))
     &            * rpdel(i,k) )
         end do
      end do

      call vdiff(lat     ,um1    ,vm1     ,told   ,qold    ,
     &           pmidm1  ,pintm1 ,rpdel   ,rpdeli ,tdt     ,
     &           thm     ,zm     ,taux    ,tauy   ,shflx   ,
     &           cflx    ,up1    ,vp1     ,thp    ,qstar   ,
     &           pblht   ,ustar  ,kvh     ,kvm    ,tpert   ,
     &           qpert   ,cgs    ,dtpblcg ,dqpblcg)
!
! Compute temperature after vertical diffusion
!
      do k = 1, plev
        do i = 1, plon
          tstar(i,k) = thp(i,k)/psomc(i,k)
        end do
      end do

      return
      end
