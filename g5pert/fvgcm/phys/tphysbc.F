#include <params.h>
      subroutine tphysbc(lat     ,clat    ,tdt     ,hbuf    ,
     $                   oro     ,phis    ,pmidm1  ,pintm1  ,
     $                   pdelm1  ,pmlnm1  ,pilnm1  ,omga    ,pblht   ,
     $                   psm1    ,tpert   ,tvm     ,tm1     ,qm1     ,
     $                   ts      ,qpert   ,precl   ,precc   ,precsl  ,
     $                   precsc  ,zm      ,asdir   ,asdif   ,
     $                   aldir   ,aldif   ,snowh   ,qrs     ,qrl     ,
     $                   flwds   ,fsns    ,lwup    ,qstar   ,tstar   ,
     $                   srfrad  ,sols    ,soll    ,solsd   ,solld ,
     $                   precst  ,absnxt  ,abstot  ,emstot  ,
     $                   hattr,   wspd)

C-----------------------------------------------------------------------
C
C Tendency physics BEFORE coupling to land, sea, and ice models.
C     o Apply dry convection, 
C     o moist convection
C     o large-scale condensation
C     o cloud calculations (cloud fraction, emissivity, etc.)
C     o radiation calculations
C Correct invalid (negative or too small) constituent mixing 
C ratios if necessary and report action in unit 6 output.
C Pass surface fields for separate surface flux calculations
C Dump appropriate fields to history file.
C
C----------------------------Code History-------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C                    T. Acker, March 1996
C Reviewed:          J. Kiehl, J. Hack, August 1992
C Include Ice:       J. Kiehl, January 1993
C Reviewed:          J. Kiehl, April 1996
C Merge Aphys:       J. Truesdale June 18 1997 
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <pagrid.h>
C------------------------------Commons----------------------------------
#include <comcon.h>
C-----------------------------------------------------------------------
#include <comctl.h>
C-----------------------------------------------------------------------
C SJL 1/3/99  #include <comhst.h>
C-----------------------------------------------------------------------
#include <comtim.h>
C------------------------------Diagnostics------------------------------
#include <diag.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real wspd(plon,plev)             ! wind speed (m/s) for ALT_CPS

#if (defined R4BYTE)
C 4-byte
      real*4 absnxt(plond,plev,4)      ! Nearest layer absorptivities
      real*4 abstot(plond,plevp,plevp) ! Non-adjacent layer absorptivites
      real*4 emstot(plond,plevp)       ! Total emissivity
#else
      real   absnxt(plond,plev,4)      ! Nearest layer absorptivities
      real   abstot(plond,plevp,plevp) ! Non-adjacent layer absorptivites
      real   emstot(plond,plevp)       ! Total emissivity
#endif

      real tvm(plond,plev)      ! Virtual temperature
      integer lat               ! Latitude index (S->N)
      real clat                 ! Current latitude(radians)
      real tdt                  ! 2 delta t (model time increment)

      integer hattr(5,pdiag)    ! History field attributes
      real*4  hbuf(*)           ! History buffer
      real oro(plond)           ! Sfc type flag (ocean,land,sea ice)
      real phis(plond)          ! Surface geopotential
      real pmidm1(plond,plev)   ! Midpoint pressures
      real pintm1(plond,plevp)  ! Interface pressures
      real pdelm1(plond,plev)   ! Layer pressure depth
      real pmlnm1(plond,plev)   ! ln(pmidm1)
      real pilnm1(plond,plevp)  ! ln(pintm1)
      real omga(plond,plev)     ! Vertical pressure velocity
      real pblht(plond)         ! Planetary boundary layer height
      real psm1(plond)          ! Surface pressure
      real tpert(plond)         ! Thermal temperature excess

      real tm1(plond,plev)      ! Temperature
      real qm1(plond,plev,pcnst+pnats) ! Specific humidity

      real ts(plond)                ! Surface temperature
      real qpert(plond,pcnst)       ! Thermal humidity & constituent excess
      real precl(plond)             ! Large-scale precipitation rate
      real precc(plond)             ! Convective-scale preciptn rate
      real precsl(plond)            ! L.S. precipitation rate (snow)
      real precsc(plond)            ! C.S. precipitation rate (snow)
      real dtcond(plond,plev)       ! dT/dt due to convective adjustments
      real dqcond(plond,plev,pcnst) ! dq/dt due to convective adjustments
      real cmfdt(plond,plev)        ! dT/dt due to moist convection
      real cmfdq(plond,plev)        ! dq/dt due to moist convection
      real cmfdqr(plond,plev)       ! dq/dt due to moist convective rainout 
      real cmfmc(plond,plev)        ! Moist convection cloud mass flux Total
      real cmfmc2(plond,plev)       ! Moist convection cloud mass flux Hack
      real cmfsl(plond,plev)        ! Moist convection lw stat energy flux
      real cmflq(plond,plev)        ! Moist convection total water flux
      real zmdt(plond,plev)         ! Zhang convective temperature tendency
      real zmdq(plond,plev)         ! Zhang convective moisture tendency
      real asdir(plond)         ! Albedo: shortwave, direct
      real asdif(plond)         ! Albedo: shortwave, diffuse
      real aldir(plond)         ! Albedo: longwave, direct
      real aldif(plond)         ! Albedo: longwave, diffuse

      logical deep_conv
      data    deep_conv / .true. /
C
C Output arguments
 
C SJL
      real tstar(plond,plev)    ! adjusted temperature
      real qstar(plond,plev,pcnst+pnats)      ! adjusted q fields
      real zm(plond,plev)       ! Height of midpoints (above surface)
C zm is not computed here
C SJL

      real snowh(plond)         ! Snow depth (liquid water equivalent)
      real qrs(plond,plev)      ! Shortwave heating rate
      real qrl(plond,plev)      ! Longwave  heating rate
      real flwds(plond)         ! Surface longwave down flux
      real fsns(plond)          ! Surface solar absorbed flux
      real lwup(plond)          ! Surface longwave up flux
      real srfrad(plond)        ! Net surface radiative flux (watts/m**2)
      real sols(plond)          ! Direct beam solar rad. onto srf (sw)
      real soll(plond)          ! Direct beam solar rad. onto srf (lw)
      real solsd(plond)         ! Diffuse solar radiation onto srf (sw)
      real solld(plond)         ! Diffuse solar radiation onto srf (lw)
C
C---------------------------Local workspace-----------------------------
C SJL: moved from output
      real precst(plond)        ! Total snow precip. rate
C SJL
C
C Most variables with appended 2 are a second copy of similar quantities
C declared above to accommodate calls to two convection schemes
C
      real coslat               ! Cosine latitude
      real cnt(plond)           ! Top level of convective activity
      real cnb(plond)           ! Bottom level of convective activity
      real qc(plond,plev)       ! dq/dt due to rainout terms
      real zi(plond,plevp)      ! Height above sfc in interface
      real precc2(plond)        ! Convective-scale preciptn rate
      real cnt2(plond)          ! Top level of convective activity
      real cnb2(plond)          ! Bottom level of convective activity
      real qc2(plond,plev)      ! dq/dt due to rainout terms
      real cmfdqr2(plond,plev)  ! dq/dt due to moist convective rainout 
      real cmfsl2(plond,plev)   ! Moist convection lw stat energy flux
      real cmflq2(plond,plev)   ! Moist convection total water flux
      real tpert2(plond)        ! Perturbation T
      real qpert2(plond,pcnst)  ! Perturbation q
      real dlf(plond,plev)      ! Detraining cld H20 from convection
      real pflx(plond,plev)     ! Conv rain flux thru out btm of lev
C JDC ADDED
      real htlcl(plond)         ! height above surface at LCL level
      real htmmse(plond)        ! height above surface at maximum moist static energy level
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

      real clwp(plond,plev)     ! Presribed cloud liq. h2o path 
      real rel(plond,plev)      ! Liquid cloud particle effective radius
      real rei(plond,plev)      ! Ice effective drop size (microns)
      real fice(plond,plev)     ! Fractional ice content within cloud
      real cld(plond,plevp)     ! Fraction cloud cover
      real effcld(plond,plevp)  ! Effective cloud=cld*emis
      real emis(plond,plev)     ! Cloud longwave emissivity
      real clc(plond)           ! Total convective cloud (cloud scheme)
      real rpdeli(plond,plev)   ! 1./(pmid(k+1)-pmid(k))
      real rpdel(plond,plev)    ! 1./(pint(k+1)-pint(k))
      real cltot(plond)         ! Diagnostic total cloud cover
      real cllow(plond)         !       "     low  cloud cover
      real clmed(plond)         !       "     mid  cloud cover
      real clhgh(plond)         !       "     hgh  cloud cover
      real coszrs(plond)        ! Cosine solar zenith angle
C        
      real actlwp(plond,plev)   ! Actual cloud liquid h2o path 
c      real totlwp(plond)       ! Actual column cloud liquid h2o path 
c      real sum                 ! Sum of actlwp
      real cloudup(plond,plev)  ! Cloud fraction during omega < 0.
C                                 
      real tpw(plond)           ! Total precipitable water
      real hl (plond)           ! Liquid water scale height

c GSL ... For averaging edge values to layer
      real cldtmp(plond,plev)   ! Temporary array
      real wk1(plond)           ! Temporary array
C SJL
      real fac
      real rtdt                 ! 1./tdt      SJL
      real  rrg
      integer  i,k,m            ! Lon, lev, constituent indices
!
!--- Rain re-evaporation section -----------------
      integer  k1
      real sumflx(plond)        ! flux integral
      real est(plond,plev)      ! Saturation vapor pressure
      real qsat(plond,plev) 
      real pflxtmp              ! Conv rain flux thru out btm of lev
      real ke
      real rlovrcp
      real dpovrg
      real envevap
      real tcrit

!-----------------------------------------------------------------------

      rtdt = 1./tdt
 
C Setup q and t accumulation fields
C
      do m=1,pcnst+pnats
        do k=1,plev
          do i=1,plon
            qstar(i,k,m) = max(0., qm1(i,k,m))
          end do
        end do
      end do

      do k=1,plev
        do i=1,plon
          tstar(i,k) = tm1(i,k)
        end do
      end do
C
C Zero out precip and convective fields before accumulating terms
C
      call resetr(precl  ,plond      ,0.)
      call resetr(precc  ,plond      ,0.)
      call resetr(precsl ,plond      ,0.)
      call resetr(precsc ,plond      ,0.)
      call resetr(qc     ,plond*plev ,0.)
      call resetr(cmfdt  ,plond*plev ,0.)
      call resetr(cmfdq  ,plond*plev ,0.)
      call resetr(zmdt   ,plond*plev ,0.)
      call resetr(zmdq   ,plond*plev ,0.)
      call resetr(cmfdqr ,plond*plev ,0.)
      call resetr(cmfmc  ,plond*plev ,0.)
      call resetr(cmfsl  ,plond*plev ,0.)
      call resetr(cmflq  ,plond*plev ,0.)
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

C
C Precompute frequently used variables
C
      do k=1,plev
         do i=1,plon
            rpdel(i,k) = 1./pdelm1(i,k)
         end do
      end do
C
      coslat = cos(clat)
C
C Calculate geopotential height for cmfmca (requires virtual temperature)
C
C
C Virtual temperatures, potential temperatures, and model level heights
C for use in radiation, surface, and diffusion calculations
C
C SJL tvm  already computed in the driver
C     call virtem(tstar   ,qstar   ,zvir    ,tvm )

      call zint(pilnm1  ,pmlnm1  ,rair    ,gravit  ,tvm      ,
     $          zm      ,zi      )      

      rrg = rair / gravit
      do k=1,plev
         do i=1,plon
! SJL: finite-volume algorithm
            zm(i,k) = zi(i,k+1) + rrg*tvm(i,k)*
     &          (1. - pintm1(i,k)*(pilnm1(i,k+1)-pilnm1(i,k))
     &            * rpdel(i,k) )
         end do
      end do
C
C Since the PBL doesn't pass constituent perturbations, they 
C are zeroed here for input to the moist convection routine
C
      if( pcnst .ge. 2) then
      do m=2,pcnst
        do i=1,plon
          qpert(i,m) = 0.0
        end do
      end do
      endif
 
! --------- Moist convection -------------

       if(deep_conv) then

#if defined (ALT_CPS)
!
! NCEP's SAS
!
! Determine TOP level for SAS
      do k=1,plev
         do i=1,plon
            if(pmidm1(i,k) > 2000. ) then   ! limit convection to 20 mb
               k1 = k
               go to 111
            endif
         enddo
      enddo 
111   k1 = plev - k1 + 1
      call sas_drv(plon, plat, plev, k1, tdt, coslat,
     &             qstar(1,1,1),tstar, psm1, pmidm1, wspd, 
     &             pintm1,  pilnm1,  rpdel, oro, omga,
     &             rearth, gravit, rair, rh2o, cpair, latvap, tmelt,
     &             precc, cnt, cnb)
#else
! Begin with Zhang-McFarlane (1996) convection parameterization
       call conv_ccm(tstar   ,qstar    ,precc   ,cnt     ,cnb     ,
     $	             pblht   ,zm       ,phis    ,zi      ,zmdq    ,
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
C
C Convert mass flux from reported mb/s to kg/m^2/s
C
      fac = 100./gravit
      do k=1,plev
        do i=1,plon
          cmfmc(i,k) = cmfmc(i,k) * fac
        end do
      end do
 
c Add production of rain by conv_ccm to qc.  PJR says the vertical loop 
c indices are inconsistent.
 
      do k=2,plev
        do i=1,plon
          cmfdqr(i,k) = (pflx(i,k)-pflx(i,k-1))*gravit*rpdel(i,k)
        end do
      end do

C JDC ADDED
c
c write out diagnostics from Zhang-McFarlane (1996) convection parameterization
c
      call outfld(hattr(1,iRNEVPDQ), cldtmp ,plond,lat,hbuf)
      do k=k1,plev
        do i=1,plon
          cldtmp(i,k) = -cldtmp(i,k) * rlovrcp
        end do
      end do
      call outfld(hattr(1,iRNEVPDT), cldtmp ,plond,lat,hbuf)
C
      do i = 1, plond
c Convert precip to mm/day
        wk1(i)  = precc(i) * 86400000.0
      end do

      call outfld(hattr(1,iZMDT),    zmdt     ,plond, lat, hbuf)
      call outfld(hattr(1,iZMDQ),    zmdq     ,plond, lat, hbuf)
      call outfld(hattr(1,iZMDQR),   cmfdqr   ,plond, lat, hbuf)
      call outfld(hattr(1,iZMCME),   cme      ,plond, lat, hbuf)
      call outfld(hattr(1,iZMDLF),   dlf      ,plond, lat, hbuf)
      call outfld(hattr(1,iZMPFLX),  pflx     ,plond, lat, hbuf)
      call outfld(hattr(1,iZMEVP),   evp      ,plond, lat, hbuf)
      call outfld(hattr(1,iZMMD),    mdn      ,plond, lat, hbuf)
      call outfld(hattr(1,iZMMU),    mup      ,plond, lat, hbuf)
      call outfld(hattr(1,iZMQL),    qlup     ,plond, lat, hbuf)
      call outfld(hattr(1,iZMDU),    dup      ,plond, lat, hbuf)
      call outfld(hattr(1,iZMEU),    eup      ,plond, lat, hbuf)
      call outfld(hattr(1,iZMED),    edn      ,plond, lat, hbuf)
      call outfld(hattr(1,iCAPE),    cape2d   ,plond, lat, hbuf)
      call outfld(hattr(1,iZMEPS),   zmeps    ,plond, lat, hbuf)

C 2D variables

      call outfld(hattr(1,iHTLCL),   htlcl    ,plond, lat, hbuf)
      call outfld(hattr(1,iHTMMSE),  htmmse   ,plond, lat, hbuf)
      call outfld(hattr(1,iZMMB),    mbup     ,plond, lat, hbuf)
      call outfld(hattr(1,iCAPEMX),  cape     ,plond, lat, hbuf)
      call outfld(hattr(1,iZMPR),    wk1      ,plond, lat, hbuf)

#endif 
! ------- End Zhang convection scheme -------------------------

      endif          ! end deep convection check

 
      do i=1,plon
        tpert2(i) = 0.
      end do
      do m=1,pcnst
        do i=1,plon
           qpert2(i,m) = qpert(i,m)
        end do
      end do
C
C Call Hack (1994) convection scheme to deal with shallow/mid-level convection
C Begin by zeroing local copies of mass flux, energy fluxes, etc.
C
      call resetr(cmfmc2  ,plond*plev ,0.)
      call resetr(cmfdqr2 ,plond*plev ,0.)
      call resetr(cmfsl2  ,plond*plev ,0.)
      call resetr(cmflq2  ,plond*plev ,0.)
      call resetr(qc2     ,plond*plev ,0.)
      call resetr(cmfetr  ,plond*plev ,0.)
      call resetr(cmfdtr  ,plond*plev ,0.)
      call resetr(precc2  ,plond      ,0.)
C
      call cmfmca(lat     ,nstep   ,tdt     ,pmidm1  ,pdelm1  ,
     $            rpdel   ,zm      ,tpert2  ,qpert2  ,phis    ,
     $            pblht   ,tstar   ,qstar   ,cmfdt   ,cmfdq   ,
     $            cmfmc2  ,cmfdqr2 ,cmfsl2  ,cmflq2  ,precc2  ,
     $            qc2     ,cnt2    ,cnb2    ,cmfetr  ,cmfdtr  )
 
!
!--------------------------------------------------------------
! Release of latent heat due to ice-phase physics added here by SJL
!     do k=k1,plev
!       do i=1,plon
!       if(tstar(i,k)<tcrit .and. tstar(i,plev)<tcrit) then
!         tstar(i,k) = tstar(i,k) + hice*max(0., qc2(i,k))
!       endif
!       end do
!     end do
!--------------------------------------------------------------
!
c Merge shallow/mid-level output with prior results from Zhang-McFarlane
 
      if(deep_conv) then
        do i=1,plon
          precc(i) = precc(i) + precc2(i)
          if (cnt2(i) .lt. cnt(i)) cnt(i) = cnt2(i)
          if (cnb2(i) .gt. cnb(i)) cnb(i) = cnb2(i)
        end do
      else
        do i=1,plon
          precc(i) = precc2(i)
          cnt(i) = cnt2(i)
          cnb(i) = cnb2(i)
        end do
      endif
c
      do k=1,plev
        do i=1,plon
C cmfmc is in [kg/m2/s]
          cmfmc(i,k)  = cmfmc(i,k)  + cmfmc2(i,k)
C convert [kg/m2/s] to [pa/s]
          cldtmp(i,k) = cmfmc(i,k)  * gravit
          cmfmc2(i,k) = cmfmc2(i,k) * gravit
          cmfetr(i,k) = cmfetr(i,k) * gravit
          cmfdtr(i,k) = cmfdtr(i,k) * gravit
          cmfdqr(i,k) = cmfdqr(i,k) + cmfdqr2(i,k)
          cmfsl(i,k)  = cmfsl(i,k)  + cmfsl2(i,k)
          cmflq(i,k)  = cmflq(i,k)  + cmflq2(i,k)
        end do 
      end do

c Convert precip to mm/day
      do i=1, plon
          wk1(i) = precc(i) * 86400000.
      enddo

      call outfld(hattr(1,iCMFDT),    cmfdt,    plond, lat, hbuf)
      call outfld(hattr(1,iCMFDQ),    cmfdq,    plond, lat, hbuf)
      call outfld(hattr(1,iCMFMC2),   cmfmc2,   plond, lat, hbuf)
      call outfld(hattr(1,iCMFMC),    cldtmp,   plond, lat, hbuf)
      call outfld(hattr(1,iCMFDQR2),  cmfdqr2,  plond, lat, hbuf)
      call outfld(hattr(1,iCMFETR),   cmfetr,   plond, lat, hbuf)
      call outfld(hattr(1,iCMFDTR),   cmfdtr,   plond, lat, hbuf)
      call outfld(hattr(1,iPRECC),    wk1,      plond, lat, hbuf)

C
C Large-scale (stable) condensation procedure
C Evaporation of stratiform precipitation is included in this call
C
      call resetr(qc     ,plond*plev ,0.)
      call resetr(evapl  ,plond*plev ,0.)

      call cond(tdt    ,pmidm1  ,pdelm1  ,tstar   ,qstar   ,
     $          qc     ,precl   ,evapl)
C
C Compute rate of temperature change due to convective processes
C
      do k=1,plev
        do i=1,plon
          dtcond(i,k) = (tstar(i,k) - tm1(i,k))*rtdt
        end do
      end do
C
        do k=1,plev
          do i=1,plon
            qstar(i,k,1) = max(1.0E-12, qstar(i,k,1))
          end do
        end do

C  Compute rate of constituent change due to convective processes
C
      do m=1,pcnst
        do k=1,plev
          do i=1,plon
            dqcond(i,k,m) = (qstar(i,k,m) - qm1(i,k,m))*rtdt
          end do
        end do
      end do
C
C determine whether precipitation, prec, is frozen (snow) or not
C
      do i=1,plon
        if (tstar(i,plev  ).lt.tmelt .and. ts(i).lt.tmelt .and. 
     $      tstar(i,plev-1).lt.tmelt) then
          precsc(i) = precc(i)
          precsl(i) = precl(i)
        else
          precsc(i) = 0.
          precsl(i) = 0.
        end if
      end do
C
      do i=1, plon
          wk1(i) = precl(i) * 86400000.
      enddo
C
      call outfld(hattr(1,iPRECL),     wk1,    plond, lat, hbuf)
      call outfld(hattr(1,iEVAPL),     evapl,  plond, lat, hbuf)
      call outfld(hattr(1,iDQRL),      qc,     plond, lat, hbuf)
      call outfld(hattr(1,iDTCOND),    dtcond, plond, lat, hbuf)
      call outfld(hattr(1,iDQCOND),    dqcond, plond, lat, hbuf)
C
      if (dosw .or. dolw) then

C Moved here by  SJL 
C Virtual temperatures, potential temperatures, and model level heights
C for use in radiation, surface, and diffusion calculations
 
      call virtem(tstar    ,qstar  ,zvir    ,tvm     )

C SJL 8/24/98
C zm is no longer needed as input  to cldint with the new (finite-volume)
C algorithm for zint.F
C
C SJL 
C 1./delta-p needed by cloud routines
         do k=1,plev-1
            do i=1,plon
               rpdeli(i,k) = 1./(pmidm1(i,k+1) - pmidm1(i,k))
            end do
         end do
C
C Compute cloud properties for input to radiation
         call cldint(pmidm1  ,rpdeli  ,tstar   ,qstar   ,omga    ,
     $               cnt     ,cnb     ,pilnm1  ,pmlnm1  ,tvm     ,
     $               zm      ,cld     ,clc     ,clwp    ,emis    ,
     $               effcld  ,cltot   ,cllow   ,clmed   ,clhgh   ,
     $               oro     ,rel     ,rei     ,fice    ,pdelm1  ,
     $               tpw     ,hl      ,cmfmc   ,psm1    ,snowh,
     &               coslat  ,rhclr   ,convcld ,stratcld )
C
C Dump cloud field information to history tape buffer (diagnostics)
C
C GSL 1-4-99: output CLOUD and EFFCLD with PLEV instead of PLEVP levels
C SEN 3-22-01: removed erroneous average for CLOUD and EFFCLD 
C
         do k = 1, plev
           do i = 1, plon
             relhum(i,k) = rhclr(i,k)*(1.0-cld(i,k))+cld(i,k)
             cldtmp(i,k) = cld(i,k)
             relhum(i,k) = 100.0 * relhum(i,k)
             rhclr(i,k)  = 100.0 * rhclr(i,k)
           end do
         end do
         call outfld(hattr(1,iCLOUD), cldtmp ,plond,lat,hbuf)
         do k = 1, plev
           do i = 1, plon
             cldtmp(i,k) = effcld(i,k)
           end do
         end do
         call outfld(hattr(1,iEFFCLD),   cldtmp,   plond, lat, hbuf)
         call outfld(hattr(1,iCNVCLD),   clc,      plond, lat, hbuf)
         call outfld(hattr(1,iCLDTOT),   cltot,    plond, lat, hbuf)
         call outfld(hattr(1,iCLDLOW),   cllow,    plond, lat, hbuf)
         call outfld(hattr(1,iCLDMED),   clmed,    plond, lat, hbuf)
         call outfld(hattr(1,iCLDHGH),   clhgh,    plond, lat, hbuf)
         call outfld(hattr(1,iLWSH),     hl,       plond, lat, hbuf)
         call outfld(hattr(1,iRELHUM),   relhum,   plond, lat, hbuf)
         call outfld(hattr(1,iRHCLR),    rhclr,    plond, lat, hbuf)
         call outfld(hattr(1,iCONVCLD),  convcld,  plond, lat, hbuf)
         call outfld(hattr(1,iSTRATCLD), stratcld, plond, lat, hbuf)
C
C Special diagnostic cloud water fields:
C
         call outfld(hattr(1,iSETLWP),   clwp,     plond, lat, hbuf)
C
C Outfield special monthly cloud diagnostics
C      
C Commented out by SJL
C SJL
C        if (ninavg(1).eq.'Q') then
c           do k=1,plev
c              do i=1,plon
c                 cloudup(i,k) = 0.
c                 if (omga(i,k).lt.0.) then
c                    cloudup(i,k) = cld(i,k)
c                 end if
c              end do
c           end do
c           call outfld(hattr(1,iCLOUDUP),cloudup ,plond,lat,hbuf)
c        end if
C SJL
C
C Compute actual liquid water path
C
         do k=1,plev
           do i=1,plon
             actlwp(i,k) = cld(i,k)*clwp(i,k)
           end do
         end do
         call outfld(hattr(1,iCLDLWP),actlwp ,plond,lat,hbuf)
c         do i=1,plon
c           sum = 0.
c           do k=1,plev
c             sum = sum + actlwp(i,k)
c           end do
c           totlwp(i) = sum
c         end do
c         call outfld(hattr(1,iTOTLWP),totlwp ,plond,lat,hbuf)
C
C Cosine solar zenith angle for current time step
C
         call zenith (calday  ,dodiavg ,clat    ,coszrs  )
C
C Complete radiation calculations
C
         call radctl(hbuf    ,clat    ,coslat  ,lat     ,lwup    ,
     $               pmidm1  ,pintm1  ,pmlnm1  ,pilnm1  ,tstar   ,
     $               qstar   ,cld     ,effcld  ,clwp    ,coszrs  ,
     $               asdir   ,asdif   ,aldir   ,aldif   ,fsns    ,
     $               qrs     ,qrl     ,flwds   ,rel     ,rei     ,
     $               fice    ,sols    ,soll    ,solsd   ,solld   ,
     $               absnxt  ,abstot  ,emstot  ,hattr)
 
         call outfld(hattr(1,iQRL),qrl,plon,lat,hbuf)
         call outfld(hattr(1,iQRS),qrs,plon,lat,hbuf)
         call outfld(hattr(1,iALDIF),aldif,plon,lat,hbuf)
         call outfld(hattr(1,iALDIR),aldir,plon,lat,hbuf)
         call outfld(hattr(1,iASDIF),asdif,plon,lat,hbuf)
         call outfld(hattr(1,iASDIR),asdir,plon,lat,hbuf)      
      end if

C Compute net surface radiative flux for use by surface temperature code.
C Note that units have already been converted to mks in RADCTL.  Since
C fsns and flwds are in the buffer, array values will be carried across
C timesteps when the radiation code is not invoked.
C
      do i=1,plon
         srfrad(i) = fsns(i) + flwds(i)
      end do
      call outfld(hattr(1,iSRFRAD),srfrad,plond,lat,hbuf)
C
C Save atmospheric fields to force surface models
C
C
C Calculate total snow precipitation rate:
C
      do i=1,plon
         precst(i) = precsc(i) + precsl(i) 
      end do

C
C Add radiative heating rates to temperature
C

      do k = 1,plev
         do i = 1,plon
            tstar(i,k) = tstar(i,k) + tdt*(qrs(i,k)+qrl(i,k))
         end do
      end do
 
      return
 
      end
