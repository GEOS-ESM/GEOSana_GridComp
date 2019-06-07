#include <params.h>
      subroutine tphysac(lat     ,hbuf    ,ztodt   ,pmidm1  ,pintm1  ,
     $                   pdelm1  ,pilnm1  ,prsk,   prske,   qq,
     &                   t2      ,fu      ,fv      ,
     $                   ustar   ,vstar   ,tstar   ,qstar   ,psomc   ,
     $                   thm     ,zm      ,pblh    ,qpert   ,
     $                   tpert   ,clat    ,shflx   ,taux    ,
     $                   tauy    ,cflx    ,sgh     ,lhflx   ,oro     ,
     $                   ts      ,tref    ,precl, precc   ,hattr,   gwdf,  qbeg)
C-----------------------------------------------------------------------
C
C Tendency physics after coupling to land, sea, and ice models.
C Call physics subroutines and compute the following: 
C     o 10 m reference height quantities
C     o Radon surface flux and decay (optional)
C     o Vertical diffusion and planetary boundary layer
C     o Rayleigh friction
C     o Multiple gravity wave drag
C Dump appropriate fields to history file.
C
C----------------------------Code History-------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C                    T. Acker, March 1996
C Reviewed:          J. Kiehl, J. Hack, August 1992
C Reviewed:          J. Kiehl April 1996
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
C------------------------------Diagnostics------------------------------
#include <diag.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer lat                 ! Model latitude index
      integer hattr(5,pdiag)      ! History field attributes
      real*4  hbuf(*)             ! History buffer
      real ztodt                  ! Two times model timestep (2 delta-t)
      real pmidm1(plond,plev)     ! Midpoint pressures
      real pintm1(plond,plevp)    ! Interface pressures
      real pdelm1(plond,plev)     ! Layer pressure depth
      real pilnm1(plond,plevp)    ! ln(pintm1)
      real t2(plond,plev)         ! Temperature tendency
      real q2(plond,plev,pcnst)   ! q tendency
      real fu(plond,plev)         ! Accumulated u tendency (from tphys)
      real fv(plond,plev)         ! Accumulated v tendency (from tphys)
      real tstar(plond,plev)        ! Temperature
      real qm1(plond,plev,pcnst+pnats)! Specific humidity and constituents
      real psomc(plond,plev)      ! (psm1/pmidm1)**cappa
      real thm(plond,plev)        ! Potential temperature 
      real zm(plond,plev)         ! Height of midpoints (above surface) 

!-----CMT : added variables for the 2.5level turbulence scheme-----------
      real prsk(plond,plev)      ! p**k at midpoints
      real prske(plond,plevp)    ! pe**k at edge levels
      real qq(plond,plev)        ! turbulent kinetic energy
      logical  qbeg

      real clat                   ! Current latitude(radians)
      real oro(plond)             ! Land/ocean/sea ice flag
      real ts(plond)              ! Land/ocean/sea surface temperature
      real tref(plond)            ! Surface air temp. (not skin temp)
      real precl(plond)
      real precc(plond)           ! convective precip
C
C Output arguments
C
      real pblh(plond)            ! Planetary boundary layer height
      real qpert(plond,pcnst)     ! Moisture/constit. perturbation (PBL)
      real tpert(plond)           ! Temperature perturbation (PBL)
      real shflx(plond)           ! Sensible heat flux (w/m^2)
      real taux(plond)            ! X surface stress (zonal)
      real tauy(plond)            ! Y surface stress (meridional)
      real cflx(plond,pcnst)      ! Surface constituent flux (kg/m^2/s)
      real sgh(plond)             ! Std. deviation of orography for gwd
      real lhflx(plond)           ! Latent heat flux (w/m^2)
C
C---------------------------Local workspace-----------------------------
C
      real ustar(plond,plev)        ! U horizontal wind component
      real vstar(plond,plev)        ! V horizontal wind component
      real qstar(plond,plev,pcnst)
      real up1(plond,plev)        ! U-wind after vdiff & rfric
      real vp1(plond,plev)        ! V-wind after vdiff & rfric
      real dtv(plond,plev)        ! T tendency from vertical diffusion
      real duv(plond,plev)        ! U tendency from vertical diffusion
      real dvv(plond,plev)        ! V tendency from vertical diffusion
      real rpdeli(plond,plev)     ! 1./(pmid(k+1)-pmid(k))
      real rpdel(plond,plev)      ! 1./(pint(k+1)-pint(k))
      real srcsnk(plond,plev,pcnst) ! Src/sink for spec. hum. & constit.
      integer i,k                 ! Longitude, level indices
      integer m                   ! Constituent index
      integer npts                ! No. of ocean points of sea ice points
      integer indx(plond)         ! Long. index array of ocn or sea ice pts
      integer nptssi              ! Number of sea ice points
      integer indxsi(plond)       ! Longitude index array of sea ice points
      integer nptsoc              ! Number of ocean points
      integer indxoc(plond)       ! Longitude index array of ocean points
      integer nptsexc             ! Number of excess points
      integer indxexc(plond)      ! Longitude index array of excess points

! NRL-GWD
      integer kbot
      real rdt
      real hdt
      real pkm(plon,plev)
      real ptv(plon,plev)
      real hprime(plon)
      real utgw(plon,plev)
      real vtgw(plon,plev)
      real ttgw(plon,plev)
C
C--------------------- Statement functions -----------------------------
C
      logical seaice
      logical ocean
      logical land
      logical gwdf     !Bowen
      seaice(i) = nint(oro(i)).eq.2
      ocean(i)  = nint(oro(i)).eq.0
      land(i)   = nint(oro(i)).eq.1  
C
C-----------------------------------------------------------------------

C Precompute frequently used variables  
C
      do k=1,plev
         do i=1,plon
            rpdel(i,k) = 1./pdelm1(i,k)
         end do
      end do
      do k=1,plev-1
         do i=1,plon
            rpdeli(i,k) = 1./(pmidm1(i,k+1) - pmidm1(i,k))
         end do
      end do
C
C Set up index arrays of ocean, sea ice. Set tref over
C land to 1.e36. Note, the 2m reference temperature over 
C land is output in the lsm history tape.
C
      npts   = plon
      nptssi = 0
      nptsoc = 0
      do i=1,plon
        indx(   i) = i
        if (ocean(i)) then
          nptsoc = nptsoc + 1
          indxoc(nptsoc) = i
        else if (seaice(i)) then
          nptssi = nptssi + 1
          indxsi(nptssi) = i
        else if (land(i)) then
           tref(i) = 1.e36
        endif
      end do
C
C Check if latent heat flux over ocean or sea ice exceeds the total
C moisture content of the lowest model layer, thereby creating 
C negative moisture.
C
C SJL: qneg4 was modified to prevent multitasking problem on sgi

      call qneg4('TPHYSAC '   ,lat           ,indx  ,npts    ,ztodt   ,
     $           qstar(1,plev,1),rpdel(1,plev) ,shflx ,lhflx   ,cflx  ,
     $           indxexc ,nptsexc)
C
C Determine 10m reference temperature over ocean
C
      call trefoce (indxoc, nptsoc, pmidm1(1,plev), ustar(1,plev),
     $              vstar(1,plev),
     $              tstar(1,plev), qstar(1,plev,1), thm(1,plev), 
     $                zm(1,plev),      ts,
     $              tref )
C
C Determine 10m reference temperature over sea ice
C
      call trefsice (indxsi, nptssi, pmidm1(1,plev), ustar(1,plev),
     $                vstar(1,plev),
     $              tstar(1,plev), qstar(1,plev,1), thm(1,plev), 
     $                zm(1,plev),       ts,
     $              tref )
      if ( pcnst.ge.2 ) then
C SJL
C Initialize source/sink array to zero for all constituents.
C
      call resetr(srcsnk,plndlv*pcnst,0.)
        if ( settrace.eq.'test' .or. settrace.eq.'TEST' ) then
c
c Add radon surface flux (from WCRP).
c
          m = 2
C SJL     call rnsfwcrp( clat, oro, cflx(1,m) )                !SJL
c
c Radon decay.
c
C SJL     call rndecay( qstar(1,1,m), ztodt, srcsnk(1,1,m) )   !SJL
          if ( pcnst.ge.4 ) then
c
c Zero strat tracer in bottom layer.
c
            m = 4
            call resetr( qstar(1,plev,m), plon, 0. )
          end if
        end if

C Add contribution of src/snk tendency to constituent fields
C (excluding moisture).
 
      do m=2,pcnst
         do k=1,plev
            do i=1,plon
               qstar(i,k,m) = qstar(i,k,m) + ztodt*srcsnk(i,k,m)
            end do
         end do
      end do
      end if
 
#if defined (ALT_PBL)
!

!----vertical diffusion with the 2.5level turbulence scheme
!----march 2001
!----COMPUTES REST OF TURBULENCE DAO level 2.5 Turbulence ------

      call turb(lat, hbuf, ztodt, pmidm1, pintm1,
     $          pdelm1,ustar,vstar,tstar, qstar,
     $          taux,tauy,shflx, cflx,fu,fv,t2,
     $          duv,dvv,dtv,up1,vp1,
     $          qbeg,prsk,prske,qq,hattr,pblh)

!----MODIFIED VERSION, ONLY CALCULATES PBL,QPERT,TPERT----
!----Computes pblh(pblh height) and temperature and humidity
!    excess with NCAR diffusion scheme.
 
! Update thm for pbldif_simpl
      do k=1,plev
         do i=1,plon
            thm(i,k) = tstar(i,k)*psomc(i,k)
         end do
      end do
c***************************************************************
c    In this version pbl is caluclated in turb.F with MH alg.
c    while tpert and qpert are still done according to NCAR
c
c      call pbldif_simpl_tq(lat, hbuf, thm  ,qstar  ,zm   ,up1,
c     $            vp1, tstar     ,pmidm1  ,cflx    ,shflx ,
c     $            taux    ,tauy   ,pblh    ,tpert ,
c     $            qpert, hattr)
c*************************************************************


      call pbldif_simpl(lat, hbuf, thm  ,qstar  ,zm   ,up1,
     $            vp1, tstar     ,pmidm1  ,cflx    ,shflx ,
     $            taux    ,tauy   ,pblh    ,tpert ,
     $            qpert, hattr)

#else
! NCAR CCM3: Vertical diffusion/pbl calculation
 
      call vdintr(hbuf    ,pmidm1  ,pintm1  ,psomc   ,thm     ,
     $            zm      ,rpdel   ,rpdeli  ,ustar   ,vstar   ,
     $            tstar   ,taux    ,tauy    ,shflx   ,cflx    ,
     $            ztodt   ,lat     ,qstar   ,fu      ,fv      ,
     $            t2      ,dtv     ,duv     ,dvv     ,up1     ,
     $            vp1     ,pblh    ,tpert   ,qpert   ,hattr)

#endif
 
! Gravity wave drag

      if (gwdf)    then
#if defined (ALT_GWD)

! ------------ Begine ATL-GWD section ---------------------------
! NRL-GWD settings:

      do i=1,plon
!        hprime(i) = (pintm1(i,plevp)/100000.)**cappa
      enddo

      do k=1,plev
         do i=1,plon
            pkm(i,k) = (pmidm1(i,k)/100000.)**cappa
!           pkm(i,k) = hprime(i) / psomc(i,k)
            ptv(i,k) = tstar(i,k)*(1.+zvir*qstar(i,k,1))/pkm(i,k)
         enddo
      enddo

      do i=1,plon
	 hprime(i) = 0.5*sgh(i)
      enddo

      kbot = plev - 1

      do k=1,plev
         do i=1,plon
            ustar(i,k) = up1(i,k)
            vstar(i,k) = vp1(i,k)
         enddo
      enddo

! The following is a very simple GWD code from NRL based on Palmer et al. 1986
      call nrl_gwd(plon,   plev,  kbot, ztodt, cpair, gravit,
     &             pmidm1, pdelm1, pkm, ptv, hprime, up1, vp1)
!
! To be implemented: GFDL-GWD (Pierrehumbert and Stern)

      rdt = 1.  / ztodt
      hdt = 0.5 * ztodt
      do k=1,plev
         do i=1,plon
            utgw(i,k) = (up1(i,k) - ustar(i,k) ) * rdt
            vtgw(i,k) = (vp1(i,k) - vstar(i,k) ) * rdt
            ttgw(i,k) = -((ustar(i,k)+hdt*utgw(i,k))*utgw(i,k) +
     &                    (vstar(i,k)+hdt*vtgw(i,k))*vtgw(i,k) )
     &                  / (cpair*(1.0 + cpvir*qstar(i,k,1)))
               t2(i,k) = t2(i,k) + ttgw(i,k)
            tstar(i,k) = tstar(i,k) + ttgw(i,k)*ztodt
              fu(i,k) = fu(i,k) + utgw(i,k)
              fv(i,k) = fv(i,k) + vtgw(i,k)
         enddo
      enddo
      call outfld (hattr(1,iUTGW),   utgw, plon, lat, hbuf)
      call outfld (hattr(1,iVTGW),   vtgw, plon, lat, hbuf)
      call outfld (hattr(1,iTTMGW),  ttgw, plon, lat, hbuf)

! ------------ end ATL-GWD section -----------------------------

#else
! ------------ NCAR GWD --------------------------------------
      call mgwintr(up1  ,vp1   ,tstar   ,qstar   ,sgh     ,
     $            pmidm1  ,pintm1  ,pdelm1  ,rpdel   ,pilnm1  ,
     $            zm      ,pblh    ,ztodt   ,hbuf    ,lat     ,
     $            clat    ,fu      ,fv      ,t2      ,dtv     ,
     $            precl, precc,  hattr)
! ------------ NCAR GWD --------------------------------------
#endif
      endif

! SJL 06/09/2001: (up1, vp1) are updated and (fu,fv) accumulated
      do k=1,plev
         do i=1,plon
            ustar(i,k) = up1(i,k)
            vstar(i,k) = vp1(i,k)
         enddo
      enddo
 
! Dump diffusive tendencies to history tape buffer
 
      call outfld (hattr(1,iDTV),dtv,plond,lat,hbuf)
 
      return
      end
