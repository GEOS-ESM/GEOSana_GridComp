#include <params.h>
      subroutine vdintr(hbuf    ,pmidm1  ,pintm1  ,psomc   ,thm     ,
     $                  zm      ,rpdel   ,rpdeli  ,um1     ,vm1     ,
     $                  tm1     ,taux    ,tauy    ,shflx   ,cflx    ,
     $                  ztodt   ,lat     ,qm1     ,fu      ,fv      ,
     $                  t2      ,dtv     ,duv     ,dvv     ,up1     ,
     $                  vp1     ,pblh    ,tpert   ,qpert   ,hattr)
C-----------------------------------------------------------------------
C
C interface routine for vertical diffusion and pbl scheme
C
C calling sequence:
C
C    vdinti        initializes vertical diffustion constants
C    pbinti        initializes pbl constants
C     .
C     .
C    vdintr        interface for vertical diffusion and pbl scheme
C      vdiff       performs vert diff and pbl
C        pbldif    boundary layer scheme
C        mvdiff    diffuse momentum
C        qvdiff    diffuse constituents
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski
C Standardized:      J. Rosinski, June 1992
C Reviewed:          P. Rasch, B. Boville, August 1992
C Reviewed:          P. Rasch, April 1996
C Reviewed:          B. Boville, April 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C------------------------------Commons----------------------------------
#include <comtrcnm.h>
C-----------------------------------------------------------------------
#include <comvd.h>
C------------------------------Diagnostics------------------------------
#include <diag.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
C SJL
      integer hattr(5,pdiag)      ! History field attributes
      real*4  hbuf(*)             ! History buffer
      real pmidm1(plond,plev)     ! midpoint pressures
      real pintm1(plond,plevp)    ! interface pressures
      real psomc(plond,plev)      ! (psm1/pmidm1)**cappa
      real thm(plond,plev)        ! potential temperature midpoints
      real zm(plond,plev)         ! midpoint geopotential height above sfc
      real rpdel(plond,plev)      ! 1./pdel (thickness between interfaces)
      real rpdeli(plond,plev)     ! 1./pdeli (thickness between midpoints)
      real um1(plond,plev)        ! u-wind input
      real vm1(plond,plev)        ! v-wind input
      real tm1(plond,plev)        ! temperature input
      real taux(plond)            ! x surface stress (n)
      real tauy(plond)            ! y surface stress (n)
      real shflx(plond)           ! surface sensible heat flux (w/m2)
      real cflx(plond,pcnst)      ! surface constituent flux (kg/m2/s)
      real ztodt                  ! 2 delta-t

      integer lat                 ! latitude index
C
C Input/output arguments
C
      real qm1(plond,plev,pcnst)  ! initial/final constituent field
      real fu(plond,plev)         ! cumulative u tendency
      real fv(plond,plev)         ! cumulative v tendency
      real t2(plond,plev)         ! cumulative T tendency
C
C Output arguments
C
      real dtv(plond,plev)        ! temperature tendency (heating)
      real duv(plond,plev)        ! u-wind tendency
      real dvv(plond,plev)        ! v-wind tendency
      real up1(plond,plev)        ! u-wind after vertical diffusion
      real vp1(plond,plev)        ! v-wind after vertical diffusion
      real pblh(plond)            ! planetary boundary layer height
      real tpert(plond)           ! convective temperature excess
      real qpert(plond,pcnst)     ! convective humidity and constituent excess
C JDC ADDED
      real dtpblcg(plond,plev)    ! T tendency due to counter gradient term
      real dqpblcg(plond,plev)    ! q tendency due to counter gradient term
C
C---------------------------Local storage-------------------------------
C
      integer i,k,m               ! longitude,level,constituent indices

      real dqv(plond,plev,pcnst)  ! constituent diffusion tendency
      real ustar(plond)           ! surface friction velocity
      real kvh(plond,plevp)       ! diffusion coefficient for heat
      real kvm(plond,plevp)       ! diffusion coefficient for momentum
      real denom                  ! denominator of expression
      real qp1(plond,plev,pcnst)  ! constituents after vdiff
      real thp(plond,plev)        ! potential temperature after vdiff
      real cgs(plond,plevp)       ! counter-gradient star (cg/flux)
      real rztodt                 ! 1./ztodt
C SJL
      real tmp                    ! 
C
C-----------------------------------------------------------------------
C
C Call vertical diffusion code. No setup work is required.
C
      if(ntopfl .ne. 1) then
      do k=1,ntopfl
        do i=1,plon
          dtv(i,k) = 0.
          up1(i,k) = um1(i,k)
          vp1(i,k) = vm1(i,k)
        enddo
      enddo
      endif
      call vdiff(lat     ,um1     ,vm1     ,tm1     ,qm1     ,
     $           pmidm1  ,pintm1  ,rpdel   ,rpdeli  ,ztodt   ,
     $           thm     ,zm      ,taux    ,tauy    ,shflx   ,
     $           cflx    ,up1     ,vp1     ,thp     ,qp1     ,
     $           pblh    ,ustar   ,kvh     ,kvm     ,tpert   ,
     $           qpert   ,cgs     ,dtpblcg ,dqpblcg )
C
C Convert the diffused fields back to diffusion tendencies.
C Add the diffusion tendencies to the cummulative physics tendencies,
C except for constituents. The diffused values of the constituents
C replace the input values.

C
      rztodt = 1./ztodt
      do k=ntopfl,plev
        do i=1,plon
          duv(i,k) = (up1(i,k) - um1(i,k))*rztodt
          dvv(i,k) = (vp1(i,k) - vm1(i,k))*rztodt
          fu(i,k)  = fu(i,k) + duv(i,k)
          fv(i,k)  = fv(i,k) + dvv(i,k)
          denom    = cpair*(1. + cpvir*qm1(i,k,1))
C SJL
C The following frictional heating rate has a 1st order truncation
C error as compared to the exact discrete form.
C
c         dtv(i,k) = -(um1(i,k)*duv(i,k) + vm1(i,k)*dvv(i,k))/denom +
c    $                (thp(i,k)/psomc(i,k) - tm1(i,k))*rztodt
c         t2(i,k)  = t2(i,k) + dtv(i,k)

          tmp = tm1(i,k)          ! Save original value
C The exact form of temperature "increment" due to frictional
C heating:
C tm1 will be the final adjusted temperature.
C
          tm1(i,k) = thp(i,k)/psomc(i,k) 
     &            + ( (up1(i,k)+um1(i,k))*(um1(i,k)-up1(i,k)) +
     &                (vp1(i,k)+vm1(i,k))*(vm1(i,k)-vp1(i,k))
     &                 ) / (2.*denom)
 
C Total heating rate due to PBL diffusion:
          dtv(i,k) =  (tm1(i,k)-tmp) * rztodt
C Total (cumulative) heating rate:
          t2(i,k)  = t2(i,k) + dtv(i,k)
C JDC ADDED
C tendency of potential temperature --> tendency of temperature
          dtpblcg(i,k) = dtpblcg(i,k) / psomc(i,k)
        end do
        do m=1,pcnst
          do i=1,plon
            dqv(i,k,m) = (qp1(i,k,m) - qm1(i,k,m))*rztodt
            qm1(i,k,m) = qp1(i,k,m)
          end do
        end do
      end do
C
C Save the pbl output variables on the history file. The other
C vertical diffusion variables are saved later, after Rayleigh 
C friction tendencies have been added.
C
      call outfld(hattr(1,iPBLH), pblh ,plond,lat,hbuf)
      call outfld(hattr(1,iTPERT),tpert,plond,lat,hbuf)
      call outfld(hattr(1,iQPERT),qpert,plond,lat,hbuf)
      call outfld(hattr(1,iUSTAR),ustar,plond,lat,hbuf)
C
C If you want to uncomment the following outfld calls you must modify
C the argument list to vdiff
C
      call outfld(hattr(1,iCGS),      cgs(1,2),   plond, lat, hbuf)
      call outfld(hattr(1,iDTPBLCG),  dtpblcg,    plond, lat, hbuf)
      call outfld(hattr(1,iDQPBLCG),  dqpblcg,    plond, lat, hbuf)
      call outfld(hattr(1,iKVH),      kvh(1,2),   plond, lat, hbuf)
      call outfld(hattr(1,iKVM),      kvm(1,2),   plond, lat, hbuf)
      call outfld(hattr(1,iDUV),      duv,        plond, lat, hbuf)
      call outfld(hattr(1,iDVV),      dvv,        plond, lat, hbuf)
      call outfld(hattr(1,iVD01),     dqv(1,1,1), plond, lat, hbuf)
c      do m=2,pcnst
cc        call outfld(vdiffnam(m),dqv(1,1,m),plond,lat,hbuf)
c         call outfld(hattr(1,iVD01+m-1),dqv(1,1,m),plond,lat,hbuf)
c      end do
C
      return
      end
 
