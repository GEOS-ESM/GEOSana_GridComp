#include <params.h>
      subroutine pbldif_simpl(lat, hbuf, th  ,q   ,z   ,u ,
     $                  v, t    ,pmid    ,cflx    ,shflx   ,
     $                  taux    ,tauy    ,pblh    ,tpert   ,
     $                  qpert, hattr )
C------------------------------------------------------------------------
c---------Modified version, computes only pblh,tpert and qpert.
c         This subroutine is called from vdiff_simpl.(Version used
c         when the new turbulence code is used in order to
c         supply pblh, tpert and qpert to the convective scheme.
C------------------------------------------------------------------------
C 
C Atmospheric boundary layer computation.
C
C Nonlocal scheme that determines eddy diffusivities based on a
C diagnosed boundary layer height and a turbulent velocity scale;
C also, countergradient effects for heat and moisture, and constituents
C are included, along with temperature and humidity perturbations which 
C measure the strength of convective thermals in the lower part of the 
C atmospheric boundary layer.
C
C For more information, see Holtslag, A.A.M., and B.A. Boville, 1993:
C Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
C Model. J. Clim., vol. 6., p. 1825--1842.
c
c Updated by Holtslag and Hack to exclude the surface layer from the
c definition of the boundary layer Richardson number. Ri is now defined
c across the outer layer of the pbl (between the top of the surface
c layer and the pbl top) instead of the full pbl (between the surface and
c the pbl top). For simiplicity, the surface layer is assumed to be the
c region below the first model level (otherwise the boundary layer depth 
c determination would require iteration).
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Boville
C Standardized:      J. Rosinski, June 1992
C Reviewed:          B. Boville, P. Rasch, August 1992
C Reviewed:          B. Boville, P. Rasch, April 1996
C
C Modified for boundary layer height diagnosis: Bert Holtslag, june 1994
C >>>>>>>>>  (Use ricr = 0.3 in this formulation)
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
#include <compbl.h>
C-----------------------------------------------------------------------
#include <comvd.h>
C------------------------------Diagnostics------------------------------
#include <diag.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer hattr(5,pdiag)      ! History field attributes
      real*4  hbuf(*)             ! History buffer

      real th(plond,plev)          ! potential temperature [K]
      real q(plond,plev,pcnst)     ! specific humidity [kg/kg]
      real z(plond,plev)           ! height above surface [m]
      real u(plond,plev)           ! windspeed x-direction [m/s]
      real v(plond,plev)           ! windspeed y-direction [m/s]
      real t(plond,plev)           ! temperature (used for density)
      real pmid(plond,plev)        ! midpoint pressures
      real cflx(plond,pcnst)       ! surface constituent flux (kg/m2/s)
      real shflx(plond)            ! surface heat flux (W/m2)
      real taux(plond)             ! surface u stress (N)
      real tauy(plond)             ! surface v stress (N)
C
      integer lat                  ! latitude index
c
C Output arguments
C
      real pblh(plond)             ! boundary-layer height [m]
      real tpert(plond)            ! convective temperature excess
      real qpert(plond)            ! convective humidity excess
C
C---------------------------Local parameters----------------------------
C
      real tiny                    ! lower bound for wind magnitude
      parameter (tiny=1.e-36)
C
C---------------------------Local workspace-----------------------------
C
      integer i                    ! longitude index
      integer k                    ! level index
      integer m                    ! constituent index

      real ustar(plond)            ! surface friction velocity [m/s]

      real heatv(plond)            ! surface virtual heat flux
      real thvsrf(plond)           ! sfc (bottom) level virtual temperature
      real thvref(plond)           ! reference level virtual temperature
      real tkv                     ! model level potential temperature
      real therm(plond)            ! thermal virtual temperature excess
      real phiminv(plond)          ! inverse phi function for momentum
      real wm(plond)               ! turbulent velocity scale for momentum
      real vvk                     ! velocity magnitude squared
      real khfs(plond)             ! surface kinematic heat flux [mK/s]
      real kqfs(plond,pcnst)       ! sfc kinematic constituent flux [m/s]
c      real vv2(plond)              ! velocity magnitude squared 
      real rino(plond,plev)        ! bulk Richardson no. from level to ref lev
      real tlv(plond)              ! ref. level pot tmp + tmp excess
      real rrho(plond)             ! 1./bottom level density (temporary)
      real obklen(plond)           ! Obukhov length
      real ustr                    ! unbounded ustar
      real term                    ! intermediate calculation
      real fac                     ! interpolation factor
      real pblmin                  ! min pbl height due to mechanical mixing

      logical unstbl(plond)        ! pts w/unstbl pbl (positive virtual ht flx)
      logical check(plond)         ! True=>chk if Richardson no.>critcal
C
C Compute kinematic surface fluxes
C
      do i=1,plon
        rrho(i) = rair*t(i,plev)/pmid(i,plev)
        ustr = sqrt(sqrt(taux(i)**2 + tauy(i)**2)*rrho(i))
        ustar(i) = max(ustr,0.01)
        khfs(i) = shflx(i)*rrho(i)/cpair
        kqfs(i,1) = cflx(i,1)*rrho(i)
      end do
      do m=2,pcnst
        do i=1,plon
          kqfs(i,m)= cflx(i,m)*rrho(i)
        end do
      end do
C
C Compute various arrays for use later:
C
      do i=1,plon
        thvsrf(i) = th(i,plev)*(1.0 + 0.61*q(i,plev,1))
        heatv(i)  = khfs(i) + 0.61*th(i,plev)*kqfs(i,1)
        wm(i)     = 0.
        therm(i)  = 0.
        qpert(i)  = 0.
        tpert(i)  = 0.
        obklen(i) = -thvsrf(i)*ustar(i)**3/
     $             (g*vk*(heatv(i) + sign(1.e-10,heatv(i))))
      end do
C
C
C >>>> Define first a new factor fac=100 for use in Richarson number
C      Calculate virtual potential temperature first level
C      and initialize pbl height to z1
C
      fac = 100.
C
      do i=1,plon
        thvref(i) = th(i,plev)*(1.0 + 0.61*q(i,plev,1))
        pblh(i) = z(i,plev)
        check(i) = .true.
C
C Initialization of lowest level Ri number 
C (neglected in initial Holtslag implementation)
C
C--      vv2(i) = fac*ustar(i)**2
C--      vv2(i) = max(vv2(i),tiny)
C--      rino(i,plev) = g*(thvsrf(i) - thvref(i))*z(i,plev)/
C--  $                  (thvref(i)*vv2(i))
        rino(i,plev) = 0.0
      end do
C
C PBL height calculation:
C Search for level of pbl. Scan upward until the Richardson number between
C the first level and the current level exceeds the "critical" value.
C
      do k=plev-1,plev-npbl+1,-1
        do i=1,plon
          if (check(i)) then
            vvk = (u(i,k) - u(i,plev))**2 + (v(i,k) - v(i,plev))**2
     $           + fac*ustar(i)**2
            vvk = max(vvk,tiny)
            tkv = th(i,k)*(1. + 0.61*q(i,k,1))
            rino(i,k) = g*(tkv - thvref(i))*(z(i,k)-z(i,plev))
     $                 /(thvref(i)*vvk)
            if(rino(i,k).ge.ricr) then
              pblh(i) = z(i,k+1) + (ricr - rino(i,k+1))/
     $                 (rino(i,k) - rino(i,k+1))*(z(i,k) - z(i,k+1))
              check(i) = .false.
            end if
          end if
        end do
      end do
C
C Set pbl height to maximum value where computation exceeds number of
C layers allowed
C
      do i=1,plon
        if (check(i)) pblh(i) = z(i,plevp-npbl)
      end do
C
C Improve estimate of pbl height for the unstable points.
C Find unstable points (virtual heat flux is positive):
C
      do i=1,plon
        if (heatv(i) .gt. 0.) then
          unstbl(i) = .true.
          check(i) = .true.
        else
          unstbl(i) = .false.
          check(i) = .false.
        end if   
      end do
C
C For the unstable case, compute velocity scale and the
C convective temperature excess:
C
      do i=1,plon
        if (check(i)) then
          phiminv(i) = (1. - binm*pblh(i)/obklen(i))**onet
          wm(i)= ustar(i)*phiminv(i)
          therm(i) = heatv(i)*fak/wm(i)       
C--         rino(i,plev) = -g*therm(i)*z(i,plev)/(thvref(i)*vv2(i))
          rino(i,plev) = 0.0
          tlv(i) = thvref(i) + therm(i)
        end if
      end do
C
C Improve pblh estimate for unstable conditions using the
C convective temperature excess:
C
      do k=plev-1,plev-npbl+1,-1
        do i=1,plon
          if (check(i)) then
            vvk = (u(i,k) - u(i,plev))**2 + (v(i,k) - v(i,plev))**2
     $           + fac*ustar(i)**2
            vvk = max(vvk,tiny)
            tkv = th(i,k)*(1. + 0.61*q(i,k,1))
            rino(i,k) = g*(tkv - tlv(i))*(z(i,k)-z(i,plev))
     $                 /(thvref(i)*vvk)
            if(rino(i,k).ge.ricr) then
              pblh(i) = z(i,k+1) + (ricr - rino(i,k+1))/
     $                 (rino(i,k) - rino(i,k+1))*(z(i,k) - z(i,k+1))
              check(i) = .false.
            end if
          end if
        end do
      end do
C
C Points for which pblh exceeds number of pbl layers allowed;
C set to maximum
C
      do i=1,plon
        if (check(i)) pblh(i) = z(i,plevp-npbl)
      end do
C
C PBL height must be greater than some minimum mechanical mixing depth
C Several investigators have proposed minimum mechanical mixing depth
C relationships as a function of the local friction velocity, u*.  We 
C make use of a linear relationship of the form h = c u* where c=700.
C The scaling arguments that give rise to this relationship most often 
C represent the coefficient c as some constant over the local coriolis
C parameter.  Here we make use of the experimental results of Koracin 
C and Berkowicz (1988) [BLM, Vol 43] for wich they recommend 0.07/f
C where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
C latitude value for f so that c = 0.07/f = 700.
C
      do i=1,plon
        pblmin  = 700.0*ustar(i)
        pblh(i) = max(pblh(i),pblmin)
      end do
C
C pblh is now available; do preparation for diffusivity calculation:
C
      do i=1,plon
C
C Do additional preparation for unstable cases only, set temperature
C and moisture perturbations depending on stability.
C
        if (unstbl(i)) then
          phiminv(i) = (1. - binm*pblh(i)/obklen(i))**onet
          wm(i)      = ustar(i)*phiminv(i)
          tpert(i)   = max(khfs(i)*fak/wm(i),0.)   
          qpert(i)   = max(kqfs(i,1)*fak/wm(i),0.)    
        else
          tpert(i)   = max(khfs(i)*fak/ustar(i),0.) 
          qpert(i)   = max(kqfs(i,1)*fak/ustar(i),0.) 
        end if
      end do


      call outfld(hattr(1,iPBLH), pblh ,plond,lat,hbuf)
      call outfld(hattr(1,iTPERT),tpert,plond,lat,hbuf)
      call outfld(hattr(1,iQPERT),qpert,plond,lat,hbuf)

      return
      end
