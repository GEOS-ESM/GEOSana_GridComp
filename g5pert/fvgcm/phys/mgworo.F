#include <params.h>
      subroutine mgworo (u, v, t, sgh, pm, pi, dpm, zm, nm, pblh,
     $     kldv, kldvmn, ksrc, ksrcmn, rdpldv, tau, ubi, ubm, xv, yv)
C-----------------------------------------------------------------------
C 
C Orographic source for multiple gravity wave drag parameterization.
c 
c The stress is returned for a single wave with c=0, over orography.
c For points where the orographic variance is small (including ocean),
c the returned stress is zero. 
C 
C-----------------------------------------------------------------------
      implicit none
C-----------------------------------------------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <mgw.h>
C-----------------------------------------------------------------------
C 
C Input variables
C 
      real 
     $     u(plond,plev),               ! midpoint zonal wind
     $     v(plond,plev),               ! midpoint meridional wind
     $     t(plond,plev),               ! midpoint temperatures
     $     sgh(plond),                  ! standard deviation of orography
     $     pm(plond,plev),              ! midpoint pressures
     $     pi(plond,0:plev),            ! interface pressures
     $     dpm(plond,plev),             ! midpoint delta p (pi(k)-pi(k-1))
     $     zm(plond,plev),              ! midpoint heights
     $     nm(plond,plev),              ! midpoint Brunt-Vaisalla frequency
     $     pblh(plond)                  ! planetary boundary layer height
C 
C Output variables
C 
      integer
     $     kldv(plond),                 ! index of top interface of low level stress 
     $                                  ! divergence region
     $     kldvmn,                      ! min value of kldv
     $     ksrc(plond),                 ! index of top interface of source region
     $     ksrcmn                       ! min value of ksrc

      real
     $     rdpldv(plond),               ! 1/dp across low level divergence region
     $     tau(plond,-pgwv:pgwv,0:plev), ! wave Reynolds stress
     $     ubi(plond,0:plev),           ! projection of wind at interfaces
     $     ubm(plond,plev),             ! projection of wind at midpoints
     $     xv(plond),                   ! unit vectors of source wind (x)
     $     yv(plond)                    ! unit vectors of source wind (y)
C 
C Local workspace
C 
      integer
     $     i,k                          ! loop indexes

      real
     $     sghmax,                      ! max orographic sdv to use
     $     tauoro(plond),               ! c=0 stress from orography
     $     zldv(plond),                 ! top of the low level stress divergence region
     $     nsrc(plond),                 ! b-f frequency averaged over source region
     $     psrc(plond),                 ! interface pressure at top of source region
     $     rsrc(plond),                 ! density averaged over source region
     $     usrc(plond),                 ! u wind averaged over source region
     $     vsrc(plond)                  ! v wind averaged over source region
C 
C---------------------------------------------------------------------------
c +
c Average the basic state variables for the wave source over the depth of
c the orographic standard deviation. Here we assume that the apropiate
c values of wind, stability, etc. for determining the wave source are 
c averages over the depth of the atmosphere pentrated by the typical mountain.
c Reduces to the bottom midpoint values when sgh=0, such as over ocean.
c 
c Also determine the depth of the low level stress divergence region, as
c the max of the boundary layer depth and the source region depth. This
c can be done here if the stress magnitude does not determine the depth,
c otherwise it must be done below.
c -
      k = plev
      do i = 1, plon
         ksrc(i) = k-1
         kldv(i) = k-1
         psrc(i) = pi(i,k-1)
         rsrc(i) = pm(i,k)/(r*t(i,k)) * dpm(i,k)
         usrc(i) = u(i,k) * dpm(i,k)
         vsrc(i) = v(i,k) * dpm(i,k)
         nsrc(i) = nm(i,k)* dpm(i,k)
      end do
      do k = plev-1, plev/2, -1
         do i = 1, plon
!SJL, 06/15/2001: bug fix
!           if (0.5*sgh(i) .gt. sqrt(zm(i,k)*zm(i,k+1))) then
            if (sgh(i) .gt. sqrt(zm(i,k)*zm(i,k+1))) then
               ksrc(i) = k-1
               kldv(i) = k-1
               psrc(i) = pi(i,k-1)
               rsrc(i) = rsrc(i)
     $              + pm(i,k) / (r*t(i,k))* dpm(i,k)
               usrc(i) = usrc(i) + u(i,k) * dpm(i,k)
               vsrc(i) = vsrc(i) + v(i,k) * dpm(i,k)
               nsrc(i) = nsrc(i) + nm(i,k)* dpm(i,k)
            elseif (pblh(i) .gt. sqrt(zm(i,k)*zm(i,k+1))) then
               kldv(i) = k-1
            end if
         end do
      end do
      do i = 1, plon
         rsrc(i) = rsrc(i) / (pi(i,plev) - psrc(i))
         usrc(i) = usrc(i) / (pi(i,plev) - psrc(i))
         vsrc(i) = vsrc(i) / (pi(i,plev) - psrc(i))
         nsrc(i) = nsrc(i) / (pi(i,plev) - psrc(i))

         ubi(i,plev) = sqrt (usrc(i)**2 + vsrc(i)**2)
         xv(i) = usrc(i) / ubi(i,plev)
         yv(i) = vsrc(i) / ubi(i,plev)
      end do
C 
C Project the local wind at midpoints onto the source wind.
C 
      do k = 1, plev
         do i = 1, plon
            ubm(i,k) = u(i,k) * xv(i) + v(i,k) * yv(i)
         end do
      end do
C 
C Compute the interface wind projection by averaging the midpoint winds.
C Use the top level wind at the top interface.
C 
      do i = 1, plon
         ubi(i,0) = ubm(i,1)
      end do
      do k = 1, plev-1
         do i = 1, plon
            ubi(i,k) = 0.5 * (ubm(i,k) + ubm(i,k+1))
         end do
      end do
c +
c Determine the orographic c=0 source term following McFarlane (1987).
c Set the source top interface index to plev, if the orographic term is zero.
c -
      do i = 1, plon
         if ((ubi(i,plev) .gt. orovmin) .and.
     $        (sgh(i) .gt. orohmin)) then
            sghmax = orofc2 * (ubi(i,plev) / nsrc(i))**2
            tauoro(i) = oroeko2 * min(sgh(i)**2, sghmax) *
     $           rsrc(i) * nsrc(i) * ubi(i,plev)
         else
            tauoro(i) = 0.
            ksrc(i) = plev
            kldv(i) = plev
         end if
      end do
C 
C Set the phase speeds and wave numbers in the direction of the source wind.
C Set the source stress magnitude (positive only, note that the sign of the 
C stress is the same as (c-u).
C 
      do i = 1, plon
         tau(i,0,plev) = tauoro(i)
      end do
c +
c Find the top interface of the low level stress divergence region according
c to the maximum depth of three criterion.
c 1. source region depth
c 2. planetary boundary layer depth
c 3. 10 * (u_*) / N where u_* is defined from the gravity wave stresss
c = sqrt(tau/rho) using source region values
c -
cif (kbot .lt. plev) then
cdo i = 1, plon
ckldv(i) = kbot
cend do
celse
cdo i = 1, plon
czldv(i) = max (pblh(i), sgh(i)
czldv(i) = max (zdlv(i),
c           zldvcon * sqrt(tau(i,0,k)/rsrc(i)) / nsrc(i))
cend do
ckldv(i) = plev-1
cdo k = plev-1, plev/2, -1
cdo i = 1, plon
cif (zldv(i) .gt. sqrt(zm(i,k)*zm(i,k+1))) then
ckldv(i)  = k-1
cend do
cend do
cend if
c +
c Determine the min value of kldv and ksrc for limiting later loops
c and the pressure at the top interface of the low level stress divergence
c region.
c -
      ksrcmn = plev
      kldvmn = plev
      do i = 1, plon
         ksrcmn = min(ksrcmn, ksrc(i))
         kldvmn = min(kldvmn, kldv(i))
         if (kldv(i) .ne. plev) then
            rdpldv(i) = 1. / (pi(i,kldv(i)) - pi(i,plev))
         end if
      end do
      if (fracldv .le. 0.) kldvmn = plev

      return
      end

