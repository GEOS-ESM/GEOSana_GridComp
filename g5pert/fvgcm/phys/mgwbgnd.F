#include <params.h>
      subroutine mgwbgnd (u, v, t, pm, pi, dpm, rdpm, piln, rlat,
     $     kldv, kldvmn, ksrc, ksrcmn, rdpldv, tau, ubi, ubm, xv, yv,
     $     ngwv, kbot, precl, precc)
C-----------------------------------------------------------------------
C 
C Driver for multiple gravity wave drag parameterization.
C 
C The parameterization is assumed to operate only where water vapor 
C concentrations are negligible in determining the density.
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
      integer
     $     kbot,                        ! index of bottom (source) interface
     $     ngwv                         ! number of gravity waves to use
      real 
     $     u(plond,plev),               ! midpoint zonal wind
     $     v(plond,plev),               ! midpoint meridional wind
     $     t(plond,plev),               ! midpoint temperatures
     $     pm(plond,plev),              ! midpoint pressures
     $     pi(plond,0:plev),            ! interface pressures
     $     dpm(plond,plev),             ! midpoint delta p (pi(k)-pi(k-1))
     $     rdpm(plond,plev),            ! 1. / (pi(k)-pi(k-1))
     $     piln(plond,0:plev),          ! ln(interface pressures)
     $     rlat                         ! current latitude in radians
      real precl(plond)
      real precc(plond)
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
     $     i,k,l                        ! loop indexes

      real
     $     tauback,                     ! background stress at c=0
     $     usrc(plond),                 ! u wind averaged over source region
     $     vsrc(plond),                 ! v wind averaged over source region
     $     al0,                         ! Used in lat dependence of GW spec. 
     $     dlat0,                       ! Used in lat dependence of GW spec.
     $     flat_gw,                     ! The actual lat dependence of GW spec.
     $     pi_g                         ! 3.14........
      real facl
      real facc
      real facw
  
C---------------------------------------------------------------------------
C Determine the source layer wind and unit vectors, then project winds.
C---------------------------------------------------------------------------
C 
C Just use the source level interface values for the source
c wind speed and direction (unit vector).
C 
      k = kbot
      do i = 1, plon
         ksrc(i) = k
         kldv(i) = k
         usrc(i) = 0.5*(u(i,k+1)+u(i,k))
         vsrc(i) = 0.5*(v(i,k+1)+v(i,k))
         ubi(i,kbot) = sqrt (usrc(i)**2 + vsrc(i)**2)
         xv(i) = usrc(i) / ubi(i,k)
         yv(i) = vsrc(i) / ubi(i,k)
      end do
C 
C Project the local wind at midpoints onto the source wind.
C 
      do k = 1, kbot
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
      do k = 1, kbot-1
         do i = 1, plon
            ubi(i,k) = 0.5 * (ubm(i,k) + ubm(i,k+1))
         end do
      end do
c +
c-----------------------------------------------------------------------
c Gravity wave sources
c-----------------------------------------------------------------------
c -
c +
c Determine the background stress at c=0
c -
      tauback = taubgnd * tauscal

c	Include dependence on latitude:
c 	The lat function was obtained by RR Garcia as 
c	currently used in his 2D model
c       [Added by F. Sassi on May 30, 1997]

        pi_g  = 4.*atan(1.)
        al0   = 40.*pi_g/180.
        dlat0 = 10.*pi_g/180.

        flat_gw =  0.5*(1.+tanh( (rlat-al0)/dlat0))
     $           + 0.5*(1.+tanh(-(rlat+al0)/dlat0)) 

!
! Source from convective activities added by S.-J. Lin
!
#if defined (GW_PRECIP)
        facl = 86400000./2.       ! conversion factor to mm/day normalized by  3
        facc = 86400000./40.      ! conversion factor to mm/day normalized by 60
#else
    	tauback = tauback * max(flat_gw, 0.2)
#endif

  
C Set the phase speeds and wave numbers in the direction of the source wind.
C Set the source stress magnitude (positive only, note that the sign of the 
C stress is the same as (c-u).
  
      do l = 0, ngwv
         facw = exp(-(c(l)/30.)**2)
         do i = 1, plon
#if defined (GW_PRECIP)
            tau(i, l,kbot) = facw * tauback * max( 0.1, 
     &                min(2.0, precl(i)*facl + precc(i)*facc) )
#else
            tau(i, l,kbot) = facw * tauback
#endif
            if( l .ne. 0 ) tau(i,-l,kbot) = tau(i, l,kbot)
         end do
      end do

c +
c Determine the min value of kldv and ksrc for limiting later loops
c and the pressure at the top interface of the low level stress divergence
c region.
c -
      ksrcmn = plev
      kldvmn = plev

      return
      end
