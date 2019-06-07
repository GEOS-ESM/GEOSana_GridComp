#include <params.h>
      subroutine mgwintr (u, v, t, q, sgh, pm, pi, dpm, rdpm, piln,
     &                    zm, pblh, dt, hbuf, lat, rlat, ut, vt, tt,
     &                    dtv, precl, precc, hattr)

C-----------------------------------------------------------------------
C
C Interface for multiple gravity wave drag parameterization.
C
C-----------------------------------------------------------------------
      implicit none
C-----------------------------------------------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <mgw.h>
#include <diag.h>
C-----------------------------------------------------------------------
C
C Input variables
C
      integer
     $     lat                          ! latitude index
      real 
     $     u(plond,plev),               ! midpoint zonal wind
     $     v(plond,plev),               ! midpoint meridional wind
     $     t(plond,plev),               ! midpoint temperatures
     $     q(plond,plev),               ! midpoint specific humidities
     $     sgh(plond),                  ! standard deviation of orography
     $     pm(plond,plev),              ! midpoint pressures
     $     pi(plond,0:plev),            ! interface pressures
     $     piln(plond,0:plev),          ! ln(interface pressures)
     $     dpm(plond,plev),     ! midpoint delta p (pi(k)-pi(k-1))
     $     rdpm(plond,plev),    ! 1. / (pi(k)-pi(k-1))
     $     zm(plond,plev),              ! midpoint heights
     $     pblh(plond),         ! planetary boundary layer height
     $     dt,                          ! time step
     $     dtv(plond,plev),             ! t tendency from vertical diffusion
     $     rlat                         ! current latitude in radians
      real precl(plond)
      real precc(plond)

      integer hattr(5,pdiag)            ! History field attributes
      real*4  hbuf(*)                   ! History buffer
C
C Output variables
C
      real
     $     ut(plond,plev),              ! zonal wind tendency
     $     vt(plond,plev),              ! meridional wind tendency
     $     tt(plond,plev)               ! temperature tendency
C
C Local workspace
C
      integer
     $     i,k                          ! loop indexes
      integer
     $     kldv(plond),                 ! index of top interface of low level stress 
     $                                  ! divergence region
     $     kldvmn,                      ! min value of kldv
     $     ksrc(plond),                 ! index of top interface of source region
     $     ksrcmn                       ! min value of ksrc
      real
     $     ni(plond,0:plev),            ! interface Brunt-Vaisalla frequency
     $     nm(plond,plev),              ! midpoint Brunt-Vaisalla frequency
     $     rdpldv(plond),               ! 1/dp across low level divergence region
     $     rhoi(plond,0:plev),          ! interface density
     $     tau(plond,-pgwv:pgwv,0:plev), ! wave Reynolds stress
     $     tau0x(plond),                ! c=0 sfc. stress (zonal)
     $     tau0y(plond),                ! c=0 sfc. stress (meridional)
     $     ti(plond,0:plev),            ! interface temperature
     $     ttgw(plond,plev),            ! temperature tendency
     $     ubi(plond,0:plev),           ! projection of wind at interfaces
     $     ubm(plond,plev),             ! projection of wind at midpoints
     $     utgw(plond,plev),            ! zonal wind tendency
     $     vtgw(plond,plev),            ! meridional wind tendency
     $     xv(plond),                   ! unit vectors of source wind (x)
     $     yv(plond)                    ! unit vectors of source wind (y)
      real dudt(plond,plev)
      real dvdt(plond,plev)
      real dtdt(plond,plev)

C
C-----------------------------------------------------------------------------
c+
c Profiles of background state variables
c-
      call mgwprof (u, v, t, pm, pi, rhoi, ni, ti, nm)
C
C-----------------------------------------------------------------------------
C Non-orographic backgound gravity wave spectrum
C-----------------------------------------------------------------------------
c+
c Determine the wave source for a background spectrum at ~100 mb
c-
      call mgwbgnd (u, v, t, pm, pi, dpm, rdpm, piln, rlat,
     $     kldv, kldvmn, ksrc, ksrcmn, rdpldv, tau, ubi, ubm, xv, yv,
     $     PGWV, kbotbg, precl, precc)
c+
c Solve for the drag profile
c-
      call mgwdrag (PGWV, kbotbg, ktopbg, u, v, t, pi, dpm, rdpm, piln,
     $     rhoi, ni, ti, nm, dt,
     $     kldv, kldvmn, ksrc, ksrcmn, rdpldv, tau, ubi, ubm, xv, yv,
     $     utgw, vtgw, tau0x, tau0y)
c+
c Scale Tendencies
c-
      do k =1 , plev
	do i = 1 , plon
	   dudt(i,k) = utgw(i,k)*effgw
	   dvdt(i,k) = vtgw(i,k)*effgw
        enddo
      enddo
c+
c Compute the temperature tendency from energy conservation.
c Add the tendencies to the model arrays.
c-
      do k = 1, plev
         do i = 1, plon
            dtdt(i,k) = -(u(i,k)*dudt(i,k) + v(i,k)*dvdt(i,k))
     $                  /(cp*(1.0 + cpvir*q(i,k)))

            ut(i,k) = ut(i,k) + dudt(i,k)
            vt(i,k) = vt(i,k) + dvdt(i,k)
            tt(i,k) = tt(i,k) + dtdt(i,k)
            dtv(i,k)= dtv(i,k)+ dtdt(i,k)
! ------------------
! SJL: update u,v, t
! ------------------
!           u(i,k) = u(i,k) + dudt(i,k)*dt
!           v(i,k) = v(i,k) + dvdt(i,k)*dt
!           t(i,k) = t(i,k) + dtdt(i,k)*dt
!-------------------------------------------
         end do
      end do
 
C-----------------------------------------------------------------------------
C Orographic stationary gravity wave
C-----------------------------------------------------------------------------
c+
c Determine the orographic wave source
c-
      call mgworo (u, v, t, sgh, pm, pi, dpm, zm, nm, pblh,
     $     kldv, kldvmn, ksrc, ksrcmn, rdpldv, tau, ubi, ubm, xv, yv)
c+
c Solve for the drag profile
c-
      call mgwdrag (0, kbotoro, ktoporo, u, v, t, pi, dpm, rdpm, piln,
     $     rhoi, ni, ti, nm, dt,
     $     kldv, kldvmn, ksrc, ksrcmn, rdpldv, tau, ubi, ubm, xv, yv,
     $     utgw, vtgw, tau0x, tau0y)
c+
c Scale Tendencies
c-
      do k =1 , plev
	do i = 1 , plon
	   utgw(i,k)=utgw(i,k)*effgw
	   vtgw(i,k)=vtgw(i,k)*effgw
        enddo
      enddo
c+
c Compute the temperature tendency from energy conservation.
c Add the tendencies to the model arrays.
c-
      do k = 1, plev
         do i = 1, plon
            ttgw(i,k) = -(u(i,k)*utgw(i,k) + v(i,k)*vtgw(i,k))
     $           /(cp*(1.0 + cpvir*q(i,k)))

            ut(i,k) = ut(i,k) + utgw(i,k)
            vt(i,k) = vt(i,k) + vtgw(i,k)
            tt(i,k) = tt(i,k) + ttgw(i,k)
            dtv(i,k)= dtv(i,k)+ ttgw(i,k)
!-------------------
! Update tendencies:
!-------------------
            utgw(i,k) = utgw(i,k) + dudt(i,k)
            vtgw(i,k) = vtgw(i,k) + dvdt(i,k)
            ttgw(i,k) = ttgw(i,k) + dtdt(i,k)
! ------------------
! SJL: update u,v, t
! ------------------
            u(i,k) = u(i,k) + utgw(i,k)*dt
            v(i,k) = v(i,k) + vtgw(i,k)*dt
            t(i,k) = t(i,k) + ttgw(i,k)*dt
         end do
      end do
 
! Write output fields to history file
 
      call outfld (hattr(1,iUTGW),  utgw, plond, lat, hbuf)
      call outfld (hattr(1,iVTGW),  vtgw, plond, lat, hbuf)
      call outfld (hattr(1,iTTMGW), ttgw, plond, lat, hbuf)

      call outfld (hattr(1,iTAUGWX), tau0x, plond, lat, hbuf)
      call outfld (hattr(1,iTAUGWY), tau0y, plond, lat, hbuf)

      return
      end
