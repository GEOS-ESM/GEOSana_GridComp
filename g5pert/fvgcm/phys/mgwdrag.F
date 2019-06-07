#include <params.h>
      subroutine mgwdrag (ngwv, kbot, ktop, u, v, t, pi, dpm, rdpm, 
     $     piln, rhoi, ni, ti, nm, dt,
     $     kldv, kldvmn, ksrc, ksrcmn, rdpldv, tau, ubi, ubm, xv, yv,
     $     ut, vt, tau0x, tau0y)
 
! Solve for the drag profile from the multiple gravity wave drag
! parameterization.
! 1. scan up from the wave source to determine the stress profile
! 2. scan down the stress profile to determine the tendencies
!     => apply bounds to the tendency
!          a. from wkb solution
!          b. from computational stability constraints
!     => adjust stress on interface below to reflect actual bounded tendency
!-----------------------------------------------------------------------
      implicit none
#include <pmgrid.h>
#include <mgw.h>
 
! Input variables
 
      integer
     $     kbot,                        ! index of bottom (source) interface
     $     ktop,                        ! index of top interface of gwd region
     $     ngwv                         ! number of gravity waves to use
      integer
     $     kldv(plond),                 ! index of top interface of low level stress 
     $                                  ! divergence region
     $     kldvmn,                      ! min value of kldv
     $     ksrc(plond),                 ! index of top interface of source region
     $     ksrcmn                       ! min value of ksrc
      real 
     $     u(plond,plev),               ! midpoint zonal wind
     $     v(plond,plev),               ! midpoint meridional wind
     $     t(plond,plev),               ! midpoint temperatures
     $     pi(plond,0:plev),            ! interface pressures
     $     dpm(plond,plev),             ! midpoint delta p (pi(k)-pi(k-1))
     $     rdpm(plond,plev),            ! 1. / (pi(k)-pi(k-1))
     $     piln(plond,0:plev),          ! ln(interface pressures)
     $     rhoi(plond,0:plev),          ! interface density
     $     ni(plond,0:plev),            ! interface Brunt-Vaisalla frequency
     $     ti(plond,0:plev),            ! interface temperature
     $     nm(plond,plev),              ! midpoint Brunt-Vaisalla frequency
     $     dt,                          ! time step
     $     rdpldv(plond),               ! 1/dp across low level divergence region
     $     tau(plond,-pgwv:pgwv,0:plev), ! wave Reynolds stress
     $     ubi(plond,0:plev),           ! projection of wind at interfaces
     $     ubm(plond,plev),             ! projection of wind at midpoints
     $     xv(plond),                   ! unit vectors of source wind (x)
     $     yv(plond)                    ! unit vectors of source wind (y)
C
C Output variables
C
      real
     $     ut(plond,plev),              ! zonal wind tendency
     $     vt(plond,plev),              ! meridional wind tendency
     $     tau0x(plond),                ! c=0 sfc. stress (zonal)
     $     tau0y(plond)                 ! c=0 sfc. stress (meridional)
C
C Local workspace
C
      integer
     $     i,k,l                        ! loop indexes

      real
     $     d(plond),                    ! "total" diffusivity 
     $     dsat(plond,-pgwv:pgwv),      ! saturation diffusivity
     $     dscal,                       ! fraction of dsat to use
     $     mi,                          ! imaginary part of vertical wavenumber
     $     taudmp,                      ! stress after damping
     $     taumax(plond),               ! max(tau) for any l
     $     tausat(plond,-pgwv:pgwv),    ! saturation stress
     $     ubmc(plond,-pgwv:pgwv),      ! (ub-c)
     $     ubmc2,                       ! (ub-c)**2
     $     ubt(plond,plev),             ! ubar tendency
     $     ubtl,                        ! ubar tendency from wave l
     $     ubtlsat                      ! saturation tendency
C
C Initialize gravity wave drag tendencies to zero
C
      do k=1,plev
         do i=1,plond
            ut(i,k) = 0.
            vt(i,k) = 0.
         end do
      end do
C
C---------------------------------------------------------------------------
C Compute the stress profiles and diffusivities
C---------------------------------------------------------------------------
C
C Loop from bottom to top to get stress profiles      
C
      do k = kbot-1, ktop, -1
C
C Determine the absolute value of the saturation stress and the diffusivity
C for each wave.
C Define critical levels where the sign of (u-c) changes between interfaces.
C
         do i = 1, plon
            d(i) = dback
         end do
         do l = -ngwv, ngwv
            do i = 1, plon
               ubmc(i,l) = ubi(i,k) - c(l)
               tausat(i,l) = abs (efkw * rhoi(i,k) * ubmc(i,l)**3
     $              / (2.*ni(i,k)) )
               if (tausat(i,l) .le. taumin)
     $              tausat(i,l) = 0.0
!-----------------------------------------------------------
! SJL          if (ubmc(i,l) / (ubi(i,k+1) - c(l)) .le. 0.0)
!-----------------------------------------------------------
! The following avoids potentially singularity
               if (ubmc(i,l) * (ubi(i,k+1) - c(l)) .le. 0.0)
     $              tausat(i,l) = 0.0
               dsat(i,l) = (ubmc(i,l) / ni(i,k))**2 *
     $              (efkw * ubmc(i,l)**2 / (rog * ti(i,k) * ni(i,k))
     $              - alpha(k))
               dscal = min (1.0, tau(i,l,k+1) / (tausat(i,l)+taumin))
               d(i) = max( d(i), dscal * dsat(i,l))
            end do
         end do
C
C Compute stress for each wave. The stress at this level is the min of 
C the saturation stress and the stress at the level below reduced by damping.
C The sign of the stress must be the same as at the level below.
C
         do l = -ngwv, ngwv
            do i = 1, plon
               ubmc2 = max(ubmc(i,l)**2, ubmc2mn)
               mi = ni(i,k) / (2. * kwv * ubmc2) *
     $              (alpha(k) + ni(i,k)**2/ubmc2 * d(i))
               taudmp = tau(i,l,k+1)
     $              * exp(-2.*mi*rog*t(i,k+1)*(piln(i,k+1)-piln(i,k)))
               if (taudmp .le. taumin) taudmp = 0.
               tau(i,l,k) = min (taudmp, tausat(i,l))
            end do
         end do
c+
c The orographic stress term does not change across the source region
c Note that k ge ksrcmn cannot occur without an orographic source term
c-
         if (k .ge. ksrcmn) then
            do i = 1, plon
               if (k .ge. ksrc(i)) then
                  tau(i,0,k) = tau(i,0,plev) 
               end if
            end do
         end if
c+
c Require that the orographic term decrease linearly (with pressure) 
c within the low level stress divergence region. This supersedes the above
c requirment of constant stress within the source region.
c Note that k ge kldvmn cannot occur without an orographic source term, since
c kldvmn=plev then and k<=plev-1
c-
         if (k .ge. kldvmn) then
            do i = 1, plon
               if (k .ge. kldv(i)) then
                  tau(i,0,k) = min (tau(i,0,k), tau(i,0,plev)  *
     $                 (1. - fracldv * (pi(i,k)-pi(i,plev)) *
     $                 rdpldv(i)))
               end if
            end do
         end if
c+
c Apply lower bounds to the stress if ngwv > 0.
c-
         if (ngwv .ge. 1) then
c+
c Determine the max value of tau for any l
c-
            do i = 1, plon
               taumax(i) = tau(i,-ngwv,k)
            end do
            do l = -ngwv+1, ngwv
               do i = 1, plon
                  taumax(i) = max(taumax(i), tau(i,l,k))
               end do
            end do
            do i = 1, plon
               taumax(i) = mxrange * taumax(i)
            end do
c+
c Set the min value of tau for each wave to the max of mxrange*taumax or
c mxasym*tau(-c)
c-
            do l = 1, ngwv
               do i = 1, plon
                  tau(i, l,k) = max(tau(i, l,k), taumax(i))
                  tau(i, l,k) = max(tau(i, l,k), mxasym*tau(i,-l,k))
                  tau(i,-l,k) = max(tau(i,-l,k), taumax(i))
                  tau(i,-l,k) = max(tau(i,-l,k), mxasym*tau(i, l,k))
               end do
            end do
            l = 0
            do i = 1, plon
               tau(i,l,k) = max(tau(i,l,k),
     $              mxasym * 0.5 * (tau(i,l-1,k) + tau(i,l+1,k)))
            end do

         end if

      end do
c+
c Put an upper bound on the stress at the top interface to force some stress
c divergence in the top layer. This prevents the easterlies from running
c away in the summer mesosphere, since most of the gravity wave activity
c will pass through a top interface at 75--80 km under strong easterly
c conditions. 
c-
      do l = -ngwv, ngwv
         do i = 1, plon
            tau(i,l,0) = min(tau(i,l,0), 0.5*tau(i,l,1))
         end do
      end do
C
C---------------------------------------------------------------------------
C Compute the tendencies from the stress divergence.
C---------------------------------------------------------------------------
C
C Loop over levels from top to bottom
C
      do k = ktop+1, kbot
C
C Accumulate the mean wind tendency over wavenumber.
C
         do i = 1, plon
            ubt (i,k) = 0.0
         end do
         do l = -ngwv, ngwv
            do i = 1, plon
c+
c Determine the wind tendency including excess stress carried down from above.
c-
               ubtl = g * (tau(i,l,k)-tau(i,l,k-1)) * rdpm(i,k)
c+
c Require that the tendency be no larger than the analytic solution for
c a saturated region [proportional to (u-c)^3].
c-
               ubtlsat = efkw * abs((c(l)-ubm(i,k))**3)
     $              /(2.*rog*t(i,k)*nm(i,k))
c     $              /(2.*7000.     *nm(i,k))
c     $              /(2.*8781.42   *nm(i,k))
               if (ubtl .gt. ubtlsat) ubtl = ubtlsat
c+
c Apply tendency limits to maintain numerical stability.
c 1. du/dt < |c-u|/dt  so u-c cannot change sign (u^n+1 = u^n + du/dt * dt)
c 2. du/dt < tndmax    so that ridicuously large tendency are not permitted
c-
               ubtl = min(ubtl, 0.5 * abs(c(l)-ubm(i,k)) / dt)
               ubtl = min(ubtl, tndmax)
c+
c Accumulate the mean wind tendency over wavenumber.
c-
               ubt (i,k) = ubt (i,k) + sign(ubtl, c(l)-ubm(i,k))
c+
c Redetermine the effective stress on the interface below from the wind 
c tendency. If the wind tendency was limited above, then the new stress
c will be small than the old stress and will cause stress divergence in
c the next layer down. This has the effect of smoothing large stress 
c divergences downward while conserving total stress.
c-
               tau(i,l,k) = tau(i,l,k-1) + ubtl * dpm(i,k) / g
            end do
         end do
C
C Project the mean wind tendency onto the components.
C
         do i = 1, plon
            ut(i,k) = ubt(i,k) * xv(i)
            vt(i,k) = ubt(i,k) * yv(i)
         end do
C
C End of level loop
C
      end do
c+
c-----------------------------------------------------------------------
c Project the c=0 stress in the direction of the source wind for recording
c on the output file.
c-----------------------------------------------------------------------
c-
      do i = 1, plon
         tau0x(i) = tau(i,0,kbot) * xv(i)
         tau0y(i) = tau(i,0,kbot) * yv(i)
      end do
C
      return
      end

