#include <params.h>
      subroutine mgwprof (u, v, t, pm, pi, rhoi, ni, ti, nm)
C-----------------------------------------------------------------------
C
C Compute profiles of background state quantities for the multiple
c gravity wave drag parameterization.
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
      real 
     $     u(plond,plev),               ! midpoint zonal wind
     $     v(plond,plev),               ! midpoint meridional wind
     $     t(plond,plev),               ! midpoint temperatures
     $     pm(plond,plev),              ! midpoint pressures
     $     pi(plond,0:plev)             ! interface pressures
C
C Output variables
C
      real
     $     rhoi(plond,0:plev),          ! interface density
     $     ni(plond,0:plev),            ! interface Brunt-Vaisalla frequency
     $     ti(plond,0:plev),            ! interface temperature
     $     nm(plond,plev)               ! midpoint Brunt-Vaisalla frequency
C
C Local workspace
C
      integer
     $     i,k                          ! loop indexes

      real
     $     dtdp,
     $     n2                           ! Brunt-Vaisalla frequency squared
c
C-----------------------------------------------------------------------------
C Determine the interface densities and Brunt-Vaisala frequencies.
C-----------------------------------------------------------------------------
C
C The top interface values are calculated assuming an isothermal atmosphere 
C above the top level.
C
      k = 0
      do i = 1, plon
         ti(i,k) = t(i,k+1)
         rhoi(i,k) = pi(i,k) / (r*ti(i,k))
         ni(i,k) = sqrt (g*g / (cp*ti(i,k)))
      end do
c+
c Interior points use centered differences
c-
      do k = 1, plev-1
         do i = 1, plon
            ti(i,k) = 0.5 * (t(i,k) + t(i,k+1))
            rhoi(i,k) = pi(i,k) / (r*ti(i,k))
            dtdp = (t(i,k+1)-t(i,k)) / (pm(i,k+1)-pm(i,k))
            n2 = g*g/ti(i,k) * (1./cp - rhoi(i,k)*dtdp)
            ni(i,k) = sqrt (max(n2min, n2))
         end do
      end do
c+
c Bottom interface uses bottom level temperature, density; next interface
c B-V frequency.
c-
      k = plev
      do i = 1, plon
         ti(i,k) = t(i,k)
         rhoi(i,k) = pi(i,k) / (r*ti(i,k))
         ni(i,k) = ni(i,k-1)
      end do
C
C-----------------------------------------------------------------------------
C Determine the midpoint Brunt-Vaisala frequencies.
C-----------------------------------------------------------------------------
C
      do k=1,plev
         do i=1,plon
            nm(i,k) = 0.5 * (ni(i,k-1) + ni(i,k))
         end do
      end do
C
      return
      end

