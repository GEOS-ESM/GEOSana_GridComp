#include <params.h>
      subroutine zint(piln    ,pmln    ,rair    ,gravit  ,tv      ,
     $                zm      ,zi      )
C-----------------------------------------------------------------------
C
C Compute the geopotential height at the interface points from *HEIGHT
C ABOVE THE SURFACE* at midlayer midpoints, using the supplied virtual
C temperatures and pressures.
C
C---------------------------Code history--------------------------------
C
C Original version:  D. Williamson, J. Hack
C Standardized:      L. Buja, Jun 1992, Feb 1996
C Reviewed:          J. Hack, B. Boville, Aug 1992, Apr 1996
C SJL
C Algorithm changed by SJL, 8/24/98
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Commons----------------------------------
#include <pmgrid.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real piln(plond,plevp)    ! Log interface pressures
      real pmln(plond,plev)     ! Log midpoint pressures
      real rair                 ! Gas constant for dry air
      real gravit               ! Acceleration of gravity
      real tv(plond,plev)       ! Virtual temperature
      real zm(plond,plev)       ! Height above surface at midpoints
C
C Output arguments
C
      real zi(plond,plevp)      ! Height above surface at interfaces
C
C---------------------------Local variables-----------------------------
C
      integer i,k               ! Lon, level indices
      real rog                  ! Rair / gravit
C
C-----------------------------------------------------------------------
C
C Add increment to midlayer height
C
 
      rog = rair/gravit
 
C The surface height is zero by definition.
C
      do i=1,plon
        zi(i,plevp) = 0.0
      end do
C
C SJL: finite-volume algorithm
C
      do k=plev,1,-1
        do i=1,plon
          zi(i,k) = zi(i,k+1) + rog*tv(i,k)*(piln(i,k+1) - piln(i,k))
        end do
      end do
      return
      end
