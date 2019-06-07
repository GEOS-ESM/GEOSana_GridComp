#include <params.h>
      subroutine trajmp(dt      ,upr     ,vpr     ,phipr   ,lampr   )
C-----------------------------------------------------------------------
C
C Estimate mid-point of parcel trajectory (geodesic coordinates) based
C upon horizontal wind field.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C Reviewed:          D. Williamson, march 1996
C Reviewed:          P. Rasch, april 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
      integer pgls              ! length of latitude slice
      parameter( pgls = plon*plev )
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real dt                   ! time step (seconds)
      real upr(pgls)            ! u-component of wind in local geodesic
      real vpr(pgls)            ! v-component of wind in local geodesic
C
C Input/Output arguments
C
      real phipr(pgls)          ! latitude coord of trajectory mid-point
C
C Output arguments
C
      real lampr(pgls)          ! longitude coord of traj. mid-point
C
C-----------------------------------------------------------------------
C
C  dt      Time interval that corresponds to the parcel trajectory.
C  upr     u-coordinate of velocity corresponding to the most recent
C          estimate of the trajectory mid-point (in geodesic system).
C  vpr     v-coordinate of velocity corresponding to the most recent
C          estimate of the trajectory mid-point (in geodesic system).
C  phipr   Phi value at trajectory mid-point (geodesic coordinates).
C          On entry this is the most recent estimate.
C  lampr   Lambda value at trajectory mid-point (geodesic coordinates).
C
C---------------------------Local variables-----------------------------
C
      integer i                 ! index
C
C-----------------------------------------------------------------------
C
      do i = 1,pgls
         lampr(i) = -.5*dt* upr(i) / cos( phipr(i) )
         phipr(i) = -.5*dt* vpr(i)
      end do
C
      return
      end
 
