#include <params.h>
      subroutine lcdbas(grd     ,dbas2   ,dbas3   )
C-----------------------------------------------------------------------
C
C Calculate weights used to evaluate derivative estimates at the
C inner grid points of a four point stencil based on Lagrange
C cubic polynomial through four unequally spaced points.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C Reviewed:          D. Williamson, P. Rasch, March 1996
C
c
c $Id$
c $Author$
c

C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real grd(4)               ! grid stencil
C
C Output arguments
C
      real dbas2(4)             ! derivatives at grid point 2.
      real dbas3(4)             ! derivatives at grid point 3.
C
C-----------------------------------------------------------------------
C
C  grd    Coordinate values of four points in stencil.
C  dbas2  Derivatives of the four basis functions at grid point 2.
C  dbas3  Derivatives of the four basis functions at grid point 3.
C
C---------------------------Local variables-----------------------------
C
      real x1                   !  |
      real x2                   !  |- grid values
      real x3                   !  |
      real x4                   !  |
      real x1mx2                !  |
      real x1mx3                !  |
      real x1mx4                !  |- differences of grid values
      real x2mx3                !  |
      real x2mx4                !  |
      real x3mx4                !  |
C
C-----------------------------------------------------------------------
C
      x1 = grd(1)
      x2 = grd(2)
      x3 = grd(3)
      x4 = grd(4)
      x1mx2 = x1 - x2
      x1mx3 = x1 - x3
      x1mx4 = x1 - x4
      x2mx3 = x2 - x3
      x2mx4 = x2 - x4
      x3mx4 = x3 - x4
C
      dbas2(1) =   x2mx3 * x2mx4 / ( x1mx2 * x1mx3 * x1mx4 )
      dbas2(2) =   -1./x1mx2 + 1./x2mx3 + 1./x2mx4
      dbas2(3) = - x1mx2 * x2mx4 / ( x1mx3 * x2mx3 * x3mx4 )
      dbas2(4) =   x1mx2 * x2mx3 / ( x1mx4 * x2mx4 * x3mx4 )
C
      dbas3(1) = - x2mx3 * x3mx4 / ( x1mx2 * x1mx3 * x1mx4 )
      dbas3(2) =   x1mx3 * x3mx4 / ( x1mx2 * x2mx3 * x2mx4 )
      dbas3(3) =   -1./x1mx3 - 1./x2mx3 + 1./x3mx4
      dbas3(4) = - x1mx3 * x2mx3 / ( x1mx4 * x2mx4 * x3mx4 )
C
      return
      end
 
