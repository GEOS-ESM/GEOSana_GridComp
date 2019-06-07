#include <params.h>
      subroutine lcbas(grd     ,bas1    ,bas2    )
C-----------------------------------------------------------------------
C
C Evaluate the partial Lagrangian cubic basis functions (denominator
C only ) for the grid points and gather grid values
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Olson
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, P. Rasch, August 1992
C Reviewed:          D. Williamson, P. Rasch, March 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
#include <implicit.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real grd(4)               ! grid stencil
C
C Output arguments
C
      real bas1(4)              ! grid values on stencil
      real bas2(4)              ! lagrangian basis functions
C
C---------------------------Local variables-----------------------------
C
      real x0mx1                ! |
      real x0mx2                ! |
      real x0mx3                ! |- grid value differences used in
      real x1mx2                ! |  weights
      real x1mx3                ! |
      real x2mx3                ! |
C
C-----------------------------------------------------------------------
C
      x0mx1   = grd(1) - grd(2)
      x0mx2   = grd(1) - grd(3)
      x0mx3   = grd(1) - grd(4)
      x1mx2   = grd(2) - grd(3)
      x1mx3   = grd(2) - grd(4)
      x2mx3   = grd(3) - grd(4)
C
      bas1(1) = grd(1)
      bas1(2) = grd(2)
      bas1(3) = grd(3)
      bas1(4) = grd(4)
C
      bas2(1) =  1./ ( x0mx1 * x0mx2 * x0mx3 )
      bas2(2) = -1./ ( x0mx1 * x1mx2 * x1mx3 )
      bas2(3) =  1./ ( x0mx2 * x1mx2 * x2mx3 )
      bas2(4) = -1./ ( x0mx3 * x1mx3 * x2mx3 )
C
      return
      end
 
