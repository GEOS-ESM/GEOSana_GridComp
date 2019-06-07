#include <params.h>
      real function smallest()
C-----------------------------------------------------------------------
C
C This function is used to compute machine epsilon
C
C-------------------------Code History----------------------------------
C
C Original version:  CCM2
C Standardized:      T. Acker, Feb 1996
C Reviewed:          
C 
C------------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
#include <implicit.h>
C--------------------------Local Variables------------------------------
C
      real x    ! dummy var. used in converging on machine epsilon
      real y    ! dummy var. used in converging on machine epsilon
C
C-----------------------------------------------------------------------
C
      x = 1.e0
   10 y = 1.e0 + x
C
      if (y.ne.1.e0) then
         smallest = x
         x = x/2       ! Use interval halving to find machine epsilon
         goto 10
      end if
C
      return
C
      end
 
