#include <params.h>
      subroutine copyac(a       ,b       ,length)
C-----------------------------------------------------------------------
C
C Copy data from array a to array b
C
C---------------------------Code history--------------------------------
C
C Original version:  L. Bath, Apr 1992
C Standardized:      L. Bath, Jun 1992
C                    L. Buja, Feb 1996
C
C-----------------------------------------------------------------------
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
      integer length       ! Number of words to copy
      real a(length)       ! Input array name
C
C Output arguments
C
      real b(length)       ! Output array name
C
C---------------------------Local variables-----------------------------
C
      integer i            ! Loop index
C
C-----------------------------------------------------------------------
C
      do i=1,length
         b(i) = a(i)
      end do
C
      return
      end
 
