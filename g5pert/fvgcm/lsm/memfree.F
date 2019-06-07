#include <preproc.h>
      subroutine memfree(address)
C-----------------------------------------------------------------------
C
C Interface routine to free dynamically allocated memory.  
C Replaces automatic arrays on Cray.  Isolates machine-specific code 
C for portable model.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski April, 1993
C                    T. Acker, March 1996
C                    M. Vertenstein, August 1996
C
C-----------------------------------------------------------------------
#include <vartyp.h>
C-----------------------------Arguments---------------------------------
C
C Input
C
      pointer (address,space)   ! Pointer to dynamically allocated space
      real space
C
C--------------------------Local Variables------------------------------
C
#if (defined( CRAY ) || defined( CRAY_T3E ))
      integer ier               ! Error return code
#endif
C
#if ( defined SUN ) || ( defined AIX )
      integer free
      external free
#endif
C-----------------------------------------------------------------------
#if (defined( CRAY ) || defined(  CRAY_T3E ))
CMIC$ GUARD
      call hpdeallc(address,ier,0)
CMIC$ END GUARD
C
C write error message to std. output
C
      if (ier.ne.0) then
         write(6,*)'MEMFREE: Bad return code = ',ier,
     $             ' from hpdeallc: Abort'
         if (ier.eq.-3) then
            write(6,*) address,' is outside heap bounds'
         else if (ier.eq.-4) then
            write(6,*)'Block is already free'
         else if (ier.eq.-5) then
            write(6,*) address,' is not at the beginning of a block'
         else if (ier.eq.-7) then
            write(6,*)'Control word for the next block overwritten'
         end if
         call stoprun
      end if
#endif
#if ( defined AIX )
      address = free(%val(address))
#endif
#if ( defined SUN )
      address = free(address)
#endif
#if ( defined __sgi ) || ( defined linux )
      call free(address)
#endif
C
      return
C
      end


 
