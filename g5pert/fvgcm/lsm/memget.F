#include <preproc.h>
#ifndef WORDSIZE
#define WORDSIZE 8
#endif
      subroutine memget(calledby,address,nwords)
C-----------------------------------------------------------------------
C
C Interface routine for dynamic memory allocation.  Replaces automatic
C arrays on Cray.  Isolates machine-specific code for portable model.
C
C---------------------------Code history--------------------------------
C
C Original version:  M. Vertenstein, August 1996 
C
C-----------------------------------------------------------------------
#include <vartyp.h>
C-----------------------------Arguments---------------------------------

      character*(*) calledby    ! Calling routine
      integer nwords            ! Number of machine words needed
      pointer (address,space)   ! Pointer to dynamically allocated space

#if ( defined linux )
      integer malloc
#endif

C
C--------------------------Local Variables------------------------------
C
      real space
#if (defined( CRAY ) || defined( CRAY_T3E ))
      integer ier               ! Error return code
#endif
#if ( defined SUN ) || ( defined AIX )
      integer malloc
      external malloc
#endif
C
C-----------------------------------------------------------------------
C
#if (defined( CRAY ) || defined( CRAY_T3E ))
CMIC$ GUARD
      call hpalloc(address,nwords,ier,0)
CMIC$ END GUARD
      if (ier.ne.0) then
         write(6,*)'MEMGET from ',calledby,': Bad return code = ',ier
         if (ier.eq.-1) then
            write(6,*) nwords,' is not between 1 and 2**31'
         else if (ier.eq.-2) then
            write(6,*)'No more memory available'
         end if
         call stoprun
      end if
#endif
#if ( defined AIX )
      address = malloc(%val(WORDSIZE*nwords))
      if (address.eq.0) then
         write(6,*)'MEMGET from ',calledby,' Cannot malloc ',nwords,
     $             ' words'
         call stoprun
      end if
#endif
#if ( defined SUN ) || ( defined __sgi ) || ( defined linux )
      address = malloc(WORDSIZE*nwords)
      if (address.eq.0) then
         write(6,*)'MEMGET from ',calledby,' Cannot malloc ',nwords,
     $             ' words'
         call stoprun
      end if
#endif
      return
      end

 
