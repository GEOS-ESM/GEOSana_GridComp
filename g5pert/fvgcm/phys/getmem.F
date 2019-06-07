#include <params.h>
#ifndef WORDSIZE
#define WORDSIZE 8
#endif
      subroutine getmem(calledby,nwords,address)
C-----------------------------------------------------------------------
C
C Interface routine for dynamic memory allocation.  Replaces automatic
C arrays on Cray.  Isolates machine-specific code for portable model.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski April, 1993
C                    T. Acker, March 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
#include <implicit.h>
C-----------------------------Arguments---------------------------------
C
C Input 
C
      character*(*) calledby
      integer nwords            ! Number of machine words needed
C
C Output
C
      pointer (address,space)   ! Pointer to dynamically allocated space
      real space
C
C--------------------------Local Variables------------------------------
C
#if ( defined SUN ) || ( defined AIX ) || ( defined linux )
      integer malloc
      external malloc
#endif
#if (defined __sgi)
      integer malloc
#endif
C
C-----------------------------------------------------------------------

#if ( defined AIX )
      address = malloc(%val(WORDSIZE*nwords))
      if (address.eq.0) then
         write(6,*)'GETMEM from ',calledby,' Cannot malloc ',nwords,
     $             ' words'
         call endrun
      end if
#endif

#if ( defined SUN ) || ( defined __sgi ) || ( defined linux )
      address = malloc(WORDSIZE*nwords)
      if (address.eq.0) then
         write(6,*)'GETMEM from ',calledby,' Cannot malloc ',nwords,
     $             ' words'
         call endrun
      end if
#endif
      return
      end

 
