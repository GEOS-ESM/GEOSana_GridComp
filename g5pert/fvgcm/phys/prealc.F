#include <params.h>
      subroutine prealc(fname   ,nwords  )
C-----------------------------------------------------------------------
C
C Routine to preallocate ksiz words (converted to blocks and rounded up)
C for file fname.  Called by CCM3 for history tape and regeneration
C datasets.
C
C-------------------------Code History----------------------------------
C
C Original version:  L. Bath, May 1993 using system routine fall  
C Standardized:      L. Buja, Feb 1996
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
      character*(*) fname   ! Character variable containing local
C                           !  pathname of file to preallocate
      integer nwords        ! Expected total file size in Cray words
#if ( defined CRAY )
C
C---------------------------Local variables-----------------------------
C
      integer nbyt          ! File size (bytes)
      integer istat         ! Status return from ishell
      character*10 cbyt     ! Character form of nbyt
      character*10 cbytlj   ! Left-justified cbyt
      character*130 cmd     ! String to build setf command
C
C------------------------------Externals--------------------------------
C
      integer lenchr, ishell
      external lenchr       ! Actual size of character string
      external ishell       ! System routine, execute shell command
C
C-----------------------------------------------------------------------
C
C Add 2% fudge factor for COS blocking
C
      nbyt = nint(1.02*(nwords*8))
      write(cbyt,'(i10)') nbyt
      call ljust(cbyt    ,cbytlj  )
C
      cmd = 'setf -n '//cbytlj(1:lenchr(cbytlj))//' '//
     $      fname(1:lenchr(fname))
      print *,' sending string to ishell from prealc:'
      print *,cmd
      istat = ishell(cmd)
      if (istat.ne.0) then
        print *,'PREALC: Cannot allocate ',nbyt,' bytes for file'
        print *,fname,', status word = ',istat
      endif
#else
      write(6,*)'PREALC: Dummy routine for this machine architecture'
#endif
C
      return
      end

 
