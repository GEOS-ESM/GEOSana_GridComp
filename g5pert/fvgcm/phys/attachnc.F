#include <params.h>
      subroutine attachnc(ncid   ,mspath  )
C-----------------------------------------------------------------------
C
C Attach netcdf boundary dataset to job
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski Mar 1997
C
C-----------------------------------------------------------------------
c
c $Id$
c
#include <implicit.h>
C-----------------------------------------------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <netcdf.inc>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer ncid           ! netcdf file id
      character mspath*(*)   ! Full Mass Store pathname
C
C---------------------------Local variables-----------------------------
C
      logical   lexist       ! Check for file existence
      integer   klen         ! Length of Mass Store pathname
      integer   i            ! Loop index
      character text*130     ! Printed as msread command
      character locdn*64     ! Local file name
      integer   ierr         ! Return code
C
C------------------------------Externals--------------------------------
C
      integer lenchr         ! Function: Find length of character string
      external lenchr
C
C-----------------------------------------------------------------------
C
C Build local dataset name
C
      klen = lenchr(mspath)
      do i=klen,1,-1
        if (mspath(i:i).eq.'/') go to 100
      end do
      i = 0
  100 continue
      if (klen-i.gt.64) then
        write(6,*)'ATTACHNC:Local dataset name ',mspath(i+1:klen),
     $            ' exceeds 64 characters'
        call endrun
      end if
      locdn = mspath(i+1:klen)
!
C First check if file is in current directory.
C If its not there then acquire from mass store.
C
      inquire(file=locdn,exist=lexist)
      if (.not.lexist) then
C
C Acquire the file from mass store
C
        text = 'ln -s '//mspath(1:lenchr(mspath))
     .         //' '//locdn(1:lenchr(locdn))
        write(6,*)'ATTACHNC: linking data file as follows:',text
        call system(text(1:lenchr(text)))
        ierr = 0
        if (ierr.eq.0) then
          write(6,*)' File ',locdn,' read from MSS'
        else
          write(6,*)'ATTACHNC: failed cmd=',text
          call endrun
        end if
      end if
C
C Open the file
C
      ncid = NCOPN(locdn, NCNOWRIT, ierr)
      if (ierr.eq.0) then
        write(6,*)'ATTACHNC: NCOPN returns id ',ncid,' for file ',locdn
      else
        write(6,*)'ATTACHNC: Failed to open file ',locdn
        call endrun
      end if
C
      return
      end
 
