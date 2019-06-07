#include <params.h>
      subroutine rdhdr(kunit   ,lenhi   ,lenhis  ,lenhc   ,lenhcs  ,
     $                 lenhr   ,hedi    ,hedc    ,hedr    ,mhflds  ,
     $                 mchflds ,kerr    )
C-----------------------------------------------------------------------
C
C Read the history file header.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      L. Bath, June 1992
C                    T. Acker, March 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
#include <pmgrid.h>
C-----------------------------------------------------------------------
#include <pagrid.h>
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer kunit           ! History file unit
      integer lenhi           ! Declared length of integer header hedi
      integer lenhis          ! Length of integer scalar portion
      integer lenhc           ! Declared length of character header hedc
      integer lenhcs          ! Length of character scalar portion
      integer lenhr           ! Declared length of real header hedr
C
C Output arguments
C
      integer hedi(*)         ! Array to store integer header record
      character*8 hedc(*)     ! Character header record
      real hedr(*)            ! Real header record
      integer mhflds(*)       ! Kunit portion of mflds array
      character*8 mchflds(*)  ! Kunit portion of mcflds array
      integer kerr            ! Error return
C                                
C---------------------------Local variables-----------------------------
C                                
      integer len             ! Length of integer header record
      integer jend            ! Length of vectors in record
      integer j, i            ! Indices
C
C------------------------------Externals--------------------------------
C
      external endrun         ! Abnormal termination
C
C-----------------------------------------------------------------------
C
      kerr = 0
C
C Check header lengths and then read header
C
      read(kunit,end=99)len
      if (len.gt.lenhi) then
         write(6,*)'RDHDR: integer header record too long.'
         write(6,*)'Allocated space = ',lenhi,', Actual length = ',len
         call endrun
      end if
      backspace(kunit)
C
C Calculate buffer length for integer record
C
      jend = len - lenhis
      read(kunit,end=99) (hedi(i),i=1,lenhis), (mhflds(j),j=1,jend)
      if (hedi(31).gt.lenhc) then
         write(6,*) 'RDHDR: character header record too long.'
         write(6,*) 'Allocated space = ',lenhc,', Actual length = ',
     $                hedi(31)
         call endrun
      else if (hedi(32).gt.lenhr) then
         write(6,*) 'RDHDR: real header record too long.'
         write(6,*) 'Allocated space = ',lenhr,', Actual length = ',
     $               hedi(32)
         call endrun
      end if
C
C Calculate buffer length for character record
C
      jend = 2*hedi(16)
      read(kunit,end=991) (hedc(i),i=1,lenhcs), (mchflds(j),j=1,jend)
      read(kunit,end=992) (hedr(i),i=1,hedi(32))
C
      return
   99 kerr = 1
      return
  991 kerr = 2
      return
  992 kerr = 3
      return
C
      end


 
