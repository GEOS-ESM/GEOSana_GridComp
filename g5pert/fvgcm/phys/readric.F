#include <params.h>
      subroutine readric(kdv     ,pa      ,klen    ,krec    )
C-----------------------------------------------------------------------
C
C Read in a direct-access record
C
C---------------------------Code history--------------------------------
C
C Original version:  L. Bath
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
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer kdv         ! Unit number
      integer klen        ! Length of record
      integer krec        ! Record number
      real    pa(klen)    ! Array
C
C---------------------------Local variables-----------------------------
C
      integer ierr        ! Error status from read statement
      character*44 clabel ! Error message
C
C-----------------------------------------------------------------------
C
      read(kdv,rec=krec,err=900,iostat=ierr) pa
C
      return
C
  900 write (clabel,800) krec,kdv
  800 format('READRIC:Bad read of record ',i4,' from unit ',i2)
      call iostop(ierr,kdv,0,clabel)
C
      end
 
