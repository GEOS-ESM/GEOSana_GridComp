#include <params.h>
      subroutine writeric(kdv     ,pa      ,klen    ,krec)
C-----------------------------------------------------------------------
C
C Write out a direct-access record
C
C----------------------------Code History-------------------------------
C
C Original version:  J. Rosinski
C Standardized:      J. Rosinski, June 1992
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
      integer kdv              ! unit number
      integer klen             ! length of record
      integer krec             ! record number
      real pa(klen)            ! array
C
C---------------------------Local variables-----------------------------
C
      integer ierr             ! error status from read statement
      character*46 clabel      ! error message
C
C-----------------------------------------------------------------------
C
      write(kdv,rec=krec,err=900,iostat=ierr) pa
C
      return
C
  900 write (clabel,800) krec,kdv
  800 format('WRITERIC:Bad write of record ',i4,' from unit ',i2)
      call iostop(ierr,kdv,0,clabel)
C
      end







 
