#include <params.h>
      subroutine rdharr(iu      ,len     ,flat    ,flon    ,buf)
C-----------------------------------------------------------------------
C
C Array read for boundary datasets
C
C---------------------------Code history--------------------------------
C
C Original version:  L. Bath
C Standardized:      L. Bath, June 1992
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
      integer iu           ! Logical unit number
      integer len          ! Length of array read
C
C Output arguments
C
      real flat            ! Latitude index
      real flon            ! Number of longitudes
      real buf(len)        ! Array to read into
C
C------------------------------Externals--------------------------------
C
      external endrun 
C
C-----------------------------------------------------------------------
C
      read(iu,end=99) flat,flon,buf
      return
   99 continue
      write(6,*)'RDHARR:  End of File on unit ',iu
      call endrun          ! Abnormal termination
C
      end
 
