#include <params.h>
      subroutine bnddyi(mcdate, mcsec, doy)
C-----------------------------------------------------------------------
C
C Convert date and seconds of day to floating point calendar day, for
C boundary dataset handling
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
C--------------------------Arguments------------------------------------
C
C Input arguments
C
      integer mcdate      ! Current date as yymmdd or yyyymmdd
      integer mcsec       ! Seconds of day for current date
C
C Output arguments
C
      real doy            ! Day of year
C
C--------------------------Local Variables------------------------------
C
      integer mnth        ! Month number
      integer mday        ! Day number of month
      integer jdcon(12)   ! Starting day number for each month
      save jdcon
      data jdcon/0,31,59,90,120,151,181,212,243,273,304,334/
C
C Decode month and day
C
      mnth = mod(mcdate,10000)/100
      if (mnth.lt.1 .or. mnth.gt.12) then
        write(6,*)'BNDDYI: Bad month index=', mnth
        call endrun
      end if
      mday = mod(mcdate,100)
      doy = jdcon(mnth) + mday + mcsec/86400.
      if (doy.lt.1. .or. doy.gt.366.) then
        write(6,*)'BNDDYI: bad day of year = ',doy
        call endrun
      end if
C
      return
      end
 
