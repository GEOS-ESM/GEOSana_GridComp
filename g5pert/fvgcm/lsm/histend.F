#include <preproc.h>
      subroutine histend (nstep, xmcdate, xmcsec)

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       histend.F
* purpose:           determine if end of history interval
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
      integer nstep       !time index
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <dynmem.h>
#include <lsmhis.h>
#include <lsmctl.h>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
* mdcur, mscur, mcdate, mcsec: commented because in the lsmhis common block
c     integer mdcur       !nstep: current day (0, ...)
c     integer mscur       !nstep: current seconds of current day (0 -> 86400)
c     integer mcdate      !nstep: current date (yyyymmdd format) (e.g., 00030131)
c     integer mcsec       !nstep: current seconds of current date (0 -> 86400)
      integer kda         !nstep: day (1 -> 31)
      integer kmo         !nstep: month (1 -> 12)
      integer kyr         !nstep: year (0 -> ...)
      real calday         !nstep: calendar day at greenwich (1.00 -> 365.99)
      integer mdcurm1     !mdcur for nstep-1
      integer mscurm1     !mscur for nstep-1
      integer mcdatem1    !mcdate for nstep-1
      integer mcsecm1     !mcsec for nstep-1
      integer kmom1       !kmo for nstep-1
      integer kyrm1       !kmo for nstep-1
      real caldaym1       !calday for nstep-1

! SJL
      integer xmcdate, xmcsec
      integer nhms, nhmsf, n
      nhmsf(n)  = n/3600*10000 + mod(n,3600 )/ 60*100 + mod(n, 60)

* -----------------------------------------------------------------

* -----------------------------------------------------------------
* skip nstep=0 if monthly average
* -----------------------------------------------------------------

      if (nstep.eq.0 .and. (ninavg.eq.'q' .or. ninavg.eq.'Q')) return

* -----------------------------------------------------------------
* set calendar for current time step
* -----------------------------------------------------------------

! SJL
!     call calendr(nstep  ,dtlsm  ,mdbase  ,msbase  ,mbdate  ,
!    &             mbsec  ,mdcur  ,mscur   ,mcdate  ,mcsec   ,
!    &             calday )

      mcdate = xmcdate
      mcsec  = xmcsec

      mscur = mcsec
      nhms = nhmsf( mcsec )

      call mcalday(mcdate, nhms, calday)
      mdcur = mcdate - mdbase     ! quick & dirty
                                  ! Does not give the correct value
                                  ! But it is not used if lsm history file
                                  ! is turned off
! SJL

      kda = mod(mcdate,100)
      kmo = mod(mcdate,10000)/100
      kyr = mcdate/10000

* -----------------------------------------------------------------
* set calendar for current time step - 1
* -----------------------------------------------------------------

! SJL
!     call calendr(nstep-1,dtlsm  ,mdbase  ,msbase  ,mbdate  ,
!    &             mbsec  ,mdcurm1,mscurm1 ,mcdatem1,mcsecm1 ,
!    &             caldaym1)

      mcdatem1 = mcdate
      call tick (mcdatem1, nhms, nint(-dtlsm) )
! SJL
      kmom1 = mod(mcdatem1,10000)/100
      kyrm1 = mcdatem1/10000

* -----------------------------------------------------------------
* determine if end of history interval

* the lsm calendar is based on the ccm calendar. note the
* peculiarities of this calendar. using a time step of 1200 
* sec and with a base date of 000901:

* o daily-averaged data for the first day in september (nstep = 1 
*   to 72) are written on mcdate = 00/09/02 with mscur = 0

*   nstep  mcdate kyr kmo kda mdcur mscur mcsec calday
*      1   000901   0   9   1     0  1200  1200 244.014
*      .      .     .   .   .     .    .     .     .
*     72   000902   0   9   2     1    0     0 245.000

*   in general: daily-averaged data for the first day in
*   month mm are written on mcdate = yyyy/mm/02 with mscur = 0

* o daily-averaged data for the 30th day (last day in september)
*   (nstep = 2089 to 2160) are written on mcdate = 00/10/01 mscur = 0

*   nstep   mcdate kyr kmo kda mdcur mscur mcsec calday
*   2089    000930   0   9  30    29  1200  1200 273.014
*     .       .      .   .   .     .   .     .      .
*   2159    000930   0   9  30    29 85200 85200 273.986
*   2160    001001   0  10   1    30     0     0 274.000

*   in general: daily-averaged data for the last day in
*   month mm are written on mcdate = yyyy/mm+1/01 with mscur = 0

* o the last time step in the month is defined by kmo(nstep-1) 
*   ne kmo(nstep). the appropriate month and year are based
*   on nstep - 1. e.g., last time step in year 0 is 122*72=8784

*   nstep   mcdate kyr kmo kda mdcur mscur mcsec calday
*   8783    001231   0  12  31   121 85200 85200 365.986
*   8784    010101   1   1   1   122     0     0   1.000
* -----------------------------------------------------------------

      if (ninavg.eq.'q' .or. ninavg.eq.'Q') then      !monthly average
         if (kmo .ne. kmom1) then
            ehi = .true.
         else
            ehi = .false.
         end if
      else
         if (mod(nstep,nhtfrq) .eq. 0) then
            ehi = .true.
         else
            ehi = .false.
         end if
      end if

* -----------------------------------------------------------------
* calendar info for current time sample: start of time interval
* -----------------------------------------------------------------

      if (mcdate_t .eq. 0) then
         mcdate_t = mcdate
         mcsec_t  = mcsec
         mdcur_t  = mdcur
         mscur_t  = mscur
      end if

* -----------------------------------------------------------------
* calendar info for history file: time first time sample is written
* if monthly-average: mimic daily-averaged data, where first day
* in month mm is written on mcdate = yyyy/mm/02 and mscur = 0
* -----------------------------------------------------------------

      if (ninavg.eq.'q' .or. ninavg.eq.'Q') then      !monthly average
         if (kda.eq.2 .and. mscur.eq.0) then
            nmon = kmo                                !for "nmon-nyr" hist file
            nyr = kyr                                 !for "nmon-nyr" hist file
            mcdate_f = mcdate
            mcsec_f  = mcsec
            mdcur_f  = mdcur
            mscur_f  = mscur
         end if
      else
         if (ehi .and. ntim .eq. 0) then
            nmon = kmo                                !not used
            nyr = kyr                                 !not used
            mcdate_f = mcdate
            mcsec_f  = mcsec
            mdcur_f  = mdcur
            mscur_f  = mscur
         end if
      end if

* -----------------------------------------------------------------
* note that when history file is full, mcdate, mdcur and mscur 
* will be the ending time information for the last time sample in 
* the file. it will also indicate the end of the current time sample.
* -----------------------------------------------------------------

      return
      end
 
