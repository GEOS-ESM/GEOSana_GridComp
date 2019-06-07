#include <preproc.h>
      subroutine histhan (nstep, nestep, beglatspmd, endlatspmd)

#include <vartyp.h>
#include <lsmpar.h>
#include <lsmspmd.h>
#include <basinpar.h>

* ------------------------ code history ---------------------------
* source file:       histhan.F
* purpose:           main history and restart file handler
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* this code does the following for every time step:
*   o increments field accumulation counters

* this code does the following at the end of a history interval:
*   o increments the current time sample counter: ntim <= mfilt
*   o opens a new history file if needed (i.e., when ntim = 1)
*   o wites history data to current history file
*   o resets field accumulation counters to zero

* this code does the following when the history file is full
* (i.e., ntim = mfilt) or if the last time step of the simulation
* is a history write:
*   o close history file and dispose to mass store
*   o reset ntim time sample counter to zero (only if file is full)
*   o increment nfil file counter by one (only if file is full)
*   o write restart files and dispose to mass store

* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
      integer nstep         !current time index
      integer nestep        !ending time index
      integer beglatspmd    !beg lat on proc, degenerates to 1 for nonspmd
      integer endlatspmd    !end lat on proc, degenerates to lsmlat for nonspmd
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <dynmem.h>
#include <lsmhis.h>
#include <lsmtvdyn.h>
#include <basin.h>
#include <ncd.h>
* -----------------------------------------------------------------

* ------------------------ netcdf declarations---------------------
#include <netcdf.inc>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      real frac_t           !fractional day at start of current time sample
      real frac_f           !fractional day for first time sample in file
      real frac             !current fractional day 
      real hour_t           !date fractional hour at start of current time sample
      real hour_f           !date fractional hour for first time sample in file
      real hour             !date current fractional hour 
      integer i             !loop index
      integer lencs         !integer function: length of character string
      character*80 msscom   !mass store comment
* -----------------------------------------------------------------
* skip nstep=0 if monthly average
* -----------------------------------------------------------------
      if (nstep.eq.0 .and. (ninavg.eq.'q' .or. ninavg.eq.'Q')) return

* -----------------------------------------------------------------
* increment counters
* -----------------------------------------------------------------

      do i = 1, nslflds
         if (slftyp(i) .eq. naver) then   !time average
            slfcnt(i) = slfcnt(i) + 1
         else                             !instantaneous, maximum, or minimun
            slfcnt(i) = 1
         end if
      end do

      do i = 1, nmlflds
         if (mlftyp(i) .eq. naver) then   !time average
            mlfcnt(i) = mlfcnt(i) + 1
         else                             !instantaneous, maximum, or minimun
            mlfcnt(i) = 1
         end if
      end do

      bascnt = bascnt + 1                 ! always assume time avg

* -----------------------------------------------------------------
* end of history interval?
* -----------------------------------------------------------------

      if (ehi) then

* increment current time sample counter

         ntim = ntim + 1

* create new netcdf history file
      
#if ( !defined LSMH_off)
         if (ntim .eq. 1) then
            call histnam 
	    if (masterproc) call histcrt 
            ncopnfil = .true.
         endif
#endif

* build time comment for current time sample based on
*   o calendar info at start of time interval 
*     (mcdate_t, mcsec_t, mdcur_t, mscur_t)
*   o current calendar info 
*     (mcdate  , mcsec  ,mdcur   , mscur  )

         frac_t = float(mscur_t)/86400.
         frac   = float(mscur  )/86400.
         hour_t = float(mcsec_t)/3600.
         hour   = float(mcsec  )/3600.
         if (masterproc) write (timcom,202) mdcur_t,frac_t,mdcur,frac,
     &                      hour_t,mcdate_t,hour,mcdate

* set mcdate_t to zero to indicate next time step is start of time sample

         mcdate_t = 0

* write data to netcdf history file

#if ( !defined LSMH_off)
         call histwrt (nstep,beglatspmd)
#endif

* if not monthly average, then write out beginning and end of time
* sample interval

         if ( masterproc ) then
         if (ninavg.ne.'q' .and. ninavg.ne.'Q') then
            write (6,*) '---------------------------------------'
            write (6,*) '(HISTHAN): nstep = ',nstep
            write (6,*) '(HISTHAN): wrote current time sample to ',
     $           'local history file ',locfnh(1:lencs(locfnh))
            write (6,*) '(HISTHAN): ',timcom(1:lencs(timcom))
            write (6,*) '---------------------------------------'
            write (6,*)
         endif
         endif


* reset history field accumulation counters to zero

         do i = 1, nslflds
            slfcnt(i) = 0
         end do
         do i = 1, nmlflds
            mlfcnt(i) = 0
         end do

         bascnt = 0

      endif


* -----------------------------------------------------------------
* close-dispose history file and write-dispose restart file when:
*   o history file is full (ntim = mfilt) or 
*   o have reaced end of run as specified by either ccm3 or the
*     flux coupler 
* -----------------------------------------------------------------

      if (ntim .eq. mfilt ) go to 100
      if (nstep.eq.nestep ) go to 100
      return
100   continue

* -----------------------------------------------------------------
* if current time step is end of history interval then build 
* mass store comment for first and last time samples in file using
*   o calendar info when first time sample written 
*     (mcdate_f, mdcur_f, mscur_f)
*   o current calendar info 
*     (mcdate, mdcur, mscur)
* -----------------------------------------------------------------

      if (masterproc) then
         if (ehi) then
            frac_f = float(mscur_f)/86400.
            hour_f = float(mcsec_f)/3600.
            if (ninavg.eq.'q' .or. ninavg.eq.'Q') then
               write (msscom,201) mdcur_f,frac_f,mdcur,frac,
     &              hour_f,mcdate_f,hour,mcdate
  201          format ('MONTH MEAN DAYS: ',i6.6,f4.3,'-',i6.6,f4.3,
     &              ' DATES:',f6.3,'Z ',i8.8,' -',f6.3,'Z ',i8.8) 
            else
               write (msscom,202) mdcur_f,frac_f,mdcur,frac,
     &              hour_f,mcdate_f,hour,mcdate
  202          format ('DAYS: ',i6.6,f4.3,'-',i6.6,f4.3,
     &              ' DATES:',f6.3,'Z ',i8.8,' -',f6.3,'Z ',i8.8) 
            end if
         endif
      end if

* -----------------------------------------------------------------
* if appropriate, close netcdf history file and dispose to mass store
* -----------------------------------------------------------------

* if history file is open then close it and dispose to mass store.
* note that ncopnfil is then set to false in subroutine restwrt.

      if (masterproc) then
         if (ncopnfil) then
#if ( !defined LSMH_off )
            call ncclos(ncid, rcode)
#endif
            if (masterproc) then
               write (6,*) '---------------------------------------'
               write (6,*) '(HISTHAN): nstep = ',nstep
               write (6,*) '(HISTHAN): closing local netcdf hist file ',
     &              locfnh(1:lencs(locfnh))
#if ( !defined LSMH_off )
               call putfil (locfnh, msspnh, irt, nswrps, msscom)
#endif
               if (ehi) write(6,*) '(HISTHAN): ',msscom(1:lencs(msscom))
            endif
         endif
      endif

* -----------------------------------------------------------------
* build restart file name. 
* set ntim and nfil counters to correct values for restart file. 
* write and dispose restart files.
* -----------------------------------------------------------------

      if( masterproc ) call restnam 

      if (ntim .eq. mfilt)  then
         ntim = 0              !reset number of time samples to 0
         nfil = nfil + 1       !increment file number
      end if

      call restwrt (nstep, beglatspmd, endlatspmd)

* -----------------------------------------------------------------
* end echoing
* -----------------------------------------------------------------
      if (masterproc) then
         write (6,*) '---------------------------------------'
         write (6,*)
      end if
      return
      end
 
