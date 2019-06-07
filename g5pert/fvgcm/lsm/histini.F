#include <preproc.h>
      subroutine histini (nsrest, fnameh)

#include <vartyp.h>
#include <lsmpar.h>
#include <lsmspmd.h>

* ------------------------ code history ---------------------------
* source file:       histini.F
* purpose:           initialize variables for history files
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ netcdf declarations---------------------
#include <netcdf.inc>
* -----------------------------------------------------------------

* ------------------------ input variables ------------------------
      integer nsrest           !type of run
      character*80 fnameh      !full name of lsm history file
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer lencs            !integer function: length of character string
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <dynmem.h>
#include <lsmhis.h>
#include <ncd.h>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer i             !do loop index
* -----------------------------------------------------------------

      if (masterproc) 
     $     write (6,*) 'initializing variables for main history files'

* -----------------------------------------------------------------
* initialize active history fields 
* -----------------------------------------------------------------

      call histlst 

* -----------------------------------------------------------------
* set open file flag to false
* -----------------------------------------------------------------

      ncopnfil = .false.

* -----------------------------------------------------------------
* initial or branch run: set accumulation counters to zero
* -----------------------------------------------------------------

      if (nsrest .eq. 0 .or. nsrest .eq. 3) then
         do i = 1, nslflds
            slfcnt(i) = 0
         end do
         do i = 1, nmlflds
            mlfcnt(i) = 0
         end do
      endif
         
* -----------------------------------------------------------------
* initial run or branch run: set current number of time samples in 
* history file and current history file counter. 
* note, for a continuation run these are part of the restart data. 
* -----------------------------------------------------------------

      if (nsrest .eq. 0 .or. nsrest .eq. 3) then   
         ntim = 0
         nfil = 1
      endif

* -----------------------------------------------------------------
* initial or branch run: set mcdate_t to zero to indicate that the
* current time step is the start of a history interval. for 
* continuation runs, mcdate_t will be read from the restart file.
* -----------------------------------------------------------------

      if (nsrest .eq. 0 .or. nsrest .eq. 3) then
         mcdate_t = 0
      endif

* -----------------------------------------------------------------
* restart run: set flag to obtain time dependent variable ncd id's
* this will only be invoked if current history file is not full 
* (i.e. ntim <=  mfilt), 
* -----------------------------------------------------------------

      if (nsrest .eq. 1) then   
         ncgetvid = .true.
      end if

* -----------------------------------------------------------------
* restart run: if current history file is not full (i.e. ntim <= mfilt) 
* open netcdf history file for write access
* -----------------------------------------------------------------

      if (nsrest .eq. 1 .and. ntim .ne. 0) then
         if (masterproc) then
            call getfil (fnameh, locfnh)
         endif
         if (masterproc) then
            ncid = ncopn(locfnh, NCWRITE, rcode)
            ncopnfil = .true.
            write (6,*) '(HISTINI): successfully opened netcdf ',
     $           ' history file ',locfnh(1:lencs(locfnh))
         end if
      endif

      return
      end
 
