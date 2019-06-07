#include <preproc.h>
      subroutine msspni

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       msspni.F
* purpose:           set mass store pathnames for history and restart files
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <lsmhis.h>
#include <lsmspmd.h>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer i             !do loop index
      integer lencs         !integer function: length of character string
      character*256 cap     !upper case logid
      character*256 msspn   !mass store path name
      character*1   ctmp    !character temporary
* -----------------------------------------------------------------

      if (masterproc)
     $write(6,*)'initializing mass store history and restart pathnames'

* -----------------------------------------------------------------
* build mass store path names for restart and history files. local history
* and restart files are written to the mass store if MSS is defined. 

*    local file           ->    mass store file
*    locpnh//lsmh_xxxx    ->    /logid/csm/caseid/hist/lsmh_xxxx
*    locpnr//lsmr_xxxx    ->    /logid/csm/caseid/rest/lsmr_xxxx

* logid and caseid are set from the environment variables LOGNAME and CASE.
* logid is converted to upper-case for the mass store
* -----------------------------------------------------------------

      do i = 1, 256
         logid(i:i) = ' '
         caseid(i:i) = ' '
         cap(i:i) = ' '
         msspn(i:i) = ' '
      end do

      if (masterproc) then
C SJL
      	call getenv ('LOGNAME',logid)
       	call getenv ('CASE',caseid)
c gsl   logid = 'userx'
c gsl   caseid ='lsmrun'
C SJL
      end if

      if (masterproc) then
        if (logid(1:1) .eq. ' ') then
           call stoprun
        end if
        if (caseid(1:1) .eq. ' ') then
           call stoprun
        end if
      endif  ! masterproc

      do i = 1, lencs(logid)
         cap(i:i) = logid(i:i)
         ctmp = cap(i:i)
         if (ichar(logid(i:i)).ge.97.and.ichar(logid(i:i)).le.122)
     &       cap(i:i) = char(ichar(ctmp) - 32)
      end do

C SJL
c     msspn = '/'//cap(1:lencs(cap))//'/csm/'//caseid(1:lencs(caseid))
c    &     //'/lnd'
      msspn = '.'
C SJL
      if (msspn(1:1) .eq. ' ') then
         if (masterproc) write (6,*) 'msspni error: mass store file ',
     &        'path not defined'
         call stoprun
      end if

      msspnh = msspn(1:lencs(msspn))//'/hist'
      msspnr = msspn(1:lencs(msspn))//'/rest'

      if (masterproc) write (6,*) 'mass store pathname for history ',
     &     'files = ',
     &     msspnh(1:lencs(msspnh))
      if (masterproc) write (6,*) 'mass store pathname for restart ',
     &     'files = ',
     &     msspnr(1:lencs(msspnr))
      if (masterproc) write (6,*)

      return
      end
 
