#include <preproc.h>
      subroutine relavu (iunit)

#include <vartyp.h>

* ------------------------ code history ---------------------------
* source file:       relavu.F
* purpose:           close and release fortran unit no longer in use
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* close and release fortran unit number iunit. set lsmiou(iunit), in
* lsmio common block, false. if coupled to ccm, use ccm function 
* relunit to close/release unit number.
* -----------------------------------------------------------------

* ------------------- common block variables ----------------------
#include <lsmio.h>
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
      integer iunit    !fortran unit number
* -----------------------------------------------------------------

#ifdef COUP_CCM
      call relunit(iunit)
#else
      if (.not.lsmiou(iunit)) then
         write (6,*) 'RELAVU: unit ',iunit,' is not flagged as in use'
         call stoprun
      end if
      if (iunit.lt.1 .or. iunit.gt.99) then
         write (6,*) 'RELAVU: attempt to return out of range unit'
         call stoprun
      end if
      close(iunit)
      lsmiou(iunit) = .false.
#endif
      return
      end
 
