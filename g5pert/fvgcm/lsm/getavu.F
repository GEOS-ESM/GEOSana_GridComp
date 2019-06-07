#include <preproc.h>
      integer function getavu()

#include <vartyp.h>

* ------------------------ code history ---------------------------
* source file:       getavu.F
* purpose:           get next available fortran unit number
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* get next available fortran unit number itst. set lsmiou(itst), in 
* lsmio common block, true. if coupled to ccm, use ccm function navu
* to get available unit number, in which case lsmiou is not needed.
* -----------------------------------------------------------------

* ------------------- common block variables ----------------------
#include <lsmio.h>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
#ifdef COUP_CCM
      integer navu  !ccm function to get next available fortran unit number
      external navu               ! added by SJL
#else
      integer itst  !fortran unit number
#endif
* -----------------------------------------------------------------

#ifdef COUP_CCM
      getavu = navu()
      return
      end
#else
      do itst = 1, 99
         if (.not.lsmiou(itst)) then
            getavu = itst
            lsmiou(itst) = .true.
            return
         end if
      end do
      write (6,*) 'GETAVU: ran out of fortran unit numbers'
      call stoprun
      end
#endif
 
