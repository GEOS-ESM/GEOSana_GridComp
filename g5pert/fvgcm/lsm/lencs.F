#include <preproc.h>
      integer function lencs (chrstr)

#include <vartyp.h>

* ------------------------ code history ---------------------------
* source file:       lencs.F
* purpose:           return position of right-most non-blank, non-null 
*                    character in chrstr
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
* input
      character*(*) chrstr       !input character string
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer l
* -----------------------------------------------------------------

#if (defined linux)
      lencs = 1
#else
      lencs = 0
#endif
      do l = len(chrstr),1,-1
         if (chrstr(l:l).ne.' ' .and. chrstr(l:l).ne.char(0)) then
            lencs = l
            goto 10
         end if
      end do
   10 return
      end
 
