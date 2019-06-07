#include <preproc.h>
      subroutine setarr (pa, kdim, pvalue)

#include <vartyp.h>

* ------------------------ code history ---------------------------
* source file:       setarr.F
* purpose:           set array pa(kdim) to pvalue
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
* input
      integer kdim      !dimension of array pa
      real pvalue       !value to store in pa

* output 
      real pa(kdim)     !array to set
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer j         !loop index
* -----------------------------------------------------------------
 
      do j = 1, kdim
         pa(j) = pvalue
      end do
      return
      end
 
