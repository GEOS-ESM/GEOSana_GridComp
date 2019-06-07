#include <preproc.h>
      subroutine allocbas

#include <vartyp.h>

* ------------------------ code history --------------------------------
* source file:       alloctc.F
* purpose:           dynamically allocate memory for basindyn.h
* date last revised: April 1997
* author:            M. Vertenstein
* ----------------------------------------------------------------------

* ------------------------ common block variables ----------------------
#include <dynmem.h>
#include <basinpar.h>
#include <basindyn.h>
* ----------------------------------------------------------------------

* dynamically allocate memory for lpt and kpt dependent variables 
* in basindyn.h

      call memget ('alloctc' ,pdrnbasn,lpt*ndrnmax)
      call memget ('alloctc' ,pdrnfrac,lpt*ndrnmax)
      call memget ('alloctc' ,pdrnarea,lpt        )
      call memget ('alloctc' ,prunoff ,lpt        )

      return
      end



 
