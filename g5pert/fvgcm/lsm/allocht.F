#include <preproc.h>
      subroutine allocht

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history --------------------------------
* source file:       alloctc.F
* purpose:           dynamically allocate memory for single and multi
*                    level fields in lsmhisdyn.h
* date last revised: August 1996
* author:            M. Vertenstein
* ----------------------------------------------------------------------

* ------------------------ common block variables ----------------------
#include <dynmem.h>
#include <lsmhisdyn.h>
* ----------------------------------------------------------------------

* dynamically allocate memory for lpt and kpt dependent variables 
* in lsmhisdyn.h

      call memget ('allocht' ,pslfval,kpt*mslflds    )           
      call memget ('allocht' ,pmlfval,kpt*mmlflds*msl)          

      return
      end




 
