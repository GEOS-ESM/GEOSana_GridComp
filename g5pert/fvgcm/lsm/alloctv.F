#include <preproc.h>
      subroutine alloctv

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history --------------------------------
* source file:       alloctc.F
* purpose:           dynamically allocate memory for lsmtvdyn.h
* date last revised: August 1996
* author:            M. Vertenstein
* ----------------------------------------------------------------------

* ------------------------ common block variables ----------------------
#include <dynmem.h>
#include <lsmtvdyn.h>
* ----------------------------------------------------------------------

* dynamically allocate memory for lpt and kpt dependent variables 
* in lsmtvdyn.h

      call memget ('alloctv' ,ph2osno,kpt      )
      call memget ('alloctv' ,ph2ocan,kpt      )
      call memget ('alloctv' ,ph2osoi,kpt*msl  )
      call memget ('alloctv' ,ptv    ,kpt      )
      call memget ('alloctv' ,ptg    ,kpt      )
      call memget ('alloctv' ,ptsoi  ,kpt*msl  )
      call memget ('alloctv' ,pmoz   ,kpt      )
      call memget ('alloctv' ,peah   ,kpt      )
      call memget ('alloctv' ,psoot  ,kpt      )
      call memget ('alloctv' ,phsno  ,kpt      )
      call memget ('alloctv' ,pfsno  ,kpt      )
      call memget ('alloctv' ,pfwet  ,kpt      )
      call memget ('alloctv' ,phtop  ,kpt      )
      call memget ('alloctv' ,ptlai  ,kpt      )
      call memget ('alloctv' ,ptsai  ,kpt      )
      call memget ('alloctv' ,pelai  ,kpt      )
      call memget ('alloctv' ,pesai  ,kpt      )
      call memget ('alloctv' ,pfoln  ,kpt      )
      call memget ('alloctv' ,pstemb ,kpt      )
      call memget ('alloctv' ,prootb ,kpt      )
      call memget ('alloctv' ,psoilc ,kpt      )
      call memget ('alloctv' ,pigs   ,kpt      )
      call memget ('alloctv' ,palbd  ,kpt*mband)
      call memget ('alloctv' ,palbi  ,kpt*mband)
      call memget ('alloctv' ,palbgrd,kpt*mband)
      call memget ('alloctv' ,palbgri,kpt*mband)
      call memget ('alloctv' ,pfabd  ,kpt*mband)
      call memget ('alloctv' ,pfabi  ,kpt*mband)
      call memget ('alloctv' ,pftdd  ,kpt*mband)
      call memget ('alloctv' ,pftid  ,kpt*mband)
      call memget ('alloctv' ,pftii  ,kpt*mband)
      call memget ('alloctv' ,pfsun  ,kpt      )

      return
      end



 
