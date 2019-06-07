#include <preproc.h>
      subroutine alloctc

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history --------------------------------
* source file:       alloctc.F
* purpose:           dynamically allocate memory for lsmtcdyn.h
* date last revised: August 1996
* author:            M. Vertenstein
* ----------------------------------------------------------------------

* ------------------------ common block variables ----------------------
#include <dynmem.h>
#include <lsmtcdyn.h>
* ----------------------------------------------------------------------

* dynamically allocate memory for lpt and kpt dependent variables 
* in lsmtcdyn.h

      call memget ('alloctc' ,pixy   ,lpt     )           
      call memget ('alloctc' ,pjxy   ,lpt     )          
      call memget ('alloctc' ,pkvec  ,lpt*msub)
      call memget ('alloctc' ,pwsg2g ,lpt*msub)
      call memget ('alloctc' ,pklnd  ,kpt     )          
      call memget ('alloctc' ,pivt   ,kpt     ) 
      call memget ('alloctc' ,pist   ,kpt     ) 
      call memget ('alloctc' ,pisc   ,kpt     ) 
      call memget ('alloctc' ,pwatsat,kpt     ) 
      call memget ('alloctc' ,phksat ,kpt     )  
      call memget ('alloctc' ,psmpsat,kpt     ) 
      call memget ('alloctc' ,pbch   ,kpt     )    
      call memget ('alloctc' ,pwatdry,kpt     ) 
      call memget ('alloctc' ,pwatopt,kpt     ) 
      call memget ('alloctc' ,pcsol  ,kpt     )   
      call memget ('alloctc' ,ptksol ,kpt     )  
      call memget ('alloctc' ,ptkdry ,kpt     )  
      call memget ('alloctc' ,pdzsoi ,kpt*msl )
      call memget ('alloctc' ,pzsoi  ,kpt*msl )
      call memget ('alloctc' ,proot  ,kpt*msl )
      call memget ('alloctc' ,psand  ,kpt     )   
      call memget ('alloctc' ,pclay  ,kpt     )   
      call memget ('alloctc' ,plati  ,kpt     )   
      call memget ('alloctc' ,plong  ,kpt     )
      call memget ('alloctc' ,pbegkpt,numlv   )
      call memget ('alloctc' ,pnumkpt,numlv   )

      return
      end



 
