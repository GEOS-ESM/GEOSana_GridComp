* ------------------------ code history ---------------------------
* source file:       lsmhisdyn.h
* purpose:           land surface model history and restart file 
*                    common block for variables determined with 
*                    dynamic memory allocation
* date last revised: August 1996 
* author:            M. Vertenstein
* standardized:
* reviewed:
* -----------------------------------------------------------------

* lpt and kpt will be determined at run time
* set up pointers so that can dynamically allocate memory for
* arrays dependent on these lengths

      common /lsmhis_p/ pslfval, pmlfval

      pointer (pslfval, slfval)
      pointer (pmlfval, mlfval)

      real slfval(kpt,    mslflds) !accumulated single-level field value
      real mlfval(kpt,msl,mmlflds) !accumulated multi-level  field value

* ------------------------ end lsmhisdyn.h ------------------------
 
