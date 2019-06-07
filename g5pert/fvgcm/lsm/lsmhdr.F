#include <preproc.h>
      subroutine lsmhdr 

* ------------------------ code history ------------------------------
* source file:       lsmhdr.F
* purpose:           write lsm header to standard output 
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:
* reviewed:
* --------------------------------------------------------------------

 
      write (6,*) 'LSM land surface model, version 1.0'
      write (6,*) 'This model is described in full by: Bonan, G.B. ',
     & '1996. A land surface model'
      write (6,*) '(LSM version 1.0) for ecological, hydrological, ',
     & 'and atmospheric studies:'
      write (6,*) 'technical description and users guide. NCAR ',
     & 'Technical Note NCAR/TN-417+STR.'
      write (6,*) 'National Center for Atmospheric Research, ',
     & 'Boulder, Colorado.'
      write (6,*)

      return
      end
 
