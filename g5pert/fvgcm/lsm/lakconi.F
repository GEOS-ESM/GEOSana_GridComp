#include <preproc.h>
      block data lakconi 

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       lakconi.F
* purpose:           initialize constants for lake temperature model
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

#include <lakcon.h>

* warning: array index of "3" corresponds to istdlak in soiconi
* warning: array index of "4" corresponds to istslak in soiconi

      data beta(3) /0.40/ !fraction solar rad absorbed at surface: deep lake
      data beta(4) /0.40/ !fraction solar rad absorbed at surface: shallow lake
      data za(3)   /0.60/ !base of surface absorption layer (m): deep lake
      data za(4)   /0.50/ !base of surface absorption layer (m): shallow lake
      data eta(3)  /0.10/ !light extinction coefficient (/m): deep lake
      data eta(4)  /0.50/ !light extinction coefficient (/m): shallow lake
      data p0      /1.  / !neutral value of turbulent prandtl number

      end
 
