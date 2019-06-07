#include <preproc.h>
      block data soiconi

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       soiconi.F
* purpose:           initialize "soil type" constants
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

#include <soicon.h>

* "soil" types 
*   1     soil
*   2     land ice (glacier)
*   3     deep lake
*   4     shallow lake
*   5     wetland: swamp, marsh, etc

* emissivity ground surface
      data eg /0.96, 0.97, 0.97, 0.97, 0.96/

* roughness length (m)
      data rlsoi /0.05, 0.05, 0.001, 0.001, 0.05/

      data istice  /2/  !land ice     "soil" type
      data istdlak /3/  !deep lake    "soil" type
      data istslak /4/  !shallow lake "soil" type
      data istwet  /5/  !wetland      "soil" type

      end
 
