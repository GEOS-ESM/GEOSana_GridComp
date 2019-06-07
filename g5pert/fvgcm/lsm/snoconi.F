#include <preproc.h>
      block data snoconi

#include <vartyp.h>

* ------------------------ code history ---------------------------
* source file:       snoconi.F
* purpose:           initialize snow constants 
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

#include <snocon.h>

      data rlsno /0.04    / !roughness length (m)
      data emsno /0.97    / !emissivity
      data bdsno /250.    / !bulk density (kg/m**3)
      data tksno /0.34    / !thermal conductivity (w/m/kelvin)
      data cvsno /0.525e06/ !volumeteric heat capacity (j/m**3/kelvin)
      data hsnoc /0.05    / !height when snow fully covers ground (m)

      end
 
