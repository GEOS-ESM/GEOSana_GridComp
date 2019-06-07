* ------------------------ code history ---------------------------
* source file:       snocon.h
* purpose:           snow constants
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb 1996
* reviewed:          G. Bonan, Feb 1996
* -----------------------------------------------------------------

      common /snocon_r/ rlsno,emsno,bdsno,tksno,cvsno,hsnoc

      real rlsno    !roughness length (m)
      real emsno    !emissivity
      real bdsno    !bulk density snow (kg/m**3)
      real tksno    !thermal conductivity snow (w/m/kelvin)
      real cvsno    !volumeteric heat capacity of snow (j/m**3/kelvin)
      real hsnoc    !height of snow when ground fully covered by snow (m)

* ------------------------ end snocon.h ---------------------------
 
