* ------------------------ code history ---------------------------
* source file:       lakcon.h
* purpose:           constants for lake temperature model
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb 1996
* reviewed:          G. Bonan, Feb 1996
* -----------------------------------------------------------------

      common /lakcon_r/ beta(mst),za(mst),eta(mst),p0

      real beta  !fraction solar rad absorbed at surface: depends on lake type
      real za    !base of surface absorption layer (m): depends on lake type
      real eta   !light extinction coefficient (/m): depends on lake type
      real p0    !neutral value of turbulent prandtl number

* ------------------------ end lakcon.h ---------------------------
 
