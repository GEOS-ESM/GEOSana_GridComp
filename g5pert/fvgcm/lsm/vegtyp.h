* ------------------------ code history --------------------------------
* source file:       vegtyp.h
* purpose:           subgrid plant type and fractional area for surface types
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb 1996
* reviewed:          G. Bonan, Feb 1996
* ----------------------------------------------------------------------

      common /vegtyp_i/ plant
      common /vegtyp_r/ cover 

      integer plant(29,3)       !subgrid plant types for 29 surface types
      real cover(29,3)          !subgrid weights for 29 surface types
 
* ------------------------ end vegtyp.h --------------------------------
 
