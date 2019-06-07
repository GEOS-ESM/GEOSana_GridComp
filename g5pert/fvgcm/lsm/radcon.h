* ------------------------ code history ---------------------------
* source file:       radcon.h
* purpose:           miscellaneous radiation constants
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb 1996
* reviewed:          G. Bonan, Feb 1996
* -----------------------------------------------------------------

      common /radcon_r/ albsat, albdry, albice, alblak,
     &                  omegas, betads, betais, avmuir

      real albsat(msc,mband)   !saturated soil albedos: 1=vis, 2=nir
      real albdry(msc,mband)   !dry soil albedos: 1=vis, 2=nir 
      real albice(mband)       !albedo land ice: 1=vis, 2=nir
      real alblak(mband)       !albedo frozen lakes: 1=vis, 2=nir 
      real omegas(mband)       !two-stream parameter omega for snow
      real betads              !two-stream parameter betad for snow
      real betais              !two-stream parameter betad for snow
      real avmuir              !ir inverse optical depth per unit leaf area 

* ------------------------ end radcon.h ---------------------------
 
