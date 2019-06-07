#include <preproc.h>
      block data radconi

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       radconi.F
* purpose:           initialize miscellaneous radiation constants
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

#include <radcon.h>
      integer i             ! loop index

* saturated soil albedos: 1=vis, 2=nir
      data(albsat(i,1),i=1,8)/0.12,0.11,0.10,0.09,0.08,0.07,0.06,0.05/
      data(albsat(i,2),i=1,8)/0.24,0.22,0.20,0.18,0.16,0.14,0.12,0.10/

* dry soil albedos: 1=vis, 2=nir 
      data(albdry(i,1),i=1,8)/0.24,0.22,0.20,0.18,0.16,0.14,0.12,0.10/
      data(albdry(i,2),i=1,8)/0.48,0.44,0.40,0.36,0.32,0.28,0.24,0.20/

* increase desert and semi-desert albedos for north africa and arabian 
* peninsula: set to lightest color and increase vis by 0.0333 and nir 
* by 0.0667 so broad-band is increased by 0.05. see also subroutine soialb
      data albsat(9,1),albdry(9,1) /0.1533,0.2733/
      data albsat(9,2),albdry(9,2) /0.3067,0.5467/

* albedo land ice: 1=vis, 2=nir
      data (albice(i),i=1,mband) /0.80, 0.55/

* albedo frozen lakes: 1=vis, 2=nir 
      data (alblak(i),i=1,mband) /0.60, 0.40/

* omega,betad,betai for snow 
      data (omegas(i),i=1,mband) /0.8, 0.4/
      data betads, betais /0.5, 0.5/

* average inverse optical depth per unit leaf area for ir 
      data avmuir /1./

      end
 
