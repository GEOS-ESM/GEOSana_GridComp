#include <preproc.h>
      block data phyconi

#include <vartyp.h>

* ------------------------ code history ---------------------------
* source file:       phyconi.F
* purpose:           initialize physical constants
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

#include <phycon.h>

      data grav   /9.80616  / !acceleration due to gravity (m/s**2)
      data sb     /5.67e-08 / !stefan-boltzmann constant (w/m**2/kelvin**4)
      data cpair  /1004.64  / !heat capacity dry air at const pres (j/kg/kelvin)
      data rair   /287.04   / !gas constant for dry air (j/kg/kelvin)
      data vkc    /0.40     / !von karman constant
      data hvap   /2.5104e06/ !latent heat of vaporization (j/kg)
      data hsub   /2.8440e06/ !latent heat of sublimation (j/kg)
      data hfus   /0.3336e06/ !latent heat of fusion (j/kg)
      data cwat   /4.188e06 / !specific heat capacity of water (j/m**3/kelvin)
      data cice   /2.094e06 / !specific heat capacity of ice (j/m**3/kelvin)
      data tkwat  /0.6      / !thermal conductivity of water (w/m/kelvin)
      data tkice  /2.2      / !thermal conductivity of ice (w/m/kelvin)
      data denh2o /1000.    / !density of water (kg/m**3)
      data tfrz   /273.16   / !freezing point (kelvin)

      end
 
