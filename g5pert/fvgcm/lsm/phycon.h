* ------------------------ code history ---------------------------
* source file:       phycon.h
* purpose:           physical constants 
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb 1996
* reviewed:          G. Bonan, Feb 1996
* -----------------------------------------------------------------

      common /phycon_r/ grav , sb   , cpair , rair, vkc ,
     &                  hvap , hsub , hfus  , cwat, cice,
     &                  tkwat, tkice, denh2o, tfrz

      real grav    !acceleration due to gravity (m/s**2)
      real sb      !stefan-boltzmann constant (w/m**2/kelvin**4)
      real cpair   !specific heat capacity dry air at const press (j/kg/kelvin)
      real rair    !gas constant for dry air (j/kg/kelvin)
      real vkc     !von karman constant
      real hvap    !latent heat of vaporization (j/kg)
      real hsub    !latent heat of sublimation (j/kg)
      real hfus    !latent heat of fusion (j/kg)
      real cwat    !specific heat capacity of water (j/m**3/kelvin)
      real cice    !specific heat capacity of ice (j/m**3/kelvin)
      real tkwat   !thermal conductivity of water (w/m/kelvin)
      real tkice   !thermal conductivity of ice (w/m/kelvin)
      real denh2o  !density of water (kg/m**3)
      real tfrz    !freezing point (kelvin)

* ------------------------ end phycon.h ---------------------------
 
