* ------------------------ code history ---------------------------
* source file:       vegcon.h
* purpose:           vegetation type constants 
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb 1996
* reviewed:          G. Bonan, Feb 1996
* -----------------------------------------------------------------

      common /vegcon_i/ nic, noveg

      integer nic       !value for irrigated crop 
      integer noveg     !value for not vegetated 

      common /vegcon_r/ vw, rdp, ch2op, dleaf
 
      real vw(mvt)      !btran exponent: [(h2osoi-watdry)/(watopt-watdry)]**vw 
      real rdp(mvt)     !defines root fraction decrease with depth
      real ch2op(mvt)   !maximum intercepted h2o per unit lai+sai (mm) 
      real dleaf(mvt)   !characteristic leaf dimension (m) 

      common /vegcon_r/ c3psn  , kc25   , akc   , ko25  , ako  ,
     &                  vcmx25 , avcmx  , bp    , mp    , qe25 ,
     &                  aqe    , rmf25  , rms25 , rmr25 , arm  ,
     &                  dmcf   , folnmx , tmin

      real c3psn(mvt)   !photosynthetic pathway: 0. = c4, 1. = c3
      real kc25(mvt)    !co2 michaelis-menten constant at 25c (pa) 
      real akc(mvt)     !q10 for kc25 
      real ko25(mvt)    !o2 michaelis-menten constant at 25c (pa) 
      real ako(mvt)     !q10 for ko25 
      real vcmx25(mvt)  !maximum rate of carboxylation at 25c (umol co2/m**2/s)
      real avcmx(mvt)   !q10 for vcmx25 
      real bp(mvt)      !minimum leaf conductance (umol/m**2/s) 
      real mp(mvt)      !slope of conductance-to-photosynthesis relationship 
      real qe25(mvt)    !quantum efficiency at 25c (umol co2 / umol photon)
      real aqe(mvt)     !q10 for qe25 
      real rmf25(mvt)   !leaf maintenance respiration at 25c (umol co2/m**2/s)
      real rms25(mvt)   !stem maintenance respiration at 25c (umol co2/kg bio/s)
      real rmr25(mvt)   !root maintenance respiration at 25c (umol co2/kg bio/s)
      real arm(mvt)     !q10 for maintenance respiration
      real dmcf(mvt)    !co2-to-biomass conversion factor (ug biomass/umol co2) 
      real folnmx(mvt)  !foliage nitrogen concentration when f(n)=1 (%) 
      real tmin(mvt)    !minimum temperature for photosynthesis (kelvin) 

      common /vegcon_r/ xl, rhol, rhos, taul, taus

      real xl(mvt)          !leaf/stem orientation index
      real rhol(mvt,mband)  !leaf reflectance: 1=vis, 2=nir 
      real rhos(mvt,mband)  !stem reflectance: 1=vis, 2=nir 
      real taul(mvt,mband)  !leaf transmittance: 1=vis, 2=nir 
      real taus(mvt,mband)  !stem transmittance: 1=vis, 2=nir 

      common /vegcon_r/ binvvt, z0mvt, zpdvt  , stembvt, rootbvt,
     &                  folnvt, mrp  , soilcvt, cwpvt

      real binvvt(mvt)  !1/vkc*ln(z0m/z0h)
      real z0mvt(mvt)   !momentum roughness length (m) 
      real zpdvt(mvt)   !displacement height (m) 
      real stembvt(mvt) !stem biomass (kg /m**2)
      real rootbvt(mvt) !root biomass (kg /m**2)
      real folnvt(mvt)  !foliage nitrogen concentration (%)
      real mrp(mvt)     !microbial respiration parameter (umol co2 /kg c/ s)
      real soilcvt(mvt) !soil carbon (kg c /m**2)
      real cwpvt(mvt)   !empirical canopy wind parameter

      common /vegcon_r/ tai, gai, hvt, hvb

      real tai(mvt,12)  !monthly leaf area index + stem area index, one-sided
      real gai(mvt,12)  !monthly leaf area index, one-sided
      real hvt(mvt)     !top of canopy (m)
      real hvb(mvt)     !bottom of canopy (m)

* ------------------------ end vegcon.h ---------------------------
 
