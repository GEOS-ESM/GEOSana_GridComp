* ------------------------ code history ---------------------------
* source file:       lsmtvdyn.h
* purpose:           lsmtvdyn common block for time-varying (restart) 
*                    variables with dynamic memory allocation 
* date last revised: August 1996 (M. Vertenstein)
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb 1996
* reviewed:          G. Bonan, Feb 1996
* -----------------------------------------------------------------

* lpt and kpt will be determined at run time
* set up pointers so that can dynamically allocate memory for
* arrays dependent on these lengths

      pointer ( ph2osno, h2osno )
      pointer ( ph2ocan, h2ocan )
      pointer ( ph2osoi, h2osoi )
      pointer ( ptv    , tv     )
      pointer ( ptg    , tg     )
      pointer ( ptsoi  , tsoi   )
      pointer ( pmoz   , moz    )
      pointer ( peah   , eah    )
      pointer ( psoot  , soot   )
      pointer ( phsno  , hsno   )
      pointer ( pfsno  , fsno   )
      pointer ( pfwet  , fwet   )
      pointer ( phtop  , htop   )
      pointer ( ptlai  , tlai   )
      pointer ( ptsai  , tsai   )
      pointer ( pelai  , elai   )
      pointer ( pesai  , esai   )
      pointer ( pfoln  , foln   )
      pointer ( pstemb , stemb  )
      pointer ( prootb , rootb  )
      pointer ( psoilc , soilc  )
      pointer ( pigs   , igs    )
      pointer ( palbd  , albd   )
      pointer ( palbi  , albi   )
      pointer ( palbgrd, albgrd )
      pointer ( palbgri, albgri )
      pointer ( pfabd  , fabd   )
      pointer ( pfabi  , fabi   )
      pointer ( pftdd  , ftdd   )
      pointer ( pftid  , ftid   )
      pointer ( pftii  , ftii   )
      pointer ( pfsun  , fsun   )

* main land surface variables needed for restart

      real h2osno(kpt)      !snow water (mm h2o / m**2)
      real h2ocan(kpt)      !canopy water (mm h2o / m**2)
      real h2osoi(msl,kpt)  !volumetric soil water content (0<=h2osoi<=watsat)
      real tv(kpt)          !vegetation temperature (kelvin)
      real tg(kpt)          !ground temperature (kelvin)
      real tsoi(msl,kpt)    !soil temperature (kelvin)
      real moz(kpt)         !monon-obukhov stability parameter
      real eah(kpt)         !canopy air vapor pressure (pa)
      real soot(kpt)        !soot content of snow
      real hsno(kpt)        !snow height (m)
      real fsno(kpt)        !fraction of ground covered with snow (0 to 1)
      real fwet(kpt)        !fraction of canopy that is wet (0 to 1)

      common /lsmtv_p/ ph2osno ,ph2ocan ,ph2osoi ,ptv   ,ptg   ,ptsoi ,
     &                 pmoz    ,peah    ,psoot   ,phsno ,pfsno ,pfwet

* vegetation for next time step

      real htop(kpt)        !vegetation height, top (m)
      real tlai(kpt)        !total leaf area index, one-sided
      real tsai(kpt)        !total stem area index, one-sided
      real elai(kpt)        !exposed leaf area index, one-sided
      real esai(kpt)        !exposed stem area index, one-sided
      real foln(kpt)        !foliage nitrogen (%)
      real stemb(kpt)       !stem biomass (kg /m**2)
      real rootb(kpt)       !root biomass (kg /m**2)
      real soilc(kpt)       !soil carbon (kg c /m**2)
      real igs(kpt)         !growing season index (0=off, 1=on)

      common /lsmtv_p/ phtop ,ptlai  ,ptsai  ,pelai  ,pesai , 
     &                 pfoln ,pstemb ,prootb ,psoilc ,pigs    

* albedo calculation for next time step 
 
      real albd(mband,kpt)     !surface albedo (direct)
      real albi(mband,kpt)     !surface albedo (diffuse)
      real albgrd(mband,kpt)   !ground albedo (direct)
      real albgri(mband,kpt)   !ground albedo (diffuse)
      real fabd(mband,kpt)     !flux absorbed by veg (per unit direct flux) 
      real fabi(mband,kpt)     !flux absorbed by veg (per unit diffuse flux) 
      real ftdd(mband,kpt)     !down direct flux below veg (per unit dir flux) 
      real ftid(mband,kpt)     !down diffuse flux below veg (per unit dir flux) 
      real ftii(mband,kpt)     !down diffuse flux below veg (per unit dif flux) 
      real fsun(kpt)           !sunlit fraction of canopy

      common /lsmtv_p/ palbd ,palbi ,palbgrd ,palbgri ,pfabd ,
     &                 pfabi ,pftdd ,pftid   ,pftii   ,pfsun    

* ------------------------ end lsmtvdyn.h -------------------------
 
