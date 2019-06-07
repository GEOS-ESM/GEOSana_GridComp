#include <preproc.h>
      subroutine lsm (nstep      ,npt        ,dtlsm      ,dtsoi      ,
     &    doalb      ,pgcm       ,tgcm       ,qgcm       ,firgcm     ,
     &    lati       ,long       ,solad      ,solai      ,co2gcm     ,
     &    o2gcm      ,qprecc     ,qprecl     ,thgcm      ,rhogcm     ,
     &    egcm       ,igs        ,hsno       ,fsno       ,fwet       ,
     &    ugcm       ,vgcm       ,hgcm       ,ivt        ,ist        ,
     &    isc        ,dzsoi      ,zsoi       ,tsoi       ,foln       ,
     &    tg         ,tv         ,moz        ,h2osoi     ,h2ocan     ,
     &    h2osno     ,soot       ,eah        ,fsh        ,frg        ,
     &    fire       ,albd       ,albi       ,albgrd     ,albgri     ,
     &    fabd       ,fabi       ,ftdd       ,ftid       ,ftii       ,
     &    fsun       ,htop       ,elai       ,esai       ,tlai       ,
     &    fcev       ,fctr       ,fgev       ,qsoil      ,qvege      ,
     &    taux       ,tauy       ,trad       ,tsa        ,stemb      ,
     &    fpsn       ,frm        ,fco2       ,tsai       ,rootb      ,
     &    fmicr      ,root       ,watsat     ,hksat      ,soilc      ,
     &    smpsat     ,bch        ,watdry     ,watopt     ,csol       ,
     &    tksol      ,tkdry      ,qvegt      ,kmo        ,kda        ,
     &    calday     ,hydro      ,pergro     ,conchk     ,qover      ,
     &    qdrai      ,loopi      ,kpti       ,
     &    z0m        ,z0h        ,zpd        ,wet1       ,wet        ,
     &    tskin4     ,trad4      )

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       lsm.F
* purpose:           process a "litle" vector of length npt
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* one-dimensional (vector) driver for [npt] points. calling sequence is:

* lsm:                  1-d driver
* -> surphy:            surface biophysics: energy fluxes and temperatures 
*      -> surrad:       surface solar radiation
*      -> surtem:       surface temperatures and fluxes
*           -> stomata: sunlit leaf stomatal resistance and photosynthesis
*           -> stomata: shaded leaf stomatal resistance and photosynthesis
*      -> soitem:       soil temperatures
*           -> tridia:  tridiagonal matrix solution
*      -> laktem:       lake temperatures
*           -> tridia:  tridiagonal matrix solution
* -> verh2o:            vertical (column) surface and soil hydrology 
*      -> canh2o:       canopy water
*      -> snoh2o:       snow water
*      -> infil:        infiltration for large-scale precipitation region
*      -> infil:        infiltration for convective + large-scale precip region
*      -> soih2o:       soil hydrology
*           -> tridia:  tridiagonal matrix solution
* -> surbgc:            surface biogeochem fluxes: co2 respiration and npp
* -> ecodyn:            ecosystem dynamics: phenology, vegetation, soil carbon 
*      -> phenol:       leaf phenology
* -> suralb:            albedos for next time step 
*      -> lsmzen:       cosine solar zenith angle for next time step
*      -> snoalb:       snow albedos: direct beam
*      -> snoalb:       snow albedos: diffuse
*      -> soialb:       soil/lake albedos
*      -> twostr:       absorbed, reflected, transmitted solar fluxes (vis dir)
*      -> twostr:       absorbed, reflected, transmitted solar fluxes (vis dif)
*      -> twostr:       absorbed, reflected, transmitted solar fluxes (nir dir)
*      -> twostr:       absorbed, reflected, transmitted solar fluxes (nir dif)
* -> histslf:           accumulate single-level field over history time interval
* -> histmlf:           accumulate multi-level field over history time interval
* -----------------------------------------------------------------

* ------------------- input variables -----------------------------
      integer nstep         !time step index
      integer npt           !number of points
      integer loopi         !"little" vector index (1 -> numlv)
      integer kpti          !first point in "big" kpt vec for loopi "little" vec
      integer kmo           !month (1, ..., 12)
      integer kda           !day of month (1, ..., 31)

      logical doalb          !true = surface albedo calculation time step
      logical hydro          !true = using prognostic hydrology
      logical pergro         !true = random perturbation growth test
      logical conchk         !true = energy and water conservation check

      real calday            !calendar day at greenwich (1.00, ..., 365.99)
      real dtlsm             !main lsm time step (s)
      real dtsoi             !soil hydrology time step (s)

* atmospheric forcing
* -----------------------------------------------------------------
      real pgcm(npt)         !atmospheric pressure (pa)
      real tgcm(npt)         !atmospheric temperature (kelvin)
      real thgcm(npt)        !atmospheric potential temperature (kelvin)
      real qgcm(npt)         !atmospheric specific humidity (kg/kg)
      real egcm(npt)         !atmospheric vapor pressure (pa)
      real ugcm(npt)         !atmospheric wind speed in east direction (m/s)
      real vgcm(npt)         !atmospheric wind speed in north direction (m/s)
      real co2gcm(npt)       !atmospheric co2 concentration (pa)
      real o2gcm(npt)        !atmospheric o2 concentration (pa)
      real rhogcm(npt)       !density of air (kg/m**3)
      real hgcm(npt)         !atmospheric reference height (m)
      real firgcm(npt)       !downward infrared (longwave) radiation (w/m**2)
      real solad(mband,npt)  !direct beam solar radiation (w/m**2) (vis, nir)
      real solai(mband,npt)  !diffuse solar radiation (w/m**2) (vis, nir)
      real qprecc(npt)       !convective precipitation (mm h2o /s)
      real qprecl(npt)       !large-scale precipitation (mm h2o /s)

* time invariant surface properties
* -----------------------------------------------------------------
      integer ivt(npt)       !vegetation type
      integer ist(npt)       !"soil texture" type
      integer isc(npt)       !soil color type

      real lati(npt)         !latitude,  radians (+ = NH, - = SH)
      real long(npt)         !longitude, radians (+ = EH, - = WH)
      real watsat(npt)       !saturated volumetric soil water content (porosity)
      real watdry(npt)       !water content when evapotranspiration stops
      real watopt(npt)       !optimal water content for evapotranspiration 
      real hksat(npt)        !hydraulic conductivity at saturation (mm h2o /s)
      real smpsat(npt)       !soil matrix potential at saturation (mm) 
      real bch(npt)          !clapp and hornberger "b"
      real csol(npt)         !specific heat capacity soil solids (j/m**3/kelvin)
      real tksol(npt)        !thermal conductivity soil solids (w/m/kelvin)
      real tkdry(npt)        !thermal conductivity dry soil (w/m/kelvin)

      real dzsoi(msl,npt)    !soil layer thickness (m)
      real zsoi(msl,npt)     !soil layer depth (m)
      real root(msl,npt)     !root fraction
* -----------------------------------------------------------------

* ------------------- input/output variables ----------------------

* main variables for restart
* -----------------------------------------------------------------
      real h2osno(npt)       !water equivalent snow (mm h2o / m**2)
      real h2ocan(npt)       !canopy water (mm h2o / m**2)
      real h2osoi(msl,npt)   !volumetric soil water content (0<=h2osoi<=watsat)
      real tv(npt)           !vegetation temperature (kelvin)
      real tg(npt)           !ground temperature (kelvin)
      real tsoi(msl,npt)     !soil temperature (kelvin)
      real eah(npt)          !canopy air vapor pressure (pa)
      real soot(npt)         !soot content of snow (mass fraction)
      real moz(npt)          !monin-obukhov stability parameter 
      real hsno(npt)         !snow height (m)
      real fsno(npt)         !fraction of ground covered with snow (0 to 1)
      real fwet(npt)         !fraction of canopy that is wet (0 to 1)

* output from ecodyn: vegetation structure for next time step
* -----------------------------------------------------------------
      real htop(npt)         !canopy top (m)
      real tlai(npt)         !one-sided leaf area index without burying by snow
      real tsai(npt)         !one-sided stem area index without burying by snow
      real elai(npt)         !one-sided leaf area index with    burying by snow
      real esai(npt)         !one-sided stem area index with    burying by snow
      real igs(npt)          !growing season index (0=off, 1=on)
      real foln(npt)         !foliage nitrogen (%)
      real stemb(npt)        !stem biomass (kg /m**2)
      real rootb(npt)        !root biomass (kg /m**2)
      real soilc(npt)        !soil carbon (kg c /m**2)

* output from suralb: albedos for next time step
* -----------------------------------------------------------------
      real albd(mband,npt)   !surface albedo (direct)
      real albi(mband,npt)   !surface albedo (diffuse)
      real albgrd(mband,npt) !ground  albedo (direct)
      real albgri(mband,npt) !ground  albedo (diffuse)
      real fabd(mband,npt)   !flux absorbed by veg (per unit direct flux) 
      real fabi(mband,npt)   !flux absorbed by veg (per unit diffuse flux) 
      real ftdd(mband,npt)   !downward direct flux below veg (per unit dir flux)
      real ftid(mband,npt)   !downward diffuse flux below veg(per unit dir flux)
      real ftii(mband,npt)   !downward diffuse flux below veg(per unit dif flux)
      real fsun(npt)         !sunlit fraction of canopy
* -----------------------------------------------------------------

* ------------------- surphy output variables ---------------------
      real taux(npt)         !wind stress: e-w (n/m**2)
      real tauy(npt)         !wind stress: n-s (n/m**2)
      real fira(npt)         !net infrared (longwave) rad (w/m**2) [+ = to atm]
      real fire(npt)         !emitted infrared (longwave) radiation (w/m**2) 
      real fsh(npt)          !sensible heat flux (w/m**2) [+ to atm]
      real fcev(npt)         !veg evaporation heat flux (w/m**2) [+ to atm]
      real fgev(npt)         !ground evaporation heat flux (w/m**2) [+ to atm]
      real fctr(npt)         !veg transpiration heat flux (w/m**2) [+ to atm]
      real fgr(npt)          !soil heat flux (w/m**2) [+ = into soil]
      real fsm(npt)          !snow melt heat flux (w/m**2)
      real trad(npt)         !radiative temperature (kelvin)
      real tsa(npt)          !2 m height surface air temperature (kelvin)
      real tam(npt)          !"aerodynamic" temperature, at z0m+zpd (kelvin)
      real tah(npt)          !"surface" temperature, at z0h+zpd (kelvin)
      real cn(npt)           !momentum drag coefficient (neutral)
      real cm(npt)           !momentum drag coefficient
      real ch(npt)           !sensible heat exchange coefficient
      real dtv(npt)          !change in tv, last iteration (kelvin)
      real dtg(npt)          !change in tg, last iteration (kelvin)

      real rssun(npt)        !sunlit stomatal resistance (s/m)
      real rssha(npt)        !shaded stomatal resistance (s/m)
      real psnsun(npt)       !sunlit leaf photosynthesis (umol co2 /m**2/ s) 
      real psnsha(npt)       !shaded leaf photosynthesis (umol co2 /m**2/ s)
      real laisun(npt)       !sunlit leaf area
      real laisha(npt)       !shaded leaf area
      real gcan(npt)         !canopy conductance for latent heat flux (m/s)

      real fsa(npt)          !absorbed solar radiation (w/m**2)
      real fsr(npt)          !reflected solar radiation (w/m**2)
      real ndvi(npt)         !normalized difference vegetation index
      real bevap(npt)        !soil water evaporation   factor (0 to 1)
      real btran(npt)        !soil water transpiration factor (0 to 1)

      real errsoi(npt)       !soil/lake energy conservation error (w/m**2)
      real errseb(npt)       !surface energy conservation error (w/m**2)
      real errsol(npt)       !solar radiation conservation error (w/m**2)
      real errlon(npt)       !longwave radiation conservation error (w/m**2)
! JDC modify begin
      real z0m(npt)          !roughness length, momentum (m)
      real z0h(npt)          !roughness length, sensible heat (m)
      real zpd(npt)          !displacement height (m)
      real wet1(npt)         !top soil layer wetness (h2osoi(1)/watsat)
      real wet(npt)          !root zone soil wetness (rsw/watsat)
      real tskin4(npt)       ! 4th power of skin temperature
      real trad4(npt)        ! 4th power of radiative temperature
! JDC modify end
* -----------------------------------------------------------------

* ------------------- verh2o output variables ---------------------
      real rsw(npt)          !soil water content for root zone
      real qinfl(npt)        !infiltration (mm h2o /s) 
      real qover(npt)        !surface runoff (mm h2o /s) 
      real qdrai(npt)        !sub-surface runoff (mm h2o /s) 
      real qsoil(npt)        !soil evaporation (mm h2o/s) (+ = to atm)
      real qvege(npt)        !vegetation evaporation (mm h2o/s) (+ = to atm)
      real qvegt(npt)        !vegetation transpiration (mm h2o/s) (+ = to atm)
      real qdrip(npt)        !throughfall (mm h2o /s)
      real qmelt(npt)        !snow melt (mm h2o /s)
      real qintr(npt)        !interception (mm h2o /s)
      real errh2o(npt)       !total water error (mm h2o)
* -----------------------------------------------------------------

* ------------------- surbgc output variables ---------------------
      real fpsn(npt)         !photosynthesis (umol co2 /m**2 /s)
      real frm(npt)          !total maintenance respiration (umol co2 /m**2/s)
      real frmf(npt)         !leaf maintenance respiration (umol co2 /m**2 /s)
      real frms(npt)         !stem maintenance respiration (umol co2 /m**2 /s)
      real frmr(npt)         !root maintenance respiration (umol co2 /m**2 /s)
      real frg(npt)          !growth respiration (umol co2 /m**2 /s)
      real fmicr(npt)        !microbial respiration (umol co2 /m**2 /s)
      real fco2(npt)         !net co2 flux (umol co2 /m**2 /s) [+ = to atm]
      real dmi(npt)          !total dry matter production (ug /m**2)
* -----------------------------------------------------------------

* ------------------- common block variables ----------------------
#include <snocon.h>
#include <soicon.h>
#include <vegcon.h>
* -----------------------------------------------------------------

* ------------------- local variables -----------------------------

      integer k              !do loop index
      logical lakpoi(npt)    !true = lake point
      real dum1(npt)         !derived field for history file
      real dum2(npt)         !derived field for history file

* -----------------------------------------------------------------

* tag lake points

      do k = 1, npt
         if (ist(k).eq.istdlak .or. ist(k).eq.istslak) then
            lakpoi(k) = .true.
         else
            lakpoi(k) = .false.
         end if
      end do

* -----------------------------------------------------------------
* surface biophysics: energy fluxes and temperatures
* -----------------------------------------------------------------

      call surphy (nstep   ,npt     ,dtlsm   ,conchk  ,loopi   ,
     &             pgcm    ,tgcm    ,qgcm    ,firgcm  ,thgcm   ,
     &             fire    ,solad   ,solai   ,egcm    ,rhogcm  ,
     &             ugcm    ,vgcm    ,hgcm    ,ivt     ,ist     ,
     &             lati    ,dzsoi   ,zsoi    ,tsoi    ,btran   ,
     &             tg      ,tv      ,moz     ,h2osoi  ,fwet    ,
     &             h2osno  ,eah     ,igs     ,fsa     ,fsr     ,
     &             ndvi    ,albd    ,albi    ,albgrd  ,albgri  ,
     &             fabd    ,fabi    ,ftdd    ,ftid    ,ftii    ,
     &             fsun    ,htop    ,elai    ,esai    ,lakpoi  ,
     &             fcev    ,fctr    ,fgev    ,rssun   ,rssha   ,
     &             taux    ,tauy    ,trad    ,tsa     ,tah     ,
     &             bevap   ,psnsun  ,laisun  ,psnsha  ,laisha  ,
     &             root    ,watsat  ,watdry  ,watopt  ,csol    ,
     &             tksol   ,tkdry   ,fsh     ,tam     ,fgr     ,
     &             fsm     ,cn      ,cm      ,dtv     ,fira    ,
     &             dtg     ,gcan    ,ch      ,co2gcm  ,o2gcm   ,
     &             pergro  ,errsoi  ,errseb  ,errsol  ,errlon  ,
     &             hsno    ,fsno    ,foln    ,kpti    ,
     &             z0m     ,z0h     ,zpd     ,tskin4  ,trad4   )

* -----------------------------------------------------------------
* vertical (column) surface and soil hydrology
* -----------------------------------------------------------------

      call verh2o (npt     ,ivt     ,elai    ,esai    ,dtsoi   ,
     &             dtlsm   ,qdrip   ,qsoil   ,qvege   ,tgcm    ,
     &             qinfl   ,tg      ,tv      ,qover   ,fsm     ,
     &             h2osno  ,qmelt   ,hsno    ,ist     ,h2ocan  ,
     &             fcev    ,fgev    ,fctr    ,h2osoi  ,qdrai   ,
     &             dzsoi   ,loopi   ,root    ,watsat  ,hksat   ,
     &             smpsat  ,bch     ,rsw     ,qintr   ,qprecl  ,
     &             qprecc  ,errh2o  ,hydro   ,zsoi    ,nstep   ,
     &             conchk  ,fsno    ,lakpoi  ,qvegt   ,kpti    )
!  save root zone soil wetness. gkw
!  save  first layer soil wetness
      do k = 1, npt
        wet(k) = rsw(k) /watsat(k)
        wet1(k) = h2osoi(1,k) /watsat(k)
      end do
* -----------------------------------------------------------------
* surface biogeochemical fluxes: co2 respiration and plant production
* -----------------------------------------------------------------

      call surbgc (npt     ,ivt     ,igs     ,tv      ,tsoi    ,
     &             psnsun  ,laisun  ,psnsha  ,laisha  ,tlai    ,
     &             fpsn    ,frmf    ,frms    ,frmr    ,frm     ,
     &             frg     ,fco2    ,dmi     ,fmicr   ,rsw     ,
     &             btran   ,foln    ,soilc   ,stemb   ,rootb   ,
     &             nstep   ,kpti    ,loopi   )

* -----------------------------------------------------------------
* ecosystem dynamics: phenology, vegetation, soil carbon
* -----------------------------------------------------------------

      call ecodyn (npt     ,kmo     ,kda     ,ivt     ,hsno    ,
     &             lati    ,tv      ,htop    ,elai    ,esai    ,
     &             tlai    ,tsai    ,igs     ,stemb   ,rootb   ,
     &             foln    ,soilc   ,h2ocan  ,fwet    ,doalb   ,
     &             nstep   ,kpti    ,loopi   )

* -----------------------------------------------------------------
* albedos for next time step
* -----------------------------------------------------------------

      if (doalb) then
         call suralb (npt      ,ivt      ,ist      ,isc    ,
     &                lati     ,long     ,calday   ,soot   ,
     &                elai     ,esai     ,tam      ,tg     ,
     &                tv       ,h2osno   ,h2osoi   ,fwet   ,
     &                albgrd   ,albgri   ,albd     ,albi   ,
     &                fabd     ,fabi     ,ftdd     ,ftid   ,
     &                ftii     ,fsun     ,fsno     ,kpti   ,
     &                loopi    ,nstep    )
      end if

* -----------------------------------------------------------------
* accumulate field values over history time interval
* -----------------------------------------------------------------

* single-level fields

      call histslf ('TAM     ', npt,  tam     , kpti, kpti+npt-1)
      call histslf ('TAH     ', npt,  tah     , kpti, kpti+npt-1)
      call histslf ('TV      ', npt,  tv      , kpti, kpti+npt-1)
      call histslf ('TG      ', npt,  tg      , kpti, kpti+npt-1)
      call histslf ('TSA     ', npt,  tsa     , kpti, kpti+npt-1)
      call histslf ('DTV     ', npt,  dtv     , kpti, kpti+npt-1)
      call histslf ('DTG     ', npt,  dtg     , kpti, kpti+npt-1)

      call histslf ('QINFL   ', npt,  qinfl   , kpti, kpti+npt-1)
      call histslf ('QOVER   ', npt,  qover   , kpti, kpti+npt-1)
      call histslf ('QDRAI   ', npt,  qdrai   , kpti, kpti+npt-1)
      call histslf ('QSOIL   ', npt,  qsoil   , kpti, kpti+npt-1)
      call histslf ('QVEGE   ', npt,  qvege   , kpti, kpti+npt-1)
      call histslf ('QVEGT   ', npt,  qvegt   , kpti, kpti+npt-1)
      call histslf ('QPRECC  ', npt,  qprecc  , kpti, kpti+npt-1)
      call histslf ('QPRECL  ', npt,  qprecl  , kpti, kpti+npt-1)
      call histslf ('QDRIP   ', npt,  qdrip   , kpti, kpti+npt-1)
      call histslf ('QINTR   ', npt,  qintr   , kpti, kpti+npt-1)
      call histslf ('QMELT   ', npt,  qmelt   , kpti, kpti+npt-1)
      call histslf ('BEVAP   ', npt,  bevap   , kpti, kpti+npt-1)
      call histslf ('BTRAN   ', npt,  btran   , kpti, kpti+npt-1)
      call histslf ('H2OSNO  ', npt,  h2osno  , kpti, kpti+npt-1)
      call histslf ('H2OCAN  ', npt,  h2ocan  , kpti, kpti+npt-1)

      call histslf ('FMICR   ', npt,  fmicr   , kpti, kpti+npt-1)
      call histslf ('FPSN    ', npt,  fpsn    , kpti, kpti+npt-1)
      call histslf ('FRMF    ', npt,  frmf    , kpti, kpti+npt-1)
      call histslf ('FRMS    ', npt,  frms    , kpti, kpti+npt-1)
      call histslf ('FRMR    ', npt,  frmr    , kpti, kpti+npt-1)
      call histslf ('FRG     ', npt,  frg     , kpti, kpti+npt-1)
      call histslf ('FCO2    ', npt,  fco2    , kpti, kpti+npt-1)
      call histslf ('DMI     ', npt,  dmi     , kpti, kpti+npt-1)

      call histslf ('FSA     ', npt,  fsa     , kpti, kpti+npt-1)
      call histslf ('FSR     ', npt,  fsr     , kpti, kpti+npt-1)
      call histslf ('NDVI    ', npt,  ndvi    , kpti, kpti+npt-1)
      call histslf ('FCTR    ', npt,  fctr    , kpti, kpti+npt-1)
      call histslf ('FCEV    ', npt,  fcev    , kpti, kpti+npt-1)
      call histslf ('FGEV    ', npt,  fgev    , kpti, kpti+npt-1)
      call histslf ('FSH     ', npt,  fsh     , kpti, kpti+npt-1)
      call histslf ('FGR     ', npt,  fgr     , kpti, kpti+npt-1)
      call histslf ('FSM     ', npt,  fsm     , kpti, kpti+npt-1)
      call histslf ('FIRA    ', npt,  fira    , kpti, kpti+npt-1)
      call histslf ('FIRE    ', npt,  fire    , kpti, kpti+npt-1)
      call histslf ('TAUX    ', npt,  taux    , kpti, kpti+npt-1)
      call histslf ('TAUY    ', npt,  tauy    , kpti, kpti+npt-1)
      call histslf ('CN      ', npt,  cn      , kpti, kpti+npt-1)
      call histslf ('CM      ', npt,  cm      , kpti, kpti+npt-1)
      call histslf ('ELAI    ', npt,  elai    , kpti, kpti+npt-1)
      call histslf ('ESAI    ', npt,  esai    , kpti, kpti+npt-1)
      call histslf ('RSSUN   ', npt,  rssun   , kpti, kpti+npt-1)
      call histslf ('RSSHA   ', npt,  rssha   , kpti, kpti+npt-1)

      call histslf ('ERRSOI  ', npt,  errsoi  , kpti, kpti+npt-1)
      call histslf ('ERRSEB  ', npt,  errseb  , kpti, kpti+npt-1)
      call histslf ('ERRSOL  ', npt,  errsol  , kpti, kpti+npt-1)
      call histslf ('ERRH2O  ', npt,  errh2o  , kpti, kpti+npt-1)

      do k = 1, npt
         dum1(k) = sqrt(ugcm(k)**2+vgcm(k)**2)
         dum2(k) = solad(1,k)+solai(1,k)+solad(2,k)+solai(2,k)
      end do

      call histslf ('TGCM    ', npt,  tgcm    , kpti, kpti+npt-1)
      call histslf ('PGCM    ', npt,  pgcm    , kpti, kpti+npt-1)
      call histslf ('QGCM    ', npt,  qgcm    , kpti, kpti+npt-1)
      call histslf ('ZGCM    ', npt,  hgcm    , kpti, kpti+npt-1)
      call histslf ('VMAGGCM ', npt,  dum1    , kpti, kpti+npt-1)
      call histslf ('L_DOWN  ', npt,  firgcm  , kpti, kpti+npt-1)
      call histslf ('S_DOWN  ', npt,  dum2    , kpti, kpti+npt-1)

* multi-level fields

      call histmlf ('TSOI    ', npt,  tsoi    , kpti, kpti+npt-1)
      call histmlf ('H2OSOI  ', npt,  h2osoi  , kpti, kpti+npt-1)

      return
      end
 
