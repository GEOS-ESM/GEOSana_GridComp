#include <preproc.h>
      subroutine surphy (nstep   ,npt     ,dtlsm   ,conchk  ,loopi   ,
     &                   pgcm    ,tgcm    ,qgcm    ,firgcm  ,thgcm   ,
     &                   fire    ,solad   ,solai   ,egcm    ,rhogcm  ,
     &                   ugcm    ,vgcm    ,hgcm    ,ivt     ,ist     ,
     &                   lati    ,dzsoi   ,zsoi    ,tsoi    ,btran   ,
     &                   tg      ,tv      ,moz     ,h2osoi  ,fwet    ,
     &                   h2osno  ,eah     ,igs     ,fsa     ,fsr     ,
     &                   ndvi    ,albd    ,albi    ,albgrd  ,albgri  ,
     &                   fabd    ,fabi    ,ftdd    ,ftid    ,ftii    ,
     &                   fsun    ,htop    ,elai    ,esai    ,lakpoi  ,
     &                   fcev    ,fctr    ,fgev    ,rssun   ,rssha   ,
     &                   taux    ,tauy    ,trad    ,tsa     ,tah     ,
     &                   bevap   ,psnsun  ,laisun  ,psnsha  ,laisha  ,
     &                   root    ,watsat  ,watdry  ,watopt  ,csol    ,
     &                   tksol   ,tkdry   ,fsh     ,tam     ,fgr     ,
     &                   fsm     ,cn      ,cm      ,dtv     ,fira    ,
     &                   dtg     ,gcan    ,ch      ,co2gcm  ,o2gcm   ,
     &                   pergro  ,errsoi  ,errseb  ,errsol  ,errlon  ,
     &                   hsno    ,fsno    ,foln    ,kpti    ,
     &                   z0m     ,z0h     ,zpd     ,tskin4  ,trad4   )

#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ------------------------------------
* source file:       surphy.F
* purpose:           surface biophysics: energy fluxes and temperatures
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* --------------------------------------------------------------------------

* ------------------------ notes -------------------------------------------
* the calling sequence is:

* surphy:                surface biophysics driver
*   -> surrad:           surface solar radiation
*   -> surtem:           surface temperatures and fluxes
*        -> stomata:     sunlit leaf stomatal resistance and photosynthesis
*        -> stomata:     shaded leaf stomatal resistance and photosynthesis
*   -> soitem:           soil temperatures
*        -> tridia:      tridiagonal matrix solution
*   -> laktem:           lake temperatures
*        -> tridia:      tridiagonal matrix solution

* the model distinguishes among the following temperatures:
*    o tv   = vegetation "skin" temperature that balances canopy energy budget
*    o tg   = ground "skin" temperature that balances ground energy budget
*    o trad = radiative temperature = (fire/sb)**0.25
*    o tam  = "aerodynamic" temp at height z0m+zpd (apparent sink for momentum)
*    o tah  = "surface" temp at height z0h+zpd (apparent sink for sensible heat)
*    o tsa  = 2 m height air temperature

* the model partitions the latent heat flux into three components:
*    o fcev: evaporation of intercepted water
*    o fctr: transpiration
*    o fgev: soil evaporation or snow sublimation

* the model conserves surface energy fluxes as:
*    o -fsa + fira + fsh + (fcev+fctr+fgev) + fcst + fgr + fsm = 0
*    o fsa + fsr = [solad(1)+solad(2)+solai(1)+solai(2)] 
*    o fira = -firgcm + fire
* currently canopy heat storage fcst = 0

* the energy conservation for soil/lake temperatures is:

*    cv(  1)*([tsoi(  1)] n+1 - [tsoi(  1)] n)*dz(  1)/dtlsm + ... +
*    cv(nsl)*([tsoi(nsl)] n+1 - [tsoi(nsl)] n)*dz(nsl)/dtlsm = fin
*          cv         = heat capacity (j/m**3/kelvin)
*          dz         = thickness, adjusted for snow (m)
*          dtlsm      = time step (s)
*          [tsoi] n   = old temperature (kelvin)
*          [tsoi] n+1 = new temperature (kelvin)
*          fin        = heat flux into the soil - heat flux out of soil (w/m**2)
*                     = fgr                                               [soil]
*                     = 0.4*fsa-fira-fsh-fgev-fsm + phi(1)+ ... +phi(nsl) [lake]
* --------------------------------------------------------------------------

* ------------------------ input variables ---------------------------------
      integer nstep         !time index (0, ....)
      integer npt           !number of points
      integer loopi         !"little" vector index (1 -> numlv)
      integer kpti          !first point in "big" kpt vec for loopi "little" vec

      real dtlsm             !lsm time step (s)

      logical pergro         !true if random perturbation growth test
      logical conchk         !true if want conservation checks

* atmospheric forcing

      real pgcm(npt)         !atmospheric pressure (pa)
      real tgcm(npt)         !atmospheric temperature (kelvin)
      real qgcm(npt)         !atmospheric specific humidity (kg/kg)
      real firgcm(npt)       !downward infrared (longwave) radiation (w/m**2)
      real solad(mband,npt)  !direct beam solar radiation (vis, nir) (w/m**2)
      real solai(mband,npt)  !diffuse solar radiation (vis, nir) (w/m**2)
      real ugcm(npt)         !atmospheric wind in eastward direction (m/s)
      real vgcm(npt)         !atmospheric wind in northward direction (m/s)
      real hgcm(npt)         !atmospheric reference height (m)
      real co2gcm(npt)       !atmospheric co2 concentration (pa)
      real o2gcm(npt)        !atmospheric o2 concentration (pa)
      real thgcm(npt)        !atmospheric potential temperature (kelvin)
      real egcm(npt)         !atmospheric vapor pressure (pa)
      real rhogcm(npt)       !density of air (kg/m**3)

* time invariant surface characteristics

      integer ivt(npt)       !vegetation type
      integer ist(npt)       !"soil" type

      logical lakpoi(npt)    !true: lake point

      real lati(npt)         !latitude,  radians (+ = NH, - = SH)
      real watsat(npt)       !saturated volumetric soil water content (porosity)
      real watdry(npt)       !water content when evapotranspiration stops
      real watopt(npt)       !optimal water content for evapotranspiration 
      real csol(npt)         !specific heat capacity soil solids (j/m**3/kelvin)
      real tksol(npt)        !thermal conductivity soil solids (w/m/kelvin)
      real tkdry(npt)        !thermal conductivity dry soil (w/m/kelvin)

      real dzsoi(msl,npt)    !soil layer thickness (m)
      real zsoi(msl,npt)     !soil layer depth (m)
      real root(msl,npt)     !root fraction

* land surface variables from previous time step: needed because hydrology
* and ecosystem dynamics done after surface fluxes/temperatures

      real h2osno(npt)       !water equivalent snow (mm h2o / m**2)
      real h2osoi(msl,npt)   !volumetric soil water content (0<=h2osoi<=watsat)
      real hsno(npt)         !snow height (m)
      real fsno(npt)         !fraction of ground covered with snow (0 to 1)
      real fwet(npt)         !fraction of canopy that is wet (0 to 1)
      real igs(npt)          !growing season index (0=off, 1=on)
      real foln(npt)         !foliage nitrogen (%)
      real htop(npt)         !canopy top (m)
      real elai(npt)         !leaf area index, one-sided, with burying by snow
      real esai(npt)         !stem area index, one-sided, with burying by snow

* land surface variables from previous time step: needed because albedos 
* are done after surface fluxes/temperatures for next time step

      real albd(mband,npt)   !surface albedo (direct)
      real albi(mband,npt)   !surface albedo (diffuse)
      real albgrd(mband,npt) !ground  albedo (direct)
      real albgri(mband,npt) !ground  albedo (diffuse)
      real fabd(mband,npt)   !flux absorbed by veg (per unit direct flux) 
      real fabi(mband,npt)   !flux absorbed by veg (per unit diffuse flux) 
      real ftdd(mband,npt)   !downward  direct flux below veg(per unit dir flux)
      real ftid(mband,npt)   !downward diffuse flux below veg(per unit dir flux)
      real ftii(mband,npt)   !downward diffuse flux below veg(per unit dif flux)
      real fsun(npt)         !sunlit fraction of canopy
* --------------------------------------------------------------------------

* ------------------------ input/output variables --------------------------
      real tv(npt)           !vegetation temperature (kelvin)
      real tg(npt)           !ground temperature (kelvin)
      real tsoi(msl,npt)     !soil temperature (kelvin)
      real eah(npt)          !canopy air vapor pressure (pa)
      real moz(npt)          !monin-obukhov stability parameter 
* --------------------------------------------------------------------------

* ------------------------ output variables --------------------------------
* co2 fluxes

      real rssun(npt)   !sunlit stomatal resistance (s/m)
      real rssha(npt)   !shaded stomatal resistance (s/m)
      real psnsun(npt)  !sunlit leaf photosynthesis (umol co2 /m**2/ s)
      real psnsha(npt)  !shaded leaf photosynthesis (umol co2 /m**2/ s)
      real laisun(npt)  !sunlit leaf area
      real laisha(npt)  !shaded leaf area

* hydrology

      real bevap(npt)  !soil water evaporation factor (0 to 1)
      real btran(npt)  !soil water transpiration factor (0 to 1)

* surface energy fluxes

      real fira(npt)   !net infrared (longwave) radiation (w/m**2) [+ = to atm]
      real fire(npt)   !emitted ir (w/m**2)
      real fgr(npt)    !soil heat flux, after snow melt (w/m**2) [+ = into soil]
      real fsm(npt)    !snow melt heat flux (w/m**2) 
      real taux(npt)   !zonal wind stress: e-w (kg/m/s**2)
      real tauy(npt)   !meridional wind stress: n-s (kg/m/s**2)
      real fsh(npt)    !total sensible heat flux (w/m**2) [+ = to atm]
      real fcev(npt)   !canopy evaporation heat flux (w/m**2) [+ = to atm]
      real fgev(npt)   !ground evaporation heat flux (w/m**2) [+ = to atm]
      real fctr(npt)   !transpiration heat flux (w/m**2) [+ = to atm]
      real gcan(npt)   !canopy conductance for latent heat flux (m/s)
      real cn(npt)     !momentum drag coefficient (neutral)
      real cm(npt)     !momentum drag coefficient 
      real ch(npt)     !sensible heat exchange coefficient 
      real tam(npt)    !"aerodynamic" temperature at height z0m+zpd (kelvin)
      real tah(npt)    !"surface" temperature at height z0h+zpd (kelvin)
      real tsa(npt)    !2 m surface air temperature (kelvin)
      real trad(npt)   !radiative temperature (kelvin)
      real dtv(npt)    !change in tv, last iteration (kelvin)
      real dtg(npt)    !change in tg, last iteration (kelvin)

* surface solar radiation 

      real ndvi(npt)   !normalized difference vegetation index
      real fsa(npt)    !absorbed solar radiation (w/m**2)
      real fsr(npt)    !reflected solar radiation (w/m**2)

* energy conservation errors

      real errsol(npt) !solar radiation energy conservation error (w/m**2)
      real errlon(npt) !longwave radiation energy conservation error (w/m**2)
      real errseb(npt) !surface energy conservation error (w/m**2)
      real errsoi(npt) !soil/lake temperature energy conservation error (w/m**2)

! JDC modify begin
      real tskin4(npt) ! 4th power of skin temperature
      real trad4(npt)  ! 4th power of radiative temperature
! JDC modify end
* --------------------------------------------------------------------------

* --------------- common block variables -----------------------------------
#include <phycon.h>
#include <radcon.h>
#include <snocon.h>
#include <soicon.h>
#include <vegcon.h>
* --------------------------------------------------------------------------

* --------------- local variables ------------------------------------------

      integer nband     !number of solar radiation waveband classes
      integer nsl       !number of soil layers
      integer k,j       !do loop indices
      real dum          !temporary variable
      real mpe          !prevents overflow for division by zero
      real tinc         !freezing/thawing temp range for apparent heat capacity
      real vlh          !volumetric latent heat of fusion (j/m**3)
      real tku          !unfrozen thermal conductivity (w/m/kelvin)
      real tkf          !frozen thermal conductivity (w/m/kelvin)
      real cvu          !unfrozen volumetric specific heat (j/m**3/kelvin)
      real cvf          !frozen volumetric specific heat (j/m**3/kelvin)
      real ds           !depth of snow (to maximum of 1 m)

      real tksur(npt)   !thermal conductivity of snow/soil (w/m/kelvin)
      real cvsur(npt)   !heat capacity of snow/soil (j/m**3/kelvin)
      real dzsur(npt)   !depth of snow/soil (m)
      real tk(msl,npt)  !thermal conductivity (w/m/kelvin)
      real cv(msl,npt)  !volumetric specific heat (j/m**3/kelvin)

      logical veg(npt)  !true: vegetated for surface temperature calculation

      real binv(npt)    !1/vkc*ln(z0m/z0h)
      real cw(npt)      !water vapor exchange coefficient 
      real cwp(npt)     !empirical canopy wind parameter
      real emv(npt)     !vegetation emissivity
      real emg(npt)     !ground emissivity
      real fcst(npt)    !canopy heat storage (w/m**2)
      real fsha(npt)    !shaded fraction of canopy
      real fv(npt)      !friction velocity (m/s)
      real mol(npt)     !monin-obukhov length (m)
      real parsun(npt)  !average vis solar rad absorbed by sunlit leaf (w/m**2)
      real parsha(npt)  !average vis solar rad absorbed by shaded leaf (w/m**2)
      real ram(npt,2)   !aerodynamic resistance for momentum (s/m) 
      real rah(npt,2)   !aerodynamic resistance for sensible heat (s/m) 
      real raw(npt,2)   !aerodynamic resistance for water vapor (s/m) 
      real rb(npt)      !bulk leaf boundary layer resistance (s/m)
      real rsurf(npt)   !soil surface resistance (s/m)
      real sag(npt)     !solar radiation absorbed by ground (w/m**2)
      real sav(npt)     !solar radiation absorbed by vegetation (w/m**2)
      real vai(npt)     !total leaf area index + stem area index, one sided
      real zpd(npt)     !displacement height (m)
      real z0mg(npt)    !ground roughness length, momentum (m)
      real z0hg(npt)    !ground roughness length, sensible heat (m)
      real z0m(npt)     !roughness length, momentum (m)
      real z0h(npt)     !roughness length, sensible heat (m)

      character*14 conlab(4) !labels for conservation error checks
      real crit(4)           !critical levels for conservation checks
      real conerr(4,npt)     !conservation errors

      data conlab / 'sol rad error:', 'sur flx error:',
     &              'lon rad error:', 'soi tem error:'/
      data crit / 0.10, 0.10, 0.10, 0.10/
      save conlab, crit

* --------------------------------------------------------------------------

      if (npt.eq.0) return

* --------------------------------------------------------------------------
* miscellaneous variables
* --------------------------------------------------------------------------

      tinc = 0.5
      mpe = 1.e-06
      nsl = msl
      nband = mband

      do k = 1, npt

* vegetation and ground emissivity

         emv(k) = 1. - exp(-(elai(k)+esai(k))/avmuir)
         emg(k) = eg(ist(k))*(1.-fsno(k)) + emsno*fsno(k)

* set wetted fraction of canopy to zero for random perturbation
* test because otherwise have rapid solution separation

         if (pergro) fwet(k) = 0.

* sunlit/shaded leaf areas. surface type (vegetated, not vegetated)
* for surface temperatures and fluxes

         fsha(k) = 1.-fsun(k) 
         laisun(k) = elai(k)*fsun(k)
         laisha(k) = elai(k)*fsha(k)
         vai(k) = elai(k)+ esai(k)
         if (vai(k) .gt. 0.) then
            veg(k) = .true.
         else
            veg(k) = .false.
         end if

* define binv for each point (so can set = 0 for lakes)

         if (ist(k).eq.istdlak .or. ist(k).eq.istslak) then
            binv(k) = 0.
         else
            binv(k) = binvvt(ivt(k))
         end if

      end do

* --------------------------------------------------------------------------
* roughness lengths and displacement height. these are based on vegetation
* type but should be based on vegetation structure (height, leaf area).
* canopy wind parameter. this is based on vegetation type but should be 
* based on leaf area
* --------------------------------------------------------------------------

      do k = 1, npt

* ground roughness length

         if (tg(k) .le. tfrz .and. lakpoi(k)) then               !frozen water
            z0mg(k) = fsno(k)*rlsno + (1.-fsno(k))*0.04
         else                                                    !soil or 
            z0mg(k) = fsno(k)*rlsno + (1.-fsno(k))*rlsoi(ist(k)) !unfrozen water
         end if
         z0hg(k) = z0mg(k) / exp(vkc*binv(k))

* overall roughness length and displacement height

         if (veg(k)) then
            z0m(k) = z0mvt(ivt(k))
            zpd(k) = zpdvt(ivt(k))
         else
            z0m(k) = z0mg(k)
            zpd(k) = 0.
         end if
         z0h(k) = z0m(k) / exp(vkc*binv(k))

* canopy wind parameter

         cwp(k) = cwpvt(ivt(k))

      end do

* --------------------------------------------------------------------------
* thermal properties. done here because first layer thermal 
* conductivity needed for soil heat flux in surtem
* --------------------------------------------------------------------------

      do j = 1, nsl
         do k = 1, npt
            vlh = h2osoi(j,k)*hfus*denh2o
            tku = ( tksol(k)**(1.-watsat(k)) * tkwat**h2osoi(j,k) -
     &            tkdry(k) ) * h2osoi(j,k)/watsat(k) + tkdry(k)
            tkf = ( tksol(k)**(1.-watsat(k)) * tkice**h2osoi(j,k) -
     &            tkdry(k) ) * h2osoi(j,k)/watsat(k) + tkdry(k)
            cvu = (1.-watsat(k))*csol(k) + cwat*h2osoi(j,k)
            cvf = (1.-watsat(k))*csol(k) + cice*h2osoi(j,k)
            if (tsoi(j,k) .gt. tfrz+tinc) then
               tk(j,k) = tku
               cv(j,k) = cvu
            end if
            if (tsoi(j,k) .ge. tfrz-tinc .and.
     &          tsoi(j,k) .le. tfrz+tinc) then
               tk(j,k) = tkf + (tku-tkf) *
     &                   (tsoi(j,k)-tfrz+tinc) / (2.*tinc)
               cv(j,k) = (cvf+cvu)/2. + vlh/(2.*tinc)
            end if
            if (tsoi(j,k) .lt. tfrz-tinc) then
               tk(j,k) = tkf
               cv(j,k) = cvf
            end if
         end do
      end do

* only need tk(1) for lakes

      do k = 1, npt
         if (lakpoi(k)) then
            if (tg(k) .gt. tfrz) then
               tk(1,k) = tkwat
            else
               tk(1,k) = tkice
            end if
         end if
      end do

* blend snow into first soil layer

      do k = 1, npt
         ds = min( hsno(k), 1.)
         dzsur(k) = dzsoi(1,k) + ds
         tksur(k) = tksno*tk(1,k)*dzsur(k)/(tksno*dzsoi(1,k)+tk(1,k)*ds)
         cvsur(k) = cvsno*cv(1,k)*dzsur(k)/(cvsno*dzsoi(1,k)+cv(1,k)*ds)
         tk(1,k) = tksur(k)
         cv(1,k) = cvsur(k)
      end do

* --------------------------------------------------------------------------
* beta factors: separate bevap into soil (ist=1) and non-soil. can do
* btran for all types because root=0 if not vegetated 
* --------------------------------------------------------------------------

      do k = 1, npt
         if (ist(k). eq. 1) then
            if (tsoi(1,k) .gt. tfrz) then
               bevap(k) = min( max(h2osoi(1,k)-watdry(k),0.) /
     &                    (watopt(k)-watdry(k)), 1. ) 
            else
               bevap(k) = 0.01
            end if
         else
            bevap(k) = 1.
         end if
         btran(k) = 0.
      end do

      do j = 1, nsl
         do k = 1, npt
            if (tsoi(j,k) .gt. tfrz) then
               dum = min( max(h2osoi(j,k)-watdry(k),0.) /
     &               (watopt(k)-watdry(k)), 1. )  ** vw(ivt(k))
            else
               dum = 0.01
            end if
            btran(k) = btran(k) + root(j,k)*dum
         end do
      end do

* --------------------------------------------------------------------------
* solar radiation 
* --------------------------------------------------------------------------

      call surrad (nband  ,npt    ,laisun ,fsa   ,sav    ,sag    ,
     &             parsun ,mpe    ,solad  ,solai ,albgrd ,albgri ,
     &             fabd   ,fabi   ,ftdd   ,ftid  ,ftii   ,fsr    ,
     &             ndvi   ,albd   ,albi   ,fsun  ,fsha   ,elai   ,
     &             vai    ,parsha ,laisha )

      do k = 1, npt
         errsol(k) = fsa(k) + fsr(k) - (solad(1,k) +
     &               solad(2,k) + solai(1,k) + solai(2,k))
      end do

* --------------------------------------------------------------------------
* surface energy exchange and temperatures
* --------------------------------------------------------------------------

      call surtem (npt    ,fv     ,cwp    ,emv    ,emg    ,
     &             dtlsm  ,rah    ,tsoi   ,sav    ,sag    ,
     &             raw    ,fgev   ,fctr   ,fcev   ,eah    ,
     &             rb     ,tah    ,co2gcm ,parsun ,qgcm   ,
     &             rssun  ,laisun ,fwet   ,tam    ,ivt    ,
     &             egcm   ,pgcm   ,rsurf  ,tgcm   ,fsh    ,
     &             fgr    ,fire   ,vai    ,moz    ,z0hg   ,
     &             firgcm ,fcst   ,veg    ,tksur  ,dzsur  ,
     &             tg     ,fira   ,fsm    ,htop   ,igs    ,
     &             mpe    ,h2osno ,fsno   ,cm     ,btran  ,
     &             taux   ,tauy   ,ugcm   ,vgcm   ,hgcm   ,
     &             o2gcm  ,psnsun ,bevap  ,tsa    ,rhogcm ,
     &             ch     ,cn     ,gcan   ,z0mg   ,ram    ,
     &             mol    ,tv     ,dtv    ,dtg    ,trad   ,
     &             binv   ,zpd    ,z0m    ,z0h    ,cw     ,
     &             parsha ,rssha  ,psnsha ,laisha ,thgcm  ,
     &             pergro ,foln   )

      do k = 1, npt
         errseb(k) = -fsa(k) + fira(k) + fsh(k) + fcev(k) +
     &                fctr(k) + fgev(k)+ fcst(k) + fgr(k) + fsm(k)
         errlon(k) = fire(k) - fira(k) - firgcm(k)
      end do

! JDC modify begin
! --------------------------------------------------------------------------
! diagnostics bulk surface skin temperature for DAO NASA
! --------------------------------------------------------------------------

      call diagtskin (npt     ,veg   ,emv  ,emg  ,tv  ,tg  ,
     &                tskin4  ,trad4 )

! JDC modify end

* --------------------------------------------------------------------------
* soil/lake temperatures
* --------------------------------------------------------------------------

      call soitem (npt ,nsl ,dtlsm,tsoi   ,dzsoi ,dzsur ,
     &             fgr ,tk  ,cv   ,lakpoi ,errsoi)

      call laktem (npt   ,nsl   ,errsoi,dtlsm ,nstep  ,
     &             zsoi  ,dzsoi ,tsoi  ,sag   ,lakpoi ,
     &             fira  ,fsh   ,fgev  ,lati  ,fv     ,
     &             fsm   ,tg    ,z0mg  ,ist   ,istdlak)

* --------------------------------------------------------------------------
* check for conservation of energy 
* --------------------------------------------------------------------------

      if (conchk) then
         do k = 1, npt
            conerr(1,k) = errsol(k) 
            conerr(2,k) = errseb(k) 
            conerr(3,k) = errlon(k) 
            conerr(4,k) = errsoi(k) 
            do j = 1, 4
               if (abs(conerr(j,k)) .gt. crit(j)) 
     &         write(6,1001)conlab(j),nstep,loopi,k,k+kpti-1,conerr(j,k)
            end do
         end do
1001     format (1x,a14,' nstep =',i10,' vector =',i3,
     &   ' point =',i3,' lsm point =',i6,' inbalance =',f8.2,' w/m2')
      end if

      return
      end
 
