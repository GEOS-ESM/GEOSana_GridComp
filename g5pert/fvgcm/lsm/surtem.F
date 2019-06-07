#include <preproc.h>
      subroutine surtem (npt    ,fv     ,cwp    ,emv    ,emg    ,
     &                   dtim   ,rah    ,tsoi   ,sav    ,sag    ,
     &                   raw    ,fgev   ,fctr   ,fcev   ,eah    ,
     &                   rb     ,tah    ,co2gcm ,parsun ,qgcm   ,
     &                   rssun  ,laisun ,fwet   ,tam    ,ivt    ,
     &                   egcm   ,pgcm   ,rsurf  ,tgcm   ,fsh    ,
     &                   fgr    ,fire   ,vai    ,moz    ,z0hg   ,
     &                   firgcm ,fcst   ,veg    ,tksur  ,dzsur  ,
     &                   tg     ,fira   ,fsm    ,htop   ,igs    ,
     &                   mpe    ,h2osno ,fsno   ,cm     ,btran  ,
     &                   taux   ,tauy   ,ugcm   ,vgcm   ,hgcm   ,
     &                   o2gcm  ,psnsun ,bevap  ,tsa    ,rhogcm ,
     &                   ch     ,cn     ,gcan   ,z0mg   ,ram    ,
     &                   mol    ,tv     ,dtv    ,dtg    ,trad   ,
     &                   binv   ,zpd    ,z0m    ,z0h    ,cw     ,
     &                   parsha ,rssha  ,psnsha ,laisha ,thgcm  ,
     &                   pergro ,foln   )

#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       surtem.F
* purpose:           surface temperatures and fluxes
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* use newton-raphson iteration to solve for vegetation (tv) and 
* ground (tg) temperatures that balance the surface energy budgets

* vegetated:
* -sa(1) + ir(1)[tv] + sh(1)[tv] + ev(1)[tv] + tr(1)[tv] + st(1)[tv] = 0 
* -sa(2) + ir(2)[tg] + sh(2)[tg] + ev(2)[tg] + tr(2)[tg] + st(2)[tg] = 0 

* not vegetated:
* -sa(2) + ir(2)[tg] + sh(2)[tg] + ev(2)[tg] + tr(2)[tg] + st(2)[tg] = 0 

* if vegetated: first solve for tv given tg. then solve for tg.

* total fluxes: 
* -fsa + fira + fsh + (fcev + fctr + fgev) + fcst + fgr + fsm = 0 
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
* input
      integer npt       !number of points

      real rhogcm(npt)  !density air (kg/m**3)
      real egcm(npt)    !vapor pressure air at agcm reference height (pa)
      real pgcm(npt)    !pressure at agcm reference height (pa)
      real qgcm(npt)    !specific humidity at agcm reference height (kg/kg)
      real tgcm(npt)    !temperature at agcm reference height (kelvin)
      real thgcm(npt)   !potential temperature at agcm reference height (kelvin)
      real firgcm(npt)  !atmospheric infrared (longwave) radiation (w/m**2)
      real ugcm(npt)    !wind speed in eastward direction (m/s)
      real vgcm(npt)    !wind speed in northward direction (m/s)
      real hgcm(npt)    !reference height for agcm (m)
      real co2gcm(npt)  !atmospheric co2 concentration (pa)
      real o2gcm(npt)   !atmospheric o2 concentration (pa)
      real h2osno(npt)  !water equivalent snow (kg h2o /m**2)
      real vai(npt)     !total leaf area index + stem area index, one sided
      real fsno(npt)    !fraction of ground covered with snow
      real tksur(npt)   !thermal conductivity of snow/soil (w/m/kelvin)
      real dzsur(npt)   !depth of snow/soil (m)
      real igs(npt)     !growing season index (0=off, 1=on)
      real fwet(npt)    !wetted fraction of canopy
      real htop(npt)    !top of canopy layer (m)
      real bevap(npt)   !soil water evaporation factor (0 to 1)
      real btran(npt)   !soil water transpiration factor (0 to 1)
      real sav(npt)     !solar radiation absorbed by vegetation (w/m**2)
      real sag(npt)     !solar radiation absorbed by ground (w/m**2)
      real laisun(npt)  !sunlit leaf area index, one-sided (m**2/m**2)
      real laisha(npt)  !shaded leaf area index, one-sided (m**2/m**2)
      real parsun(npt)  !par absorbed per unit sunlit lai (w/m**2)
      real parsha(npt)  !par absorbed per unit shaded lai (w/m**2)
      real cwp(npt)     !canopy wind parameter
      real binv(npt)    !1/vkc*ln(z0m/z0h)
      real zpd(npt)     !zero plane displacement (m)
      real z0m(npt)     !roughness length, momentum (m)
      real z0h(npt)     !roughness length, sensible heat (m)
      real z0hg(npt)    !roughness length, sensible heat, ground (m)
      real z0mg(npt)    !roughness length, momentum, ground (m)
      real emv(npt)     !vegetation emissivity
      real emg(npt)     !ground emissivity
      real foln(npt)    !foliage nitrogen (%)

      logical veg(npt)  !true if vegetated surface
      logical pergro    !true if random perturbation growth test

      integer ivt(npt)  !vegetation physiology type

      real tsoi(msl,npt)!soil temperature (kelvin)

      real dtim         !time step (s)
      real mpe          !prevents overflow error if division by zero
 
* input/output
      real eah(npt)     !canopy air vapor pressure (pa)
      real tv(npt)      !vegetation temperature (kelvin)
      real tg(npt)      !ground temperature (kelvin)
      real moz(npt)     !monin-obukhov stability parameter 

* output
      real cm(npt)      !momentum drag coefficient 
      real ch(npt)      !sensible heat exchange coefficient 
      real cw(npt)      !water vapor exchange coefficient 
      real cn(npt)      !momentum drag coefficient (neutral)
      real fv(npt)      !friction velocity (m/s)
      real taux(npt)    !wind stress: e-w (n/m**2)
      real tauy(npt)    !wind stress: n-s (n/m**2)
      real rb(npt)      !bulk leaf boundary layer resistance (s/m)
      real mol(npt)     !monin-obukhov length (m) 
      real dtv(npt)     !change in tv, last iteration (kelvin)
      real dtg(npt)     !change in tg, last iteration (kelvin)
      real rssun(npt)   !sunlit leaf stomatal resistance (s/m)
      real rssha(npt)   !shaded leaf stomatal resistance (s/m)
      real psnsun(npt)  !sunlit leaf photosynthesis (umol co2 /m**2/ s) [+]
      real psnsha(npt)  !shaded leaf photosynthesis (umol co2 /m**2/ s) [+]
      real fira(npt)    !total net infrared (longwave) rad (w/m**2) [+ = to atm]
      real fsh(npt)     !total sensible heat flux (w/m**2)       [+ = to atm]
      real fcev(npt)    !canopy evaporation heat flux (w/m**2)   [+ = to atm]
      real fgev(npt)    !ground evaporation heat flux (w/m**2)   [+ = to atm]
      real fctr(npt)    !transpiration heat flux (w/m**2)        [+ = to atm]
      real fcst(npt)    !canopy heat storage (w/m**2)            [+ = to atm]
      real fgr(npt)     !soil heat flux (w/m**2)                 [+ = into soil]
      real fsm(npt)     !snow melt heat flux (w/m**2)             
      real fire(npt)    !emitted ir (w/m**2)
      real trad(npt)    !radiative temperature (kelvin)
      real tam(npt)     !"aerodynamic" temperature at height z0m+zpd (kelvin)
      real tah(npt)     !"surface" temperature at height z0h+zpd (kelvin)
      real tsa(npt)     !2 m height air temperature (kelvin)
      real gcan(npt)    !canopy conductance for latent heat flux (m/s)
      real rsurf(npt)   !ground surface resistance (s/m)
      real ram(npt,2)   !aerodynamic resistance for momentum (s/m) 
      real rah(npt,2)   !aerodynamic resistance for sensible heat (s/m) 
      real raw(npt,2)   !aerodynamic resistance for water vapor (s/m) 
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <phycon.h>
#include <vegcon.h>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------

      real sa(npt,2)    !absorbed solar radiation (w/m**2)      [always +]
      real ir(npt,2)    !net infrared (longwave) radiation (w/m**2) [+ = to atm]
      real sh(npt,2)    !sensible heat flux (w/m**2)            [+ = to atm]
      real ev(npt,2)    !evaporation heat flux (w/m**2)         [+ = to atm]
      real tr(npt,2)    !transpiration heat flux (w/m**2)       [+ = to atm]
      real st(npt,2)    !canopy heat storage/soil heat (w/m**2) [+ = to veg/soi]

      real air(npt),bir(npt)  !coefficients for ir as function of ts**4
      real ash(npt),bsh(npt)  !coefficients for sh as function of ts
      real aev(npt),bev(npt)  !coefficients for ev as function of esat[ts]
      real atr(npt),btr(npt)  !coefficients for tr as function of esat[ts]
      real ast(npt),bst(npt)  !coefficients for st as function of ts
      real ata(npt),bta(npt)  !coefficients for tah as function of ts
      real aea(npt),bea(npt)  !coefficients for eah as function of esat[ts]

      real ts(npt)      !vegetation or ground temperature (kelvin)
      real dts(npt)     !change in temperature (kelvin)
      real svpts(npt)   !saturation vapor pressure at ts (pa)
      real svptg(npt)   !saturation vapor pressure at tg (pa)
      real dsvpts(npt)  !d svp / d t at ts (pa/kelvin)
      real dsvptg(npt)  !d svp / d t at tg (pa/kelvin)
      real fm(npt)      !momentum stability correction, weighted by prior iters
      real fh(npt)      !sen heat stability correction, weighted by prior iters
      real gamma(npt)   !psychrometric constant (pa/kelvin)
      real lathea(npt)  !lathent heat of vaporization/sublimation (j/kg)
      real hcan(npt)    !canopy height (m) [note: hcan >= z0mg]
      real ur(npt)      !wind speed at height hgcm (m/s)
      real tmpra(npt)   !temporary calculation for aerodynamic resistances
      real tmpram(npt)  !temporary calculation for aerodynamic resistances
      real tmprah(npt)  !temporary calculation for aerodynamic resistances
      real tmprah2(npt) !temporary calculation for aerodynamic resistances
      real tmpuc(npt)   !temporary calculation for canopy wind
      real tmprb(npt)   !temporary calculation for rb
      real mozold(npt)  !monin-obukhov stability parameter from prior iteration

      real a            !temporary calculation
      real b            !temporary calculation
      real tmp1         !temporary calculation
      real tmp2         !temporary calculation
      real tmp3         !temporary calculation
      real tmp4         !temporary calculation
      real tmp5         !temporary calculation
      real cah          !sensible heat conductance, surface air to gcm air (m/s)
      real cvh          !sensible heat conductance, veg to surface air (m/s)
      real cgh          !sensible heat conductance, ground to surface air (m/s)
      real caw          !latent heat conductance, surface air to gcm air (m/s)
      real ctw          !transpiration conductance, veg to surface air (m/s)
      real cew          !evaporation conductance, veg to surface air (m/s)
      real cgw          !latent heat conductance, ground to surface air (m/s)
      real cond         !sum of conductances (s/m)
      real uc           !wind speed at top of canopy (m/s)
      real kh           !turbulent transfer coefficient, sensible heat, (m**2/s)
      real fmnew        !stability correction factor, momentum, for current moz
      real fhnew        !stability correction factor, sen heat, for current moz
      real tvir         !temporary virtual temperature (kelvin)
      real h            !temporary sensible heat flux (w/m**2)
      real hvir         !temporary virtual sensible heat flux (w/m**2)
      real e            !temporary water vapor flux (kg/m**2/s)

      integer k,i          !index
      integer iter         !iteration index
      integer niters       !number of iterations for surface temperature
      integer niterg       !number of iterations for ground temperature
      integer nmozsgn(npt) !number of times moz changes sign

      data niters,niterg /10,3/
      save niters,niterg

* use polynomials to calculate saturation vapor pressure and derivative with
* respect to temperature: over water when t > 0 c and over ice when t <= 0 c

      real t                     !temperature 
      real esatw                 !saturation vapor pressure over water (pa)
      real esati                 !saturation vapor pressure over ice (pa)
      real dsatw                 !d(esat)/dt over water (pa/kelvin)
      real dsati                 !d(esat)/dt over ice (pa/kelvin)
      real a0,a1,a2,a3,a4,a5,a6  !coefficients for esat over water
      real b0,b1,b2,b3,b4,b5,b6  !coefficients for esat over ice
      real c0,c1,c2,c3,c4,c5,c6  !coefficients for dsat over water
      real d0,d1,d2,d3,d4,d5,d6  !coefficients for dsat over ice
      real tdc                  !kelvin to degree celsius with limit -50 to +50

      parameter (a0=6.107799961    , a1=4.436518521e-01,
     &           a2=1.428945805e-02, a3=2.650648471e-04,
     &           a4=3.031240396e-06, a5=2.034080948e-08,
     &           a6=6.136820929e-11)
 
      parameter (b0=6.109177956    , b1=5.034698970e-01,
     &           b2=1.886013408e-02, b3=4.176223716e-04,
     &           b4=5.824720280e-06, b5=4.838803174e-08,
     &           b6=1.838826904e-10)
 
      parameter (c0= 4.438099984e-01, c1=2.857002636e-02,
     &           c2= 7.938054040e-04, c3=1.215215065e-05,
     &           c4= 1.036561403e-07, c5=3.532421810e-10,
     &           c6=-7.090244804e-13)
 
      parameter (d0=5.030305237e-01, d1=3.773255020e-02,
     &           d2=1.267995369e-03, d3=2.477563108e-05,
     &           d4=3.005693132e-07, d5=2.158542548e-09,
     &           d6=7.131097725e-12)

      tdc(t,tfrz) = min( 50., max(-50.,(t-tfrz)) )
      esatw(t) = 100.*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
      esati(t) = 100.*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
      dsatw(t) = 100.*(c0+t*(c1+t*(c2+t*(c3+t*(c4+t*(c5+t*c6))))))
      dsati(t) = 100.*(d0+t*(d1+t*(d2+t*(d3+t*(d4+t*(d5+t*d6))))))

* -----------------------------------------------------------------
* initialization variables that do not depend on stability iteration
* -----------------------------------------------------------------

      do k = 1, npt
         dtv(k) = 0.
         dtg(k) = 0.

* set psychrometric constant

         if (tgcm(k) .gt. tfrz) then
            lathea(k) = hvap
         else
            lathea(k) = hsub
         end if
         gamma(k) = cpair*pgcm(k)/(0.622*lathea(k))

* initialize absorbed solar radiation and ts=tv or ts=tg

         sa(k,1) = sav(k)
         sa(k,2) = sag(k)
         if (veg(k)) then
            ts(k) = tv(k)
         else
            ts(k) = tg(k)
         end if

* saturation vapor pressure at ground temperature

         t = tdc(tg(k),tfrz) 
         if (t .gt. 0.) then
            svptg(k) = esatw(t)
         else
            svptg(k) = esati(t)
         end if

* initialize vegetated fluxes (for non-vegetated surfaces)

         ir(k,1) = 0.
         sh(k,1) = 0.
         ev(k,1) = 0.
         tr(k,1) = 0.
         st(k,1) = 0.

      end do

      do k = 1, npt

* canopy height

         hcan(k) = max( htop(k), z0mg(k) )

* wind speed at agcm reference height: ur >= 1 or else large 
* aerodynamic resistances with ur < 1 m/s

         ur(k) = max( sqrt(ugcm(k)**2+vgcm(k)**2), 1. )

* neutral drag coefficient evaluated at z=z0m+zpd

         cn(k) = ( vkc / log((hgcm(k)-zpd(k))/z0m(k)) )**2

* miscellaneous 

         nmozsgn(k) = 0
         mozold(k) = 0.
         tmp1 = exp( -cwp(k)*z0hg(k)/hcan(k) )
         tmp2 = exp( -cwp(k)*(z0h(k)+zpd(k))/hcan(k) )
         tmprah2(k) = hcan(k)*exp(cwp(k)) / cwp(k) * (tmp1-tmp2)
         tmpuc(k) = log((hcan(k)-zpd(k))/z0m(k)) / vkc
         tmprb(k) = cwp(k)*50. / (1. - exp(-cwp(k)/2.))
         tmpra(k) = 1./(vkc*vkc*ur(k)) 
         tmpram(k) = log( (hgcm(k)-zpd(k)) / z0m(k) ) 
         tmprah(k) = log( (hgcm(k)-zpd(k)) / z0h(k) ) 

      end do

* -----------------------------------------------------------------
* coefficients for ir = air + bir * ts**4
* coefficients for st = ast + bst * ts
* -----------------------------------------------------------------

      do k = 1, npt
         if (veg(k)) then
            air(k) = -emv(k) * (1.+(1.-emv(k))*(1.-emg(k))) * firgcm(k) 
     &               -emv(k)*emg(k)*sb * tg(k)**4
            bir(k) = (2.-emv(k)*(1.-emg(k))) * emv(k) * sb
            ast(k) = 0.
            bst(k) = 0.
         else
            air(k) = -emg(k)*firgcm(k) 
            bir(k) = emg(k)*sb
            bst(k) = 2.*tksur(k)/dzsur(k)
            ast(k) = -bst(k)*tsoi(1,k)
         end if
      end do

* -----------------------------------------------------------------
* begin stability iteration
* -----------------------------------------------------------------

      do iter = 1, niters

* -----------------------------------------------------------------
* evaluate stability-dependent variables using moz from prior iteration
* -----------------------------------------------------------------

         do k = 1, npt

* adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.

            if (moz(k) .lt. 0.) then
               tmp1 = (1. - 16.*moz(k))**0.25
               tmp2 = log((1.+tmp1*tmp1)/2.)
               tmp3 = log((1.+tmp1)/2.)
               fmnew = 2.*tmp3 + tmp2 - 2.*atan(tmp1) + 1.5707963
               fhnew = 2*tmp2
            else
               fmnew = -5.*moz(k)
               fhnew = fmnew
            endif

* except for first iteration, weight stability factors for previous 
* iteration to help avoid flip-flops from one iteration to the next

            if (iter .eq. 1) then
               fm(k) = fmnew
               fh(k) = fhnew
            else
               fm(k) = 0.5 * (fm(k)+fmnew)
               fh(k) = 0.5 * (fh(k)+fhnew)
            endif

* aerodynamic resistances ram, rah and raw between heights hgcm and 
* zpd+z0m, zpd+z0h, and zpd+z0w

            tmp4 = tmpram(k)-fm(k)
            tmp5 = tmprah(k)-fh(k)
            ram(k,1) = max( tmpra(k)*tmp4*tmp4, 1.)
            rah(k,1) = max( tmpra(k)*tmp4*tmp5, 1.)
            raw(k,1) = rah(k,1)

* exchange coefficients

            cm(k) = 1./(ram(k,1)*ur(k))
            ch(k) = 1./(rah(k,1)*ur(k))
            cw(k) = 1./(raw(k,1)*ur(k))

* friction velocity

            fv(k) = ur(k) * sqrt(cm(k))

* aerodynamic resistances raw and rah between heights zpd+z0h and z0hg.
* if no vegetation, rah(2)=0 because zpd+z0h = z0hg.

            kh = max ( vkc*fv(k)*(hcan(k)-zpd(k)), mpe )
            ram(k,2) = 0.
            rah(k,2) = tmprah2(k) / kh
            raw(k,2) = rah(k,2)

* bulk boundary layer resistance

            uc = max ( fv(k)*tmpuc(k), mpe )
            rb(k) = tmprb(k) * sqrt(dleaf(ivt(k))/uc)

         end do

* -----------------------------------------------------------------
* svp and d(svp)/dt evaluated at ts
* -----------------------------------------------------------------

         do k = 1, npt
            t = tdc(ts(k),tfrz) 
            if (t .gt. 0.) then
               svpts(k) = esatw(t)
               dsvpts(k) = dsatw(t)
            else
               svpts(k) = esati(t)
               dsvpts(k) = dsati(t)
            end if
         end do

         if (iter .eq. 1) then

* -----------------------------------------------------------------
* soil surface resistance. should do each iteration because raw
* changes with stability. but not really necessary. set minimum
* rsurf so that in bare ground case es = svp(tg)
* -----------------------------------------------------------------

            do k = 1, npt
               if (veg(k)) then
                  rsurf(k) = 150.*fsno(k) + (1.-fsno(k))*raw(k,2) *
     &                       (1.-bevap(k))/max(bevap(k),mpe) 
               else
                  rsurf(k) = 150.*fsno(k) + (1.-fsno(k))*raw(k,1) *
     &                       (1.-bevap(k))/max(bevap(k),mpe)
               end if
               rsurf(k) = max( rsurf(k), mpe )
            end do

* -----------------------------------------------------------------
* stomatal resistances for sunlit and shaded fractions of canopy. 
* should do each iteration to account for differences in eah, tv,
* and rb but computationally expensive and not really necessary
* -----------------------------------------------------------------

            call stomata (npt  ,ivt    ,tfrz ,mpe    ,parsun,
     &                    tv   ,svpts  ,eah  ,tgcm   ,pgcm  ,
     &                    o2gcm,co2gcm ,igs  ,btran  ,rb    ,
     &                    rssun,psnsun ,qe25 ,aqe    ,kc25  ,
     &                    ko25 ,vcmx25 ,akc  ,ako    ,avcmx ,
     &                    bp   ,mp     ,foln ,folnmx ,c3psn )

            call stomata (npt  ,ivt    ,tfrz ,mpe    ,parsha,
     &                    tv   ,svpts  ,eah  ,tgcm   ,pgcm  ,
     &                    o2gcm,co2gcm ,igs  ,btran  ,rb    ,
     &                    rssha,psnsha ,qe25 ,aqe    ,kc25  ,
     &                    ko25 ,vcmx25 ,akc  ,ako    ,avcmx ,
     &                    bp   ,mp     ,foln ,folnmx ,c3psn )
         end if

* -----------------------------------------------------------------
* coefficients to calculate eah, tah, ev, tr, sh as functions of ts=tv 
* or ts=tg. eah and tah are at a height zpd+z0h as per calculation of 
* rah(1) and raw(1)
* -----------------------------------------------------------------

* tah(k) = ata(k) + bta(k)*ts(k) 
*  sh(k) = ash(k) + bsh(k)*ts(k) 

         do k = 1, npt
            if (veg(k)) then
               cah = 1./rah(k,1) 
               cvh = 2.*vai(k)/rb(k) 
               cgh = 1./rah(k,2) 
               cond = cah + cvh + cgh 
               ata(k) = (thgcm(k)*cah + tg(k)*cgh) / cond
               bta(k) = cvh/cond
               tmp1 = rhogcm(k)*cpair*cvh
               ash(k) = -ata(k)*tmp1
               bsh(k) = (1.-bta(k))*tmp1
            else
               cah = 1./rah(k,1) 
               cvh = 0. 
               cgh = 0.
               cond = cah + cvh + cgh 
               ata(k) = 0.
               bta(k) = 1.
               tmp1 = rhogcm(k)*cpair*cah
               ash(k) = -thgcm(k)*tmp1
               bsh(k) = tmp1
            end if
         end do

* eah(k) = aea(k) + bea(k)*esat[ts(k)] 
*  ev(k) = aev(k) + bev(k)*esat[ts(k)] 
*  tr(k) = atr(k) + btr(k)*esat[ts(k)]

         do k = 1, npt
            tmp1 = rhogcm(k)*cpair/gamma(k) 
            if (veg(k)) then
               caw = 1./raw(k,1) 
               cew = fwet(k)*vai(k)/rb(k)
               ctw = (1.-fwet(k)) * (laisun(k)/(rb(k)+rssun(k)) +
     &                               laisha(k)/(rb(k)+rssha(k)))
               cgw = 1./(raw(k,2)+rsurf(k)) 
               cond = caw + cew + ctw + cgw 
               aea(k) = (egcm(k)*caw + svptg(k)*cgw) / cond
               bea(k) = (cew+ctw)/cond
               aev(k) = -aea(k)*tmp1*cew 
               bev(k) = (1.-bea(k))*tmp1*cew
               atr(k) = -aea(k)*tmp1*ctw
               btr(k) = (1.-bea(k))*tmp1*ctw
            else
               caw = 1./raw(k,1)
               cew = 0. 
               ctw = 0.
               cgw = 1./rsurf(k)
               cond = caw + cew + ctw + cgw 
               aea(k) = egcm(k)*caw/cond
               bea(k) = cgw/cond
               aev(k) = -tmp1*(egcm(k)-aea(k))*caw
               bev(k) = tmp1*bea(k)*caw
               atr(k) = 0.
               btr(k) = 0.
            end if
            gcan(k) = cew + ctw
         end do

* -----------------------------------------------------------------
* evaluate surface fluxes with current temperature and solve for dts
* -----------------------------------------------------------------

*the following directive is inserted only to correct a bug on
*the j9's that cause the system to crash with a floating
*point exception. when this bug is fixed this directive will be
*removed.

         do k = 1, npt
            if (veg(k)) then
               i = 1
            else
               i = 2
            end if
            ir(k,i) = air(k) + bir(k)*ts(k)**4
            sh(k,i) = ash(k) + bsh(k)*ts(k)
            ev(k,i) = aev(k) + bev(k)*svpts(k)
            tr(k,i) = atr(k) + btr(k)*svpts(k)
            st(k,i) = ast(k) + bst(k)*ts(k)
            b = sa(k,i)-ir(k,i)-sh(k,i)-ev(k,i)-tr(k,i)-st(k,i)
            a = 4 .*bir(k)*ts(k)**3 + bsh(k) +
     &          (bev(k)+btr(k))*dsvpts(k) + bst(k)
            dts(k) = b/a
            ir(k,i) = ir(k,i) + 4.*bir(k)*ts(k)**3*dts(k) 
            sh(k,i) = sh(k,i) + bsh(k)*dts(k)
            ev(k,i) = ev(k,i) + bev(k)*dsvpts(k)*dts(k) 
            tr(k,i) = tr(k,i) + btr(k)*dsvpts(k)*dts(k) 
            st(k,i) = st(k,i) + bst(k)*dts(k)  
         end do

* -----------------------------------------------------------------
* update vegetation/ground surface temperature, canopy air temperature, 
* canopy vapor pressure, aerodynamic temperature, and
* monin-obukhov stability parameter moz for next iteration
* -----------------------------------------------------------------

         do k = 1, npt
            ts(k) = ts(k) + dts(k)
            tah(k) = ata(k) + bta(k)*ts(k) 
            t = tdc(ts(k),tfrz) 
            if (t .gt. 0.) then
               eah(k) = aea(k) + bea(k)*esatw(t) 
            else
               eah(k) = aea(k) + bea(k)*esati(t) 
            end if
            h = -(thgcm(k)-tah(k)) * rhogcm(k)*cpair/rah(k,1)
            tam(k) = tah(k) - h/(rhogcm(k)*cpair*fv(k)) * binv(k)

* monin-obukhov stability parameter moz for next iteration.
* set denominator of mol to minimum value if h and e equal 0.

            e = -(egcm(k)-eah(k)) * rhogcm(k)*cpair/(gamma(k)*raw(k,1))
            hvir = h + 0.61*cpair*tgcm(k) * e/lathea(k)
            tvir = (1. + 0.61*qgcm(k)) * tgcm(k)
            tmp1 = vkc * (grav/tvir) * hvir/(rhogcm(k)*cpair) 
            if (abs(tmp1) .le. mpe) tmp1 = mpe
            mol(k) = -1. * fv(k)**3 / tmp1
            if (pergro) then
               moz(k) = 0.
            else
               moz(k) = min( (hgcm(k)-zpd(k))/mol(k), 1.)
            end if

* accumulate number of times moz changes sign. set moz=0 if moz 
* changes sign 4 times so that fh=0 and fm=0. need to also reset
* fm and fh from previous iterations.

            if (mozold(k)*moz(k) .lt. 0.) nmozsgn(k) = nmozsgn(k)+1
            if (nmozsgn(k) .ge. 4) then
               moz(k) = 0.
               fm(k) = 0.
               fh(k) = 0.
            end if
            mozold(k) = moz(k)
         end do

      end do

* -----------------------------------------------------------------
* end stability iteration
* -----------------------------------------------------------------

* fill in vegetation and ground temperatures and delta temperatures 
* from ts and dts. if vegetated: tg and dtg will be set to correct
* values later. if not vegetated: set tv and dtv to ground values.

* also: surface can switch from vegetated to non-vegetated if buried 
* by snow. need to set tv=tg for non-vegetated surfaces so that tv is
* available if vegetation is exposed on next time step

      do k = 1, npt
         tv(k) = ts(k)
         dtv(k) = dts(k)
         tg(k) = ts(k)
         dtg(k) = dts(k)
      end do

* -----------------------------------------------------------------
* ground fluxes and tg if vegetated
* -----------------------------------------------------------------

      do k = 1, npt
         if (veg(k)) then
            air(k) = -emg(k)*(1.-emv(k))*firgcm(k) 
     &               -emg(k)*emv(k)*sb*tv(k)**4
            bir(k) = emg(k)*sb
            bsh(k) = rhogcm(k)*cpair/rah(k,2)
            ash(k) = -tah(k)*bsh(k)
            bev(k) = rhogcm(k)*cpair / (gamma(k)*(raw(k,2)+rsurf(k)))
            aev(k) = -eah(k)*bev(k)
            bst(k) = 2.*tksur(k)/dzsur(k)
            ast(k) = -bst(k)*tsoi(1,k)
         end if
      end do

      do iter = 1, niterg
         do k = 1, npt
            if (veg(k)) then
               t = tdc(tg(k),tfrz) 
               if (t. gt. 0.) then
                  svptg(k) = esatw(t)
                  dsvptg(k) = dsatw(t)
               else
                  svptg(k) = esati(t)
                  dsvptg(k) = dsati(t)
               end if
               ir(k,2) = air(k) + bir(k)*tg(k)**4
               sh(k,2) = ash(k) + bsh(k)*tg(k)
               ev(k,2) = aev(k) + bev(k)*svptg(k)
               tr(k,2) = 0.
               st(k,2) = ast(k) + bst(k)*tg(k)
               b = sa(k,2)-ir(k,2)-sh(k,2)-ev(k,2)-tr(k,2)-st(k,2)
               a = 4.*bir(k)*tg(k)**3+bsh(k)+bev(k)*dsvptg(k)+bst(k)
               dtg(k) = b/a
               if (iter .lt. niterg) tg(k) = tg(k)+dtg(k)
            end if
         end do
      end do

      do k=1,npt
         if (veg(k)) then
            ir(k,2) = ir(k,2) + 4.*bir(k)*tg(k)**3*dtg(k) 
            sh(k,2) = sh(k,2) + bsh(k)*dtg(k)
            ev(k,2) = ev(k,2) + bev(k)*dsvptg(k)*dtg(k) 
            st(k,2) = st(k,2) + bst(k)*dtg(k)  
            tg(k) = tg(k) + dtg(k)
         end if
      end do

* -----------------------------------------------------------------
* if snow on ground and tg > tfrz: reset tg = tfrz. reevaluate ground fluxes.
* energy inbalance used to melt snow. h2osno > 0.5 prevents spurious fluxes
* -----------------------------------------------------------------

      do k= 1, npt
         if (h2osno(k) .gt. 0.5 .and. tg(k) .gt. tfrz) then
            tg(k) = tfrz
            ir(k,2) = air(k) + bir(k)*tfrz**4
            sh(k,2) = ash(k) + bsh(k)*tfrz
            ev(k,2) = aev(k) + bev(k)*610.918
            tr(k,2) = 0.
            st(k,2) = sa(k,2) - (ir(k,2)+sh(k,2)+ev(k,2)+tr(k,2))
         end if
         if (h2osno(k) .gt. 0. .and. tg(k) .ge. tfrz) then           
            fsm(k) = min( h2osno(k)*hfus/dtim, max(st(k,2),0.) )
         else
            fsm(k) = 0.
         end if
      end do

* -----------------------------------------------------------------
* total fluxes. radiative temperature. wind stresses. 
* -----------------------------------------------------------------

      do k = 1, npt
         fira(k) = ir(k,1) + ir(k,2) 
         fsh(k) = sh(k,1) + sh(k,2) 
         fcev(k) = ev(k,1) 
         fctr(k) = tr(k,1) 
         fgev(k) = ev(k,2) 
         fcst(k) = st(k,1)
         fgr(k) = st(k,2) - fsm(k)
         fire(k) = firgcm(k) + fira(k) 
         trad(k) = (fire(k)/sb)**0.25
         taux(k) = -rhogcm(k)*cm(k)*ur(k)*ugcm(k)
         tauy(k) = -rhogcm(k)*cm(k)*ur(k)*vgcm(k)
      end do

* -----------------------------------------------------------------
* 2 m height air temperature
* -----------------------------------------------------------------

      do k = 1, npt
         tsa(k) = tah(k) - fsh(k)/(rhogcm(k)*cpair*fv(k)) * 
     &            1./vkc * log((2.+z0h(k))/z0h(k))
      end do

      return
      end
 
