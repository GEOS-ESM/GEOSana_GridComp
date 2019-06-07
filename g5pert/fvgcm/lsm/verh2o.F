#include <preproc.h>
      subroutine verh2o (npt     ,ivt     ,elai    ,esai    ,dtsoi   ,
     &                   dtlsm   ,qdrip   ,qsoil   ,qvege   ,tgcm    ,
     &                   qinfl   ,tg      ,tv      ,qover   ,fsm     , 
     &                   h2osno  ,qmelt   ,hsno    ,ist     ,h2ocan  , 
     &                   fcev    ,fgev    ,fctr    ,h2osoi  ,qdrai   , 
     &                   dzsoi   ,loopi   ,root    ,watsat  ,hksat   , 
     &                   smpsat  ,bch     ,rsw     ,qintr   ,qprecl  ,
     &                   qprecc  ,errh2o  ,hydro   ,zsoi    ,nstep   ,
     &                   conchk  ,fsno    ,lakpoi  ,qvegt   ,kpti    )

#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ------------------------------
* source file:       verh2o.F
* purpose:           vertical (column) surface hydrology driver
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* --------------------------------------------------------------------

* ------------------------ notes -------------------------------------
* surface hydrology: all water fluxes and pools are per unit ground area.
* 1 kg h2o / m**2 = 1 mm h2o. processes simulated are:

*    o interception
*    o throughfall/stemflow
*    o snow accumulation and melt
*    o infiltration and surface runoff
*    o soil water and drainage
*    o irrigation

* the water vapor flux has three components:
*    o qvege: vegetation evaporation
*    o qvegt: vegetation transpiration
*    o qsoil: soil/snow  evaporation/sublimation

* water balances are:
*    o canopy water:
* [h2ocan] n+1 = [h2ocan] n + (qprecc + qprecl - qdrip + qcdew - qceva)*dtlsm

*    o snow water:
* [h2osno] n+1 = [h2osno] n + (qsnow + qfros - qsubl - qmelt)*dtlsm

*    o soil water:
* [h2osoi(1)*dzsoi(1)*1000 + ... + h2osoi(nsl)*dzsoi(nsl)*1000] n+1 =
* [h2osoi(1)*dzsoi(1)*1000 + ... + h2osoi(nsl)*dzsoi(nsl)*1000] n   +
* (qinfl - qseva - qtran - qdrai)*dtlsm

*    o total water:
* [h2ocan] n+1 + [h2osno] n+1 +
* [h2osoi(1)*dzsoi(1)*1000 + ... + h2osoi(nsl)*dzsoi(nsl)*1000] n+1 = 
* [h2ocan] n   + [h2osno] n   +
* [h2osoi(1)*dzsoi(1)*1000 + ... + h2osoi(nsl)*dzsoi(nsl)*1000] n   + 
* (qprecc + qprecl - qvege - qvegt - qsoil - qover - qdrai)*dtlsm

* water conservation check: valid for non-irrigated soil only
* because qover=qdrai=0 for non-soil surfaces,  no soil hydrology 
* for no soil (ist > 1), and some crops are irrigated

* the calling sequence is:
*  -> verh2o:           surface hydrology driver
*       -> canh2o:      canopy water
*       -> snoh2o:      snow water
*       -> infil:       infiltration for large-scale precipitation region
*       -> infil:       infiltration for convective + large-scale precip region
*       -> soih2o:      soil hydrology
*            -> tridia: tridiagonal matrix solution
* --------------------------------------------------------------------

* ------------------------ input/output variables --------------------
* input
      integer npt         !number of points
      integer ist(npt)    !soil "texture" type
      integer ivt(npt)    !plant type
      integer nstep       !time index
      integer loopi       !"little" vector index (1 -> numlv)
      integer kpti        !first point in "big" kpt vec for loopi "little" vec

      real dtlsm          !main time step (s)
      real dtsoi          !soil hydrology time step (s)

      logical hydro       !true if using prognostic hydrology
      logical conchk      !true if want conservation check
      logical lakpoi(npt) !true if lake point

      real tgcm(npt)      !air temperature (kelvin)
      real tv(npt)        !vegetation temperature (kelvin)
      real tg(npt)        !ground temperature (kelvin)
      real fcev(npt)      !canopy evaporation (w/m**2) [+ = to atm, - = to srf]
      real fctr(npt)      !transpiration (w/m**2) [+ = to atm, - = to srf]
      real fgev(npt)      !ground evaporation (w/m**2) [+ = to atm, - = to srf]
      real fsm(npt)       !snow melt heat flux (w/m**2) 
      real qprecc(npt)    !convective precipitation (mm h2o /s)
      real qprecl(npt)    !large-scale precipitation (mm h2o /s)
      real watsat(npt)    !volumetric soil water content, saturation (porosity)
      real hksat(npt)     !hydraulic conductivity at saturation (mm h2o /s) 
      real smpsat(npt)    !soil matrix potential at saturation (mm) 
      real bch(npt)       !clapp and hornberger "b"
      real elai(npt)      !leaf area index, after burying by snow
      real esai(npt)      !stem area index, after burying by snow

      real root(msl,npt)  !relative root abundance (0 to 1)
      real zsoi(msl,npt)  !soil layer depth (m)
      real dzsoi(msl,npt) !soil layer thickness (m)

* input/output
      real h2osno(npt)    !snow water per unit ground area (mm h2o)
      real h2ocan(npt)    !intercepted water per unit ground area (mm h2o)
      real h2osoi(msl,npt)!volumetric soil water content (0<=h2osoi<=watsat)

* output
      real qmelt(npt)     !snow melt (mm h2o /s) 
      real qinfl(npt)     !infiltration (mm h2o /s)
      real qsoil(npt)     !soil evaporation (mm h2o/s) [+ = to atm, - = to srf]
      real qvege(npt)     !vegetation evaporation (mm h2o /s) [+ = to atm]
      real qvegt(npt)     !vegetation transpiration (mm h2o /s) [+ = to atm]
      real qover(npt)     !surface runoff (mm h2o /s) 
      real qdrai(npt)     !sub-surface runoff (mm h2o /s) 
      real qintr(npt)     !flux of water intercepted by canopy (mm h2o /s)
      real qdrip(npt)     !flux of water at ground (throughfall) (mm h2o /s) 
      real hsno(npt)      !depth of snow pack (m)
      real rsw(npt)       !soil water content for root zone (up to 1 m depth)
      real errh2o(npt)    !total water conservation inbalance (mm h2o)
      real fsno(npt)      !fraction of ground covered by snow
* --------------------------------------------------------------------

* ------------------------ common block variables --------------------
#include <phycon.h>
#include <snocon.h>
#include <soicon.h>
#include <vegcon.h>
* --------------------------------------------------------------------

* ------------------------ local variables ---------------------------

      integer j,k         !do loop/array index
      integer nsl         !number of soil layers

      real lathea         !latent heat of sublimation/vaporization (j/kg)

      real qoverl(npt)    !qover for large-scale region (mm h2o /s)
      real qovercl(npt)   !qover for convective + large-scale region (mm h2o/s)
      real qinfll(npt)    !qinfl for large-scale region (mm h2o /s)
      real qinflcl(npt)   !qinfl for convective + large-scale region (mm h2o/s)
      real qdripl(npt)    !qdrip for large-scale region (mm h2o /s) 
      real qdripcl(npt)   !qdrip for convective + large-scale region (mm h2o/s)
      real owat(npt)      !total water at begining of time step (mm h2o)
      real nwat(npt)      !total water at end of time step (mm h2o)
      real qceva(npt)     !evaporation of intercepted water (mm h2o /s) [+]
      real qcdew(npt)     !canopy dew formation (mm h2o /s) [+]
      real qrain(npt)     !rain at ground srf, aftr interception (mm h2o/s) [+]
      real qsnow(npt)     !snow at ground srf, aftr interception (mm h2o/s) [+]
      real qsubl(npt)     !sublimation rate from snow pack (mm h2o /s) [+]
      real qseva(npt)     !ground surface evaporation rate (mm h2o/s) [+]
      real qtran(npt)     !transpiration rate (mm h2o /s) [+]
      real qsdew(npt)     !ground surface dew formation (mm h2o /s) [+]
      real qfros(npt)     !surface dew added to snow pack (mm h2o /s) [+]
      real rwat(npt)      !soil water to maximum depth of 1 m (m)
      real rz(npt)        !thickness of soil layers contributing to rwat (m)

      logical irrig(npt)  !true if soil is irrigated

      real frc            !fraction of surface receiving convective precip
      data frc /0.6/
      save frc
* --------------------------------------------------------------------

      nsl = msl

* --------------------------------------------------------------------
* tag irrigated soils. save total water for conservation check
* --------------------------------------------------------------------

      do k = 1, npt
         if (ivt(k) .eq. nic) then
            irrig(k) = .true.
         else
            irrig(k) = .false.
         end if
         owat(k) = h2ocan(k)+h2osno(k)
      end do

      do j = 1, nsl
         do k = 1, npt
            owat(k) = owat(k) + h2osoi(j,k)*dzsoi(j,k)*1000.
         end do
      end do

* --------------------------------------------------------------------
* canopy hydrology
* --------------------------------------------------------------------

      call canh2o (npt     ,ivt     ,elai    ,esai    ,dtlsm   ,
     &             hvap    ,hsub    ,tfrz    ,ch2op   ,hydro   ,
     &             tgcm    ,tv      ,fcev    ,fctr    ,qprecc  ,
     &             qprecl  ,h2ocan  ,qdrip   ,qdripl  ,qdripcl ,
     &             qintr   ,qceva   ,qcdew   ,qtran   ,frc     ,
     &             qrain   ,qsnow   )

* --------------------------------------------------------------------
* surface evaporation and dew
* --------------------------------------------------------------------

      do k = 1, npt
         if (tg(k) .gt. tfrz) then
            lathea = hvap
         else
            lathea = hsub
         endif
         qseva(k) = max( fgev(k)/lathea, 0. )
         qsdew(k) = abs( min( fgev(k)/lathea, 0. ) )
      end do

* --------------------------------------------------------------------
* snow hydrology
* --------------------------------------------------------------------

      call snoh2o (npt   ,dtlsm  ,hfus   ,tfrz   ,bdsno ,
     &             hsnoc ,hydro  ,lakpoi ,tg     ,fsm   ,
     &             qsnow ,h2osno ,qseva  ,qsdew  ,qmelt ,
     &             hsno  ,fsno   ,qsubl  ,qfros  )

      if (hydro) then

* --------------------------------------------------------------------
* infiltration and surface runoff 
* --------------------------------------------------------------------

* infiltration and surface runoff for large-scale and convective +
* large-scale fractions of surface

         call infil (npt    ,ist     ,qdripl  ,qmelt   ,qsdew  ,
     &               bch    ,smpsat  ,hksat   ,watsat  ,dzsoi  ,
     &               h2osoi ,qoverl  ,qinfll  )
         call infil (npt    ,ist     ,qdripcl ,qmelt   ,qsdew  , 
     &               bch    ,smpsat  ,hksat   ,watsat  ,dzsoi  ,
     &               h2osoi ,qovercl ,qinflcl )

* total infiltration and surface runoff

         do k = 1, npt
            qover(k) = qoverl(k)*(1.-frc) + qovercl(k)*frc
            qinfl(k) = qinfll(k)*(1.-frc) + qinflcl(k)*frc
         end do

* --------------------------------------------------------------------
* soil water and sub-surface drainage
* --------------------------------------------------------------------

         call soih2o (npt    ,nsl    ,dtlsm ,dtsoi ,ist    ,qtran ,
     &                qinfl  ,h2osoi ,root  ,dzsoi ,watsat ,bch   ,
     &                smpsat ,hksat  ,qdrai ,qseva )

* --------------------------------------------------------------------
* irrigate crops: to depth of 1 m
* --------------------------------------------------------------------

         do j = 1, nsl
            do k = 1, npt
               if (irrig(k) .and. (elai(k)+esai(k)) .gt. 0.) then
                  if (zsoi(j,k)+0.5*dzsoi(j,k) .le. 1.) then
                     h2osoi(j,k) = watsat(k)
                  end if
               end if
            end do
         end do

      end if

      if (.not. hydro) then

* --------------------------------------------------------------------
* no interactive hydrology
* --------------------------------------------------------------------

         do k = 1, npt
            qover(k) = 0.
            qinfl(k) = 0.
            qdrai(k) = 0.
         end do

      end if

* --------------------------------------------------------------------
* water fluxes to atmosphere. soil water content up to a depth of 1 m.
* --------------------------------------------------------------------

      do k = 1, npt
         qsoil(k) = qseva(k) + qsubl(k) - qsdew(k) - qfros(k)
         qvege(k) = qceva(k) - qcdew(k)
         qvegt(k) = qtran(k)
         rwat(k) = 0.
         rz(k) = 0.
      end do

      do j = 1, nsl
         do k = 1, npt
            if (zsoi(j,k)+0.5*dzsoi(j,k) .le. 1.) then
               rwat(k) = rwat(k) + h2osoi(j,k)*dzsoi(j,k)
               rz(k) = rz(k) + dzsoi(j,k)
            end if
         end do
      end do

      do k = 1, npt
         rsw(k) = rwat(k)/rz(k)
      end do

* --------------------------------------------------------------------
* water conservation check: valid for non-irrigated soil only
* because qover=qdrai=0 for non-soil surfaces and no soil hydrology 
* for non-soil surfaces
* --------------------------------------------------------------------

      do k = 1, npt
         nwat(k) = h2ocan(k)+h2osno(k)
      end do

      do j = 1, nsl
         do k = 1, npt
            nwat(k) = nwat(k) + h2osoi(j,k)*dzsoi(j,k)*1000.
         end do
      end do

      do k = 1, npt
         if (ist(k).eq.1 .and. .not.irrig(k)) then 
            errh2o(k) = nwat(k) - owat(k) - (qprecc(k)+qprecl(k) -
     &                  qvege(k)-qvegt(k)-qsoil(k)-qover(k)-qdrai(k)) *
     &                  dtlsm
         else
            errh2o(k) = 0.
         end if
      end do

* check for conservation of water. do only if hydro = true 

      if (conchk .and. hydro) then
         do k = 1, npt
           if (abs(errh2o(k)) .gt. 0.01) write (6,1001)
     &     nstep,loopi,k,k+kpti-1,errh2o(k)
         end do
1001     format (1x,'tot wat error:',' nstep =',i10,' vector =',i3,
     &   ' point =',i3,' lsm point =',i6,' inbalance =',f8.2,' mm')
      end if

      return
      end
 
