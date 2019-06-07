#include <preproc.h>
      subroutine ecodyn (npt     ,kmo     ,kda     ,ivt     ,hsno    ,
     &                   lati    ,tv      ,htop    ,elai    ,esai    ,
     &                   tlai    ,tsai    ,igs     ,stemb   ,rootb   ,
     &                   foln    ,soilc   ,h2ocan  ,fwet    ,doalb   ,
     &                   nstep   ,kpti    ,loopi   )

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       ecodyn.F
* purpose:           ecosystem dynamics: phenology, vegetation, soil carbon
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* this subroutine calculates:
*    o leaf areas (tlai, elai)
*    o stem areas (tsai, esai)
*    o wetted fraction of canopy (fwet)
*    o height (htop)
*    o growing season (igs)
*    o stem biomass (stemb)
*    o root biomass (rootb)
*    o foliage nitrogen (foln)
*    o soil carbon (soilc)

* calling sequence:
*    ecodyn:               !ecosystem dynamics driver
*      -> phenol:          !leaf phenology
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
* input
      integer npt      !number of points
      integer nstep    !time step index
      integer loopi    !"little" vector index (1 -> numlv)
      integer kpti     !first point in "big" kpt vec for loopi "little" vec
      integer kmo      !month (1, ..., 12)
      integer kda      !day (1, ..., 31)
      integer ivt(npt) !vegetation type

      logical doalb    !true = surface albedo calculation time step

      real hsno(npt)   !snow depth (m)
      real lati(npt)   !latitude [+ = NH, - = SH]
      real tv(npt)     !vegetation temperature (kelvin)
      real h2ocan(npt) !intercepted water (mm h2o)

* output
      real elai(npt)   !leaf area index, one-sided, adjusted for burying by snow
      real esai(npt)   !stem area index, one-sided, adjusted for burying by snow
      real tlai(npt)   !leaf area index, one-sided, unadjust for burying by snow
      real tsai(npt)   !stem area index, one-sided, unadjust for burying by snow
      real htop(npt)   !top of canopy (m)
      real igs(npt)    !growing season index (0=off, 1=on)
      real stemb(npt)  !stem biomass (kg /m**2)
      real rootb(npt)  !root biomass (kg /m**2)
      real foln(npt)   !foliage nitrogen (%)
      real soilc(npt)  !soil carbon (kg c /m**2)
      real fwet(npt)   !wetted fraction of canopy (0 to 1)
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <vegcon.h>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer k        !do loop index
      real ol          !thickness of canopy layer covered by snow (m)
      real fb          !fraction of canopy layer covered by snow
      real maxh2o      !maximum intercepted water (mm h2o / m**2)
      real hbot(npt)   !bottom of canopy (m)
* -----------------------------------------------------------------

* -----------------------------------------------------------------
* define growing season
* -----------------------------------------------------------------

      do k = 1, npt
         if (tv(k) .gt. tmin(ivt(k))) then
            igs(k) = 1.
         else
            igs(k) = 0.
         end if
      end do

* -----------------------------------------------------------------
* need to update elai and esai only every albedo time step so do not 
* have any inconsistency in lai and sai between suralb calls (i.e., 
* if albedos are not done every time step).
* -----------------------------------------------------------------

      if (doalb) then

* -----------------------------------------------------------------
* leaf phenology
* -----------------------------------------------------------------

         call phenol (npt  ,kmo ,kda  ,ivt  ,lati ,
     &                tai  ,gai ,tlai ,tsai )

* -----------------------------------------------------------------
* vegetation dynamics
* -----------------------------------------------------------------

         do k = 1, npt
            htop(k)  = hvt(ivt(k))     !canopy top height
            hbot(k)  = hvb(ivt(k))     !canopy bottom height
            stemb(k) = stembvt(ivt(k)) !stem biomass
            rootb(k) = rootbvt(ivt(k)) !root biomass
            foln(k)  = folnvt(ivt(k))  !foliage nitrogen
         end do

* -----------------------------------------------------------------
* adjust lai and sai for burying by snow. if exposed lai and sai
* are less than 0.05, set equal to zero to prevent numerical 
* problems associated with very small lai and sai.
* -----------------------------------------------------------------

         do k = 1, npt
            ol = min( max(hsno(k)-hbot(k),0.), htop(k)-hbot(k))
            fb = 1. - ol / max(1.e-06,htop(k)-hbot(k))
            elai(k) = tlai(k)*fb
            esai(k) = tsai(k)*fb
            if (elai(k) .lt. 0.05) elai(k) = 0.
            if (esai(k) .lt. 0.05) esai(k) = 0.
         end do

      end if

* -----------------------------------------------------------------
* wetted fraction of canopy
* -----------------------------------------------------------------

      do k = 1, npt
         maxh2o = ch2op(ivt(k)) * (elai(k)+ esai(k))
         fwet(k) = max(h2ocan(k),0.) / max(maxh2o,1.e-06)
C SJL Fix for CRAY
         fwet(k) = max( fwet(k), 0.)    !ccm3.6.6
         fwet(k) = min( fwet(k), 1.) **0.667
      end do

* -----------------------------------------------------------------
* soil carbon dynamics
* -----------------------------------------------------------------

      do k = 1, npt
         soilc(k) = soilcvt(ivt(k))
      end do

      return
      end
 
