#include <preproc.h>
      subroutine surbgc (npt    ,ivt    ,igs    ,tv     ,tsoi  ,
     &                   psnsun ,laisun ,psnsha ,laisha ,tlai  ,
     &                   fpsn   ,frmf   ,frms   ,frmr   ,frm   ,
     &                   frg    ,fco2   ,dmi    ,fmicr  ,rsw   ,
     &                   btran  ,foln   ,soilc  ,stemb  ,rootb ,
     &                   nstep  ,kpti   ,loopi  )

#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       surbgc.F
* purpose:           surface biogeochemical fluxes. currently respiration 
*                    co2 fluxes and plant production
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
* input 
      integer nstep      !time step index
      integer npt        !number of points
      integer loopi      !"little" vector index (1 -> numlv)
      integer kpti       !first point in "big" kpt vec for loopi "little" vec
      integer ivt(npt)   !vegetation physiology type

      real igs(npt)      !growing season index (0=off, 1=on)
      real rsw(npt)      !root zone soil water
      real btran(npt)    !transpiration factor (0 to 1)
      real tlai(npt)     !total leaf area, unadjusted for burying by snow
      real tv(npt)       !vegetation temperature (kelvin)
      real psnsun(npt)   !sunlit leaf photosynthesis (umol co2 /m**2/ s) 
      real psnsha(npt)   !shaded leaf photosynthesis (umol co2 /m**2/ s)
      real laisun(npt)   !sunlit leaf area
      real laisha(npt)   !shaded leaf area
      real foln(npt)     !foliage nitrogen (%)
      real stemb(npt)    !stem biomass (kg /m**2)
      real rootb(npt)    !root biomass (kg /m**2)
      real soilc(npt)    !soil carbon (kg c /m**2)

      real tsoi(msl,npt) !soil temperature (kelvin)

* output
      real fpsn(npt)     !photosynthesis (umol co2 /m**2 /s)
      real frmf(npt)     !foliage maintenance respiration (umol co2 /m**2 /s)
      real frms(npt)     !stem maintenance respiration (umol co2 /m**2 /s)
      real frmr(npt)     !root maintenance respiration (umol co2 /m**2 /s)
      real frm(npt)      !total maintenance respiration (umol co2 /m**2 /s)
      real frg(npt)      !growth respiration (umol co2 /m**2 /s)
      real fmicr(npt)    !microbial respiration (umol co2 /m**2 /s)
      real fco2(npt)     !total co2 exchange (umol co2 /m**2 /s) [+ = to atm]
      real dmi(npt)      !total dry matter production (ug /m**2 /s)
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <vegcon.h>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      real fnf           !foliage nitrogen adjustemt to respiration (<= 1)
      real tf            !temperature factor
      real fsw           !soil water factor for microbial respiration
      real fst           !soil temperature factor for microbial respiration
      real rf            !respiration reduction factor (<= 1) 

      integer i,k        !array/loop index
* -----------------------------------------------------------------

      do k = 1, npt
         i = ivt(k)

* total photosynthesis

         fpsn(k) = psnsun(k)*laisun(k) + psnsha(k)*laisha(k)

* maintenance respiration

         if (igs(k) .eq. 0.) then
            rf = 0.5
         else
            rf = 1.0
         endif
         fnf = min( foln(k)/max(1.e-06,folnmx(i)), 1.0 )
         tf = arm(i)**( (tv(k)-298.16)/10. )
         frmf(k) = rmf25(i) * tf * fnf * tlai(k) * rf * btran(k)
         frms(k) = rms25(i) * stemb(k) * tf * rf
         frmr(k) = rmr25(i) * rootb(k) * tf * rf
         frm(k) = frmf(k) + frms(k) + frmr(k)          

* growth respiration and production

         frg(k) = 0.25 * fpsn(k)                     
         dmi(k) = (fpsn(k) - frm(k) - frg(k)) * dmcf(i)     

* microbial respiration

         fsw = rsw(k) / (0.20+rsw(k)) * 0.23 / (0.23+rsw(k))
         fst = 2.0**( (tsoi(1,k)-283.16)/10. )
         fmicr(k) = fsw * fst * mrp(i)*soilc(k)

* net co2 flux

         fco2(k) = -fpsn(k) + frm(k) + frg(k) + fmicr(k)

      end do
 
      return
      end
 
