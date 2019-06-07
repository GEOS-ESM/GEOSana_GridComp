#include <preproc.h>
      subroutine surrad (nband  ,npt    ,laisun ,fsa   ,sav    ,sag    ,
     &                   parsun ,mpe    ,solad  ,solai ,albgrd ,albgri ,
     &                   fabd   ,fabi   ,ftdd   ,ftid  ,ftii   ,fsr    ,
     &                   ndvi   ,albd   ,albi   ,fsun  ,fsha   ,elai   ,
     &                   vai    ,parsha ,laisha )

#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       surrad.F
* purpose:           solar fluxes absorbed by vegetation and ground surface
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
* input
      integer nband           !number of solar radiation waveband classes
      integer npt             !number of points

      real mpe                !prevents underflow errors if division by zero

      real fsun(npt)          !sunlit fraction of canopy
      real fsha(npt)          !shaded fraction of canopy
      real elai(npt)          !leaf area, one-sided
      real vai(npt)           !leaf + stem area, one-sided
      real laisun(npt)        !sunlit leaf area index, one-sided
      real laisha(npt)        !shaded leaf area index, one-sided

      real solad(mband,npt)   !incoming direct solar radiation (w/m**2) 
      real solai(mband,npt)   !incoming diffuse solar radiation (w/m**2) 
      real fabd(mband,npt)    !flux abs by veg (per unit incoming direct flux) 
      real fabi(mband,npt)    !flux abs by veg (per unit incoming diffuse flux)
      real ftdd(mband,npt)    !down dir flux below veg (per incoming dir flux) 
      real ftid(mband,npt)    !down dif flux below veg (per incoming dir flux) 
      real ftii(mband,npt)    !down dif flux below veg (per incoming dif flux) 
      real albgrd(mband,npt)  !ground albedo (direct)
      real albgri(mband,npt)  !ground albedo (diffuse)
      real albd(mband,npt)    !overall surface albedo (direct)
      real albi(mband,npt)    !overall surface albedo (diffuse)

* output
      real parsun(npt)        !average absorbed par for sunlit leaves (w/m**2)
      real parsha(npt)        !average absorbed par for shaded leaves (w/m**2)
      real sav(npt)           !solar radiation absorbed by vegetation (w/m**2)
      real sag(npt)           !solar radiation absorbed by ground (w/m**2)
      real fsa(npt)           !total absorbed solar radiation (w/m**2)
      real fsr(npt)           !total reflected solar radiation (w/m**2)
      real ndvi(npt)          !normalized difference vegetation index
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer i               !loop counter/array index
      integer ib              !waveband number (1=vis, 2=nir)

      real abs                !absorbed solar radiation (w/m**2) 
      real rnir               !reflected solar radiation [nir] (w/m**2)
      real rvis               !reflected solar radiation [vis] (w/m**2)
      real laifra             !leaf area fraction of canopy
      real trd                !transmitted solar radiation: direct (w/m**2)
      real tri                !transmitted solar radiation: diffuse (w/m**2)
      real cad(mband,npt)     !direct beam absorbed by canopy (w/m**2)
      real cai(mband,npt)     !diffuse radiation absorbed by canopy (w/m**2)

* -----------------------------------------------------------------

* zero summed solar fluxes

      do i = 1, npt
         sag(i) = 0.
         sav(i) = 0.
         fsa(i) = 0.
      end do

* loop over nband wavebands

      do ib = 1, nband
         do i = 1, npt

* absorbed by canopy

            cad(ib,i) = solad(ib,i)*fabd(ib,i)
            cai(ib,i) = solai(ib,i)*fabi(ib,i)
            sav(i) = sav(i) + cad(ib,i) + cai(ib,i)
            fsa(i) = fsa(i) + cad(ib,i) + cai(ib,i)

* transmitted = solar fluxes incident on ground

            trd = solad(ib,i)*ftdd(ib,i)
            tri = solad(ib,i)*ftid(ib,i) + solai(ib,i)*ftii(ib,i)

* solar radiation absorbed by ground surface

            abs = trd*(1.-albgrd(ib,i)) + tri*(1.-albgri(ib,i)) 
            sag(i) = sag(i) + abs
            fsa(i) = fsa(i) + abs

         end do
      end do

      do i = 1, npt

* partion visible canopy absorption to sunlit and shaded fractions
* to get average absorbed par for sunlit and shaded leaves

         laifra = elai(i) / max(vai(i),mpe)
         if (fsun(i) .gt. 0.) then
            parsun(i) = (cad(1,i)+fsun(i)*cai(1,i)) * laifra / 
     &                  max(laisun(i),mpe)
            parsha(i) = (fsha(i)*cai(1,i))*laifra / max(laisha(i),mpe)
         else
            parsun(i) = 0. 
            parsha(i) = (cad(1,i)+cai(1,i))*laifra /max(laisha(i),mpe)
         endif

* ndvi and reflected solar radiation

         rvis = albd(1,i)*solad(1,i) + albi(1,i)*solai(1,i) 
         rnir = albd(2,i)*solad(2,i) + albi(2,i)*solai(2,i)
         fsr(i) = rvis + rnir
         ndvi(i) = (rnir-rvis) / max(rnir+rvis,mpe)
      end do

      return
      end
 
