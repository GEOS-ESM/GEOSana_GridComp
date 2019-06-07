#include <preproc.h>
      subroutine infil (npt    ,ist    ,qrain  ,qmelt  ,qsdew  , 
     &                  bch    ,smpsat ,hksat  ,watsat ,dzsoi  , 
     &                  h2osoi ,qover  ,qinfl  )
 
#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ------------------------------
* source file:       infil.F
* purpose:           infiltration and surface runoff 
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* --------------------------------------------------------------------

* ------------------------ input/output variables --------------------
* input
      integer npt       !number of points 
      integer ist(npt)  !soil texture type

      real qrain(npt)   !rain at ground surface, after interception (mm h2o /s) 
      real qsdew(npt)   !ground surface dew formation (mm h2o /s)
      real qmelt(npt)   !snow melt (mm h2o/s) 
      real watsat(npt)  !volumetric soil water content at saturation (porosity)
      real hksat(npt)   !hydraulic conductivity at saturation (mm h2o /s) 
      real smpsat(npt)  !soil matrix potential at saturation (mm) 
      real bch(npt)     !clapp and hornberger "b"

      real dzsoi(msl,npt)   !soil layer thickness (m)
      real h2osoi(msl,npt)  !volumetric soil water content

* output
      real qover(npt)   !surface runoff (mm h2o /s)
      real qinfl(npt)   !infiltration (mm h2o /s)
* --------------------------------------------------------------------

* ------------------------ local variables ---------------------------
      integer k         !do loop/array index
      real qsurf        !baseline flux liquid water (down) at surface (mm h2o/s)
      real v            !- [d psi / ds |s=1] / dz 
      real s            !h2osoi/watsat 
      real sq           !qsurf > infiltration capacity for s > sq 
      real flx1         !runoff where s >= 1
      real flx2         !runoff where s >= sq and s < 1
      real flx3         !runoff where s >= 0 and s < sq
      real c1,c2,c3,c4,c5,c6,c7 !constants
* --------------------------------------------------------------------

         do k = 1, npt
            qsurf = qmelt(k) + qsdew(k)
            v = bch(k)*smpsat(k) / (0.5*dzsoi(1,k)*1000.)
            s = h2osoi(1,k)/watsat(k)
            c1 = qsurf - hksat(k)*(1.-v)
            c2 = hksat(k)*v
            sq = max( min(c1/c2, 1.), 0. )
            c3 = max( qrain(k), 1.e-20 )
            c4 = exp(-1./s)
            c5 = exp(-sq/s)
            c6 = exp( min(80., c1/c3) )           !upper limit prevents overflow
            c7 = exp( min(80.,-sq/s-sq*c2/c3) )   !upper limit prevents overflow
c           flx1 = (qrain(k)+qsurf) * c4          !old, too much dunne runoff
            if (qsurf .le. hksat(k)) then
               flx1 = qrain(k)*c4*exp( (qsurf-hksat(k))/c3 )
            else
               flx1 = (qrain(k)+qsurf-hksat(k)) * c4
            end if
            flx2 = (qrain(k)+c1)*(c5-c4) + c2*(c4*(1.+s)-c5*(sq+s))
            flx3 = qrain(k)**2 / (qrain(k)+c2*s) * c6 * (1.-c7)
            if (ist(k). eq. 1) then               !soil
               qover(k) = flx1+flx2+flx3
            else                                  !glaciers, lakes, wetlands
               qover(k) = 0.
            end if
            qinfl(k) = qrain(k) + qsurf - qover(k)
         end do       

      return
      end
 
