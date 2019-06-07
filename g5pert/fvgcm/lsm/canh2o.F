#include <preproc.h>
      subroutine canh2o (npt     ,ivt     ,elai    ,esai    ,dtlsm   ,
     &                   hvap    ,hsub    ,tfrz    ,ch2op   ,hydro   ,
     &                   tgcm    ,tv      ,fcev    ,fctr    ,qprecc  ,
     &                   qprecl  ,h2ocan  ,qdrip   ,qdripl  ,qdripcl ,
     &                   qintr   ,qceva   ,qcdew   ,qtran   ,frc     ,
     &                   qrain   ,qsnow   )

#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ------------------------------
* source file:       canh2o.F
* purpose:           canopy hydrology 
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* --------------------------------------------------------------------

* ------------------------ input/output variables --------------------
* input
      integer npt         !number of points
      integer ivt(npt)    !plant type

      real dtlsm          !main time step (s)
      real hvap           !latent heat of vaporization (j/kg)
      real hsub           !latent heat of sublimation (j/kg)
      real tfrz           !freezing point (kelvin)
      real frc            !fraction of surface receiving convective precip
      real ch2op(mvt)     !maximum intercepted h2o per unit lai+sai (mm)

      logical hydro       !true if using prognostic hydrology

      real tgcm(npt)      !air temperature (kelvin)
      real tv(npt)        !vegetation temperature (kelvin)
      real fcev(npt)      !canopy evaporation (w/m**2) [+ = to atm, - = to srf]
      real fctr(npt)      !transpiration (w/m**2) [+ = to atm, - = to srf]
      real qprecc(npt)    !convective precipitation (mm h2o /s)
      real qprecl(npt)    !large-scale precipitation (mm h2o /s)
      real elai(npt)      !leaf area index, after burying by snow
      real esai(npt)      !stem area index, after burying by snow

* input/output
      real h2ocan(npt)    !intercepted water per unit ground area (mm h2o)

* output
      real qintr(npt)     !flux of water intercepted by canopy (mm h2o /s)
      real qdrip(npt)     !flux of water at ground (throughfall) (mm h2o /s) 
      real qdripl(npt)    !qdrip for large-scale region (mm h2o /s) 
      real qdripcl(npt)   !qdrip for convective + large-scale region (mm h2o /s)
      real qceva(npt)     !evaporation of intercepted water (mm h2o /s) [+]
      real qcdew(npt)     !canopy dew formation (mm h2o /s) [+]
      real qtran(npt)     !transpiration rate (mm h2o /s) [+]
      real qrain(npt)     !rain at ground srf, after interception (mm h2o/s) [+]
      real qsnow(npt)     !snow at ground srf, after interception (mm h2o/s) [+]
* --------------------------------------------------------------------

* ------------------------ local variables ---------------------------
      integer k         !do loop index

      real lathea       !latent heat of sublimation/vaporization (j/kg)
      real maxh2o(npt)  !maximum intercepted water (mm h2o / m**2)

      real prl          !precip over large-scale region (mm h2o /s)
      real prcl         !precip over convective + large-scale region (mm h2o /s)
      real inl          !interception for large-scale region (mm h2o /s)
      real incl         !interception for conv + large-scale region (mm h2o /s)
* --------------------------------------------------------------------

      do k = 1, npt

* maximum canopy water

         maxh2o(k) = ch2op(ivt(k)) * (elai(k)+ esai(k))

* interception for region 1-frc that receives only large scale precipitation

         prl = qprecl(k)                                  
         inl = min( 0.2*prl, (maxh2o(k)-h2ocan(k))/dtlsm ) 
         qdripl(k) = prl - inl                            

* interception for region frc that receives large scale + convective precip

         prcl = qprecc(k)/frc + qprecl(k)                     
         incl = min( 0.2*prcl, (maxh2o(k)-h2ocan(k))/dtlsm ) 
         qdripcl(k) = prcl - incl                           

* average interception and throughfall

         qintr(k) = inl*(1.-frc) + incl*frc
         qdrip(k) = qdripl(k)*(1.-frc) + qdripcl(k)*frc

* evaporation, transpiration, and dew

         if (tv(k) .gt. tfrz) then
            lathea = hvap
         else
            lathea = hsub
         endif
         qceva(k) = max( fcev(k)/lathea, 0. )
         qtran(k) = max( fctr(k)/lathea, 0. )
         qcdew(k) = abs( min( (fcev(k)+fctr(k))/lathea, 0. ) )

* canopy water balance. for convenience allow dew to bring h2ocan above 
* maxh2o or else would have to re-adjust qdrip

         if (hydro) then
            h2ocan(k) = h2ocan(k) + (qintr(k)+qcdew(k)-qceva(k))*dtlsm
         end if

* partition throughfall into rain or snow. if snow, no liquid water 
* at surface for infiltration (qdripl = qdripcl = 0)

         if (tgcm(k) .gt. tfrz+2.2) then
            qrain(k) = qdrip(k)
            qsnow(k) = 0.
         else
            qsnow(k) = qdrip(k)
            qrain(k) = 0.
            qdripl(k) = 0. 
            qdripcl(k) = 0.
         endif

      end do

      return
      end
 
