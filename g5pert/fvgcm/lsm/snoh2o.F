#include <preproc.h>
      subroutine snoh2o (npt   ,dtlsm  ,hfus   ,tfrz   ,bdsno ,   
     &                   hsnoc ,hydro  ,lakpoi ,tg     ,fsm   , 
     &                   qsnow ,h2osno ,qseva  ,qsdew  ,qmelt ,
     &                   hsno  ,fsno   ,qsubl  ,qfros  )

#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ------------------------------
* source file:       snoh2o.F
* purpose:           snow hydrology 
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* --------------------------------------------------------------------

* ------------------------ input/output variables --------------------
* input
      integer npt         !number of points

      real dtlsm          !main time step (s)
      real hfus           !latent heat of fusion (j/kg)
      real tfrz           !freezing point (kelvin)
      real bdsno          !bulk density snow (kg/m**3)
      real hsnoc          !height of snow when ground fully covered by snow (m)

      logical hydro       !true if using prognostic hydrology
      logical lakpoi(npt) !true if lake point

      real tg(npt)        !ground temperature (kelvin)
      real fsm(npt)       !snow melt heat flux (w/m**2) 
      real qsnow(npt)     !snow at ground srf, after interception (mm h2o/s) [+]

* input/output
      real h2osno(npt)    !snow water per unit ground area (mm h2o)
      real qseva(npt)     !ground srf evap, after sublimation (mm h2o/s) [+]
      real qsdew(npt)     !ground srf dew formation, after frost (mm h2o /s) [+]

* output
      real qmelt(npt)     !snow melt (mm h2o /s) 
      real hsno(npt)      !depth of snow pack (m)
      real fsno(npt)      !fraction of ground covered by snow
      real qsubl(npt)     !sublimation rate from snow pack (mm h2o /s) [+]
      real qfros(npt)     !surface dew added to snow pack (mm h2o /s) [+]
* --------------------------------------------------------------------

* ------------------------ local variables ---------------------------
      integer k           !do loop/array index
* --------------------------------------------------------------------

      do k = 1, npt

* snow melt

         qmelt(k) = fsm(k)/hfus

* sublimation. do not allow for more sublimation than there is snow
* after melt. remaining surface evaporation used for infiltration

         qsubl(k) = min( qseva(k), h2osno(k)/dtlsm-qmelt(k) )
         qseva(k) = qseva(k) - qsubl(k)

* if dew added to snow as frost, reduce surface dew accordingly
* for later use in net surface water flux for infiltration. 
* note: add dew as frost only if surface temperature less than
* approximately freezing. do not use tfrz exactly because if 
* small amount of snow on ground and melting, tg=tfrz. then
* water flux towards surface (i.e. dew) will increase snow.

         if (tg(k) .ge. (tfrz-0.1)) then
            qfros(k) = 0.
         else
            qfros(k) = qsdew(k)
         endif
         qsdew(k) = qsdew(k) - qfros(k)

* update snow pack

         if (hydro) then
            h2osno(k) = h2osno(k) +
     &                 (qsnow(k)-qmelt(k)-qsubl(k)+qfros(k))*dtlsm
         end if
         h2osno(k) = max( h2osno(k), 0. )

      end do

* no snow if lake unfrozen

      do k = 1, npt
         if (lakpoi(k) .and. tg(k).gt.tfrz) h2osno(k) = 0.
      end do

* snow height and fractional coverage

      do k = 1, npt
         hsno(k) = h2osno(k)/bdsno
         fsno(k) = min( hsno(k)/hsnoc, 1. )
      end do

      return
      end
 
