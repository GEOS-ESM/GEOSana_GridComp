#include <preproc.h>
      subroutine phenol (npt  ,kmo ,kda  ,ivt  ,lat ,
     &                   tai  ,gai ,tlai ,tsai )
 
#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       phenol.F
* purpose:           set leaf and stem areas based on day of year
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
* input
      integer npt       !number of points
      integer kmo       !month (1, ..., 12)
      integer kda       !day (1, ..., 31)
      integer ivt(npt)  !vegetation type

      real lat(npt)     !latitude [+ = NH, - = SH]

      real tai(mvt,12)  !monthly leaf area index + stem area index, one-sided
      real gai(mvt,12)  !monthly leaf area index, one-sided

* output
      real tlai(npt)    !leaf area index, one-sided
      real tsai(npt)    !stem area index, one-sided
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer k           !index
      integer ndayn       !days in current year since jan 1: 1, ..., 365
      integer ndays       !ndayn shifted 6 mon for SH
      integer it1,it2     !interpolation months
      integer ndaypm(12)  !days per month
      data ndaypm /31,28,31,30,31,30,31,31,30,31,30,31/
      save ndaypm

      real wt1,wt2        !interpolation weights
      real day            !current day of year
      real t              !current month (1.00, ..., 12.00) 
* -----------------------------------------------------------------

* ndayn = days in current year since jan 1: 1, ..., 365
* ndays = ndayn shifted 6 mon for SH: 1 -> 183; 183 -> 365; 184 -> 1; 365 -> 182
 
! SJL: note that Feb is still used as if no leap year
      ndayn = 0
      do k = 1, kmo
        ndayn = ndayn + ndaypm(k)
      end do
      ndayn = ndayn - ndaypm(kmo) + kda
      ndays = mod (ndayn-1+365/2, 365) + 1
 
      do k = 1, npt

* compute seasonal indices it1,it2,wt1,wt2

         if (lat(k) .ge. 0.) then
            day = ndayn 
         else
            day = ndays 
         end if
         t = 12. * (day-0.5)/365.
         it1 = t + 0.5
         it2 = it1 + 1
         wt1 = (it1+0.5) - t
         wt2 = 1.-wt1
         if (it1 .lt.  1) it1 = 12
         if (it2 .gt. 12) it2 = 1

* leaf and stem areas

         tlai(k) = wt1*gai(ivt(k),it1) + wt2*gai(ivt(k),it2)
         tsai(k) = wt1*tai(ivt(k),it1) + wt2*tai(ivt(k),it2) - tlai(k)

      end do

      return
      end
 
