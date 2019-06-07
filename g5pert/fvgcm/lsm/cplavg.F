#include <preproc.h>
      subroutine cplavg (taux, tauy, lhflx , shflx,
     $                   lwup, qflx)

#include <vartyp.h>
#include <lsmpar.h>
#include <basinpar.h>

* ------------------------ code history --------------------------------
* source file:       cplavg.F
* purpose:           average lsm output fluxes for flux coupler
* date last revised: September 1996 - lsm version 1
* author:            Mariana Vertenstein
* standardized:      
* reviewed:         
* ----------------------------------------------------------------------

* ------------------------ notes ---------------------------------------
* Add land surface model output fluxes to accumulators every time step.
* When icnt==ncnt, compute the average flux over the time interval.
* ----------------------------------------------------------------------

* --------------- common block variables -------------------------------
#include <lsmtc.h>
#include <cplavg.h>
#include <basin.h>
* ----------------------------------------------------------------------

* ------------------------ input/output variables ----------------------
      real taux (lsmlon,lsmlat)   ! accumulator for land flux
      real tauy (lsmlon,lsmlat)   ! accumulator for land flux
      real lhflx(lsmlon,lsmlat)   ! accumulator for land flux
      real shflx(lsmlon,lsmlat)   ! accumulator for land flux
      real lwup (lsmlon,lsmlat)   ! accumulator for land flux
      real qflx (lsmlon,lsmlat)   ! accumulator for land flux
* ----------------------------------------------------------------------

* ------------------------ local variables -----------------------------
      integer i,lat               ! longitude,level,latitude indices
      integer n                   ! basin index
      real rncnt                  ! reciprocal of ncnt
* ----------------------------------------------------------------------
#if ( ! defined COUP_CCM )
* If ncnt == 1, then no averaging is required

      if (ncnt .eq. 1) return

      rncnt = 1./ncnt

!$omp  parallel do default(shared) private(i,lat)

      do lat =1, lsmlat

* Initial call of averaging interval, copy data to accumulators

         if (icnt .eq. 1) then
            do i = 1, numlon(lat)
               tauxa (i,lat) = taux (i,lat)   
               tauya (i,lat) = tauy (i,lat)   
               lhflxa(i,lat) = lhflx(i,lat)   
               shflxa(i,lat) = shflx(i,lat)   
               lwupa (i,lat) = lwup (i,lat)  
               qflxa (i,lat) = qflx (i,lat)
            end do

* Final call of averaging interval, complete averaging and copy data back

         else if (icnt .eq. ncnt) then
            do i = 1, numlon(lat)
               taux (i,lat) = rncnt * (tauxa (i,lat) + taux (i,lat))
               tauy (i,lat) = rncnt * (tauya (i,lat) + tauy (i,lat))
               lhflx(i,lat) = rncnt * (lhflxa(i,lat) + lhflx(i,lat))
               shflx(i,lat) = rncnt * (shflxa(i,lat) + shflx(i,lat))
               lwup (i,lat) = rncnt * (lwupa (i,lat) + lwup (i,lat))
               qflx (i,lat) = rncnt * (qflxa (i,lat) + qflx (i,lat))
            end do

* Intermediate call, add data to accumulators

         else
            do i = 1, numlon(lat)
               tauxa (i,lat) = tauxa (i,lat) + taux (i,lat)
               tauya (i,lat) = tauya (i,lat) + tauy (i,lat)
               lhflxa(i,lat) = lhflxa(i,lat) + lhflx(i,lat)
               shflxa(i,lat) = shflxa(i,lat) + shflx(i,lat)
               lwupa (i,lat) = lwupa (i,lat) + lwup (i,lat)
               qflxa (i,lat) = qflxa (i,lat) + qflx (i,lat)
            end do
         end if
      end do

* Average basin runoff vector

      do n = 1,ndrn
         if (icnt .eq. 1) then
            drnveca(n) = drnvec(n)
         else if (icnt .eq. ncnt) then            
            drnvec(n)  = rncnt* (drnveca(n) + drnvec(n))
         else
            drnvec(n)  = drnveca(n) + drnvec(n)
         endif
      end do

#endif
      return
      end
 
