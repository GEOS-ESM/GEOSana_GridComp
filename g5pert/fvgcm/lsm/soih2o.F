#include <preproc.h>
      subroutine soih2o (npt    ,nsl    ,dtlsm ,dtsoi ,ist    ,qtran , 
     &                   qinfl  ,h2osoi ,root  ,dzsoi ,watsat ,bch   , 
     &                   smpsat ,hksat  ,qdrai ,qseva )
 
#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ------------------------------
* source file:       soih2o.F
* purpose:           soil hydrology and sub-surface drainage
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* --------------------------------------------------------------------

* ------------------------ notes -------------------------------------
* use tridiagonal system of equations to solve one-dimensional water balance:
 
*     d wat   
* dz  ----- = -qi + qo - s
*       dt    

* with q = -k d(psi+z)/dz = -k (d psi/dz + 1) and with s=0
* this is the Richards equation for vertical water flow

* d wat    d     d wat d psi
* ----- = -- [ k(----- ----- + 1) ]
*   dt    dz       dz  d wat

* where: wat = volume of water per volume of soil (mm**3/mm**3)
*        psi = soil matrix potential (mm)
*        dt  = time step (s)
*        z   = depth (mm)
*        dz  = thickness (mm)
*        qi  = inflow at top (mm h2o /s) (+ = up, - = down)
*        qo  = outflow at bottom (mm h2o /s) (+ = up, - = down)
*        s   = source/sink flux (mm h2o /s) (+ = loss, - = gain)
*        k   = hydraulic conductivity (mm h2o /s)

* solution: linearize k and psi about d wat and use tridiagonal system 
* of equations to solve for d wat, where for layer i 
* r_i = a_i [d wat_i-1] + b_i [d wat_i] + c_i [d wat_i+1] 

* the solution conserves water as:
* [h2osoi(1)*dzsoi(1)*1000 + ... + h2osoi(nsl)*dzsoi(nsl)*1000] n+1 =
* [h2osoi(1)*dzsoi(1)*1000 + ... + h2osoi(nsl)*dzsoi(nsl)*1000] n +
* (qinfl - qseva - qtran - qdrai)*dtlsm

* code only executed for soils (ist = 1)
* --------------------------------------------------------------------

* ------------------------ input/output variables --------------------
* input
      integer npt           !number of points
      integer nsl           !number of soil layers
      integer ist(npt)      !soil texture type

      real dtlsm            !main lsm time step (s)
      real dtsoi            !soil water time step <= dtlsm (s)

      real qtran(npt)       !transpiration water flux (mm h2o/s)
      real qinfl(npt)       !infiltration rate (mm h2o/s)
      real qseva(npt)       !ground surface evaporation rate (mm h2o/s)
      real watsat(npt)      !saturated volumetric soil water content (porosity)
      real smpsat(npt)      !soil matrix potential at saturation (mm)
      real bch(npt)         !clapp and hornberger "b" 
      real hksat(npt)       !hydraulic conductivity at saturation (mm h2o/s)

      real root(msl,npt)    !relative root abundance (0 to 1)
      real dzsoi(msl,npt)   !soil layer thickness (m)

* input/output
      real h2osoi(msl,npt)  !volumetric soil water content

* output
      real qdrai(npt)       !sub-surface runoff (mm h2o/s)
* --------------------------------------------------------------------

* ------------------------ local variables ---------------------------
      integer j,k,l         !do loop/array indices
      integer ncal          !number of points for soil water calculation
      integer ind(npt)      !index to points for soil water calculation
      integer iter          !number of iterations (dtlsm/dtsoi)

      real r(msl,npt)       !solution matrix
      real a(msl,npt)       !"a" vector for tridiagonal matrix
      real b(msl,npt)       !"b" vector for tridiagonal matrix
      real c(msl,npt)       !"c" vector for tridiagonal matrix
      real dwat(msl,npt)    !change in soil water
      real smp(msl,npt)     !soil matrix potential (mm)
      real hk(msl,npt)      !hydraulic conductivity (mm h2o/s)
      real hk2(msl,npt)     !hk**2
      real dsmpdw(msl,npt)  !d(smp)/d(h2osoi)
      real dhkdw(msl,npt)   !d(hk)/d(h2osoi)

      real s                !h2osoi/watsat
      real hydcon           !hydraulic conductivity (mm h2o/s)
      real qin              !flux of water into soil layer (mm h2o/s)
      real qout             !flux of water out of soil layer (mm h2o/s)
      real num              !used in calculating qi, qout
      real den              !used in calculating qi, qout
      real den2             !den**2 used in calculating qi, qout
      real dqidw1           !d(qin)/d(h2osoi(i-1))
      real dqidw2           !d(qin)/d(h2osoi(i))
      real dqodw1           !d(qout)/d(h2osoi(i))
      real dqodw2           !d(qout)/d(h2osoi(i+1))
      real xs(npt)          !soil water > sat or < some minimum (mm h2o)
      real axs(npt)         !amount of xs added to a soil layer
      real x                !temporary value of axs
      real newwat(npt)      !temporary value of updated h2osoi
      real watmin           !limit h2osoi >= watmin
* --------------------------------------------------------------------

* do soil water calculation for points with soil

      ncal = 0
      do k = 1, npt
         qdrai(k) = 0.
         if (ist(k) .eq. 1) then
            ncal = ncal+1
            ind(ncal) = k
         endif
      end do
      if (ncal .eq. 0) then
         return
      end if	

      do iter = 1, int(dtlsm/dtsoi)

* initialize xs

         do l = 1, ncal 
            k = ind(l)
            xs(k) = 0.
         end do

* evaluate hydraulic conductivity, soil matrix potential, 
* d(smp)/d(h2osoi), and d(hk)/d(h2osoi). limit s >= 0.05 
* when evaluating these terms. this helps prevent numerical
* problems for very small h2osoi. also limit hk >= some very
* small number for same reason.

         do j = 1, nsl
            do l = 1, ncal 
               k = ind(l)
               s = max(h2osoi(j,k)/watsat(k),0.05)
               smp(j,k) = smpsat(k) * s**(-bch(k))
               dsmpdw(j,k) = -bch(k)*smp(j,k)/h2osoi(j,k)
               hydcon = hksat(k) * s**(2.*bch(k)+3.)
               hk(j,k) = max(hydcon,1.e-10)
               dhkdw(j,k) = hydcon*(2.*bch(k)+3.)/h2osoi(j,k)
               hk2(j,k) = hk(j,k)*hk(j,k)
            end do
         end do

* set up r, a, b, and c vectors for tridiagonal solution

* node j=1

         j = 1
         do l = 1, ncal
            k = ind(l)
            num = -2.*(smp(j,k)-smp(j+1,k)) -
     &            1000.*(dzsoi(j,k)+dzsoi(j+1,k))
            den = 1000.*(dzsoi(j,k)/hk(j,k) + dzsoi(j+1,k)/hk(j+1,k))
            den2 = den*den
            qout = num/den
            dqodw1 = (-2.*den*dsmpdw(j  ,k) +
     &               num*1000.*dzsoi(j  ,k)/hk2(j  ,k)*dhkdw(j  ,k))
     &               / den2
            dqodw2 = ( 2.*den*dsmpdw(j+1,k) +
     &               num*1000.*dzsoi(j+1,k)/hk2(j+1,k)*dhkdw(j+1,k))
     &               / den2

            r(j,l) = (qseva(k)+qtran(k)*root(j,k)) - qinfl(k) - qout
            a(j,l) = 0.
            b(j,l) = dqodw1 - 1000.*dzsoi(j,k)/dtsoi
            c(j,l) = dqodw2
         end do
     
* node j=nsl
 
         j = nsl
         do l = 1, ncal
            k = ind(l)
            num = -2.*(smp(j-1,k)-smp(j,k)) -
     &            1000.*(dzsoi(j-1,k)+dzsoi(j,k))
            den = 1000.*(dzsoi(j-1,k)/hk(j-1,k) + dzsoi(j,k)/hk(j,k))
            den2 = den*den
            qin = num/den
            dqidw1 = (-2.*den*dsmpdw(j-1,k) +
     &               num*1000.*dzsoi(j-1,k)/hk2(j-1,k)*dhkdw(j-1,k))
     &               / den2
            dqidw2 = ( 2.*den*dsmpdw(j  ,k) +
     &               num*1000.*dzsoi(j  ,k)/hk2(j  ,k)*dhkdw(j  ,k))
     &               / den2

            qout = -hk(j,k)
            dqodw1 = -dhkdw(j,k)

            r(j,l) = qtran(k)*root(j,k) + qin - qout
            a(j,l) = -dqidw1
            b(j,l) = dqodw1 - dqidw2 - 1000.*dzsoi(j,k)/dtsoi
            c(j,l) = 0.
         end do

* nodes j=2 to j=nsl-1
     
         do j = 2, nsl-1 
            do l = 1, ncal
               k = ind(l)

               num = -2.*(smp(j-1,k)-smp(j,k)) -
     &               1000.*(dzsoi(j-1,k)+dzsoi(j,k))
               den = 1000.*(dzsoi(j-1,k)/hk(j-1,k) + dzsoi(j,k)/hk(j,k))
               den2 = den*den
               qin = num/den
               dqidw1 = (-2.*den*dsmpdw(j-1,k) +
     &                  num*1000.*dzsoi(j-1,k)/hk2(j-1,k)*dhkdw(j-1,k))
     &                  / den2
               dqidw2 = ( 2.*den*dsmpdw(j  ,k) +
     &                  num*1000.*dzsoi(j  ,k)/hk2(j  ,k)*dhkdw(j  ,k))
     &                  / den2

               num = -2.*(smp(j,k)-smp(j+1,k)) -
     &               1000.*(dzsoi(j,k)+dzsoi(j+1,k))
               den = 1000.*(dzsoi(j,k)/hk(j,k) + dzsoi(j+1,k)/hk(j+1,k))
               den2 = den*den
               qout = num/den
               dqodw1 = (-2.*den*dsmpdw(j  ,k) +
     &                  num*1000.*dzsoi(j  ,k)/hk2(j  ,k)*dhkdw(j  ,k))
     &                  / den2
               dqodw2 = ( 2.*den*dsmpdw(j+1,k) +
     &                  num*1000.*dzsoi(j+1,k)/hk2(j+1,k)*dhkdw(j+1,k))
     &                  / den**2

               r(j,l) = qtran(k)*root(j,k) + qin - qout
               a(j,l) = -dqidw1
               b(j,l) = dqodw1 - dqidw2 - 1000.*dzsoi(j,k)/dtsoi
               c(j,l) = dqodw2

            end do
         end do

* solve for dwat: a, b, c, r, dwat go from 1 to ncal. h2osoi = 1 to npt

         call tridia (ncal ,nsl ,a ,b ,c ,r ,dwat)

* could now update h2osoi = h2osoi + dwat except for one problem:
* need to make sure h2osoi <= watsat. if not, what to do with
* excess water? add total excess water back to soil profile to
* bring any layers less than watsat up to watsat. any remaining
* excess water is sub-surface runoff.

* update water, constraining h2osoi <= watsat. accumulate excess water

         do j = 1, nsl
            do l = 1, ncal
               k = ind(l)
               newwat(k) = h2osoi(j,k) + dwat(j,l)
               xs(k) = xs(k) + max(newwat(k)-watsat(k),0.) * 
     &                 1000.*dzsoi(j,k)
            end do
            do l = 1, ncal
               k = ind(l)
               h2osoi(j,k) = min(watsat(k),newwat(k))
            end do
         end do

* add excess water back to bring soil layers up to saturation

         do j = 1, nsl
            do l = 1, ncal
               k = ind(l)
               x = min((watsat(k)-h2osoi(j,k))*1000.*dzsoi(j,k),xs(k))
               axs(k) = x
            end do
            do l = 1, ncal
               k = ind(l)
               xs(k) = xs(k) - axs(k)
               h2osoi(j,k) = h2osoi(j,k ) + axs(k)/(1000.*dzsoi(j,k))
            end do
         end do

* sub-surface drainage (accumulate over dtlsm/dtsoi iterations)

         do l = 1, ncal
            k = ind(l)
            qdrai(k) = qdrai(k) + xs(k) +
     &                (hk(nsl,k)+dhkdw(nsl,k)*dwat(nsl,l))*dtsoi
         end do

      end do

* sub-surface drainage over time step = dtlsm

      do l = 1, ncal
         k = ind(l)
         qdrai(k) = qdrai(k)/dtlsm
      end do

* limit h2osoi >= watmin. get water needed to bring h2osoi = watmin
* from lower layer. do for all points so inner loop vectorizes

      watmin = 0.01
      do j = 1, nsl-1 
         do k = 1, npt
            if (h2osoi(j,k) .lt. watmin) then
               xs(k) = (watmin-h2osoi(j,k))*1000.*dzsoi(j,k)
            else
               xs(k) = 0.
            end if
            h2osoi(j  ,k) = h2osoi(j  ,k) + xs(k)/(1000.*dzsoi(j  ,k))
            h2osoi(j+1,k) = h2osoi(j+1,k) - xs(k)/(1000.*dzsoi(j+1,k))
         end do
      end do

      j = nsl
      do k = 1, npt
         if (h2osoi(j,k) .lt. watmin) then
            xs(k) = (watmin-h2osoi(j,k))*1000.*dzsoi(j,k)
         else
            xs(k) = 0.
         end if
         h2osoi(j,k) = h2osoi(j,k) + xs(k)/(1000.*dzsoi(j,k))
         qdrai(k) = qdrai(k) - xs(k)/dtlsm
      end do
         
      return
      end
 
