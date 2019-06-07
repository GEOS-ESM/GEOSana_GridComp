#include <preproc.h>
      subroutine soitem (npt ,nsl ,dtim ,tsoi   ,dzsoi ,dzsur ,
     &                   fgr ,tk  ,cv   ,lakpoi ,errsoi)
 
#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       soitem.F
* purpose:           soil temperatures
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------
 
* ------------------------ notes ----------------------------------
* calculate soil temperatures from one-dimensional thermal diffusion 
* equation using apparent heat capacity to account for phase change 
 
*     d     d ts       d ts
*     -- (k ----) = ca ----
*     dz     dz         dt
 
* where: k  = thermal conductivity (w/m/kelvin)
*        ca = apparent specific heat capacity (j/m**3/kelvin)
*        ts = temperature (kelvin)
*         z = depth (m)
*         t = time (s)

* use crank-nicholson method to set up tridiagonal system of equations to
* solve for ts at time n+1, where the temperature equation for layer i is
* r_i = a_i [ts_i-1] n+1 + b_i [ts_i] n+1 + c_i [ts_i+1] n+1

* the solution conserves energy as

* cv(  1)*([ts(  1)] n+1 - [ts(  1)] n)*dz(  1)/dt + ... +
* cv(nsl)*([ts(nsl)] n+1 - [ts(nsl)] n)*dz(nsl)/dt = fgr

* where 
* cv       = heat capacity (j/m**3/kelvin)
* dz       = thickness (m)
* dt       = time step (s)
* [ts] n   = old temperature (kelvin)
* [ts] n+1 = new temperature (kelvin)
* fgr      = heat flux into the soil (w/m**2)
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
* input
      integer npt           !number of points
      integer nsl           !number of soil layers

      real dzsoi(msl,npt)   !soil node thickness (m)
      real cv(msl,npt)      !specific heat capacity (j/m**3/kelvin)
      real tk(msl,npt)      !thermal conductivity (w/m/kelvin)

      real dtim             !time step (s)
      real fgr(npt)         !heat flux into soil (w/m**2) [+ = into soil]
      real dzsur(npt)       !thickness first soil layer + snow (m)

      logical lakpoi(npt)   !true if lake point

* input/output
      real tsoi(msl,npt)    !soil temperature (kelvin)
      real errsoi(npt)      !soil temperature energy conservation error (w/m**2)
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer j,k,ik        !do loop or array index
      integer nsoi          !number of soil points 
      integer isoi(npt)     !index to soil points 
     
      real a(msl,npt)       !"a" vector for tridiagonal matrix
      real b(msl,npt)       !"b" vector for tridiagonal matrix
      real c(msl,npt)       !"c" vector for tridiagonal matrix
      real r(msl,npt)       !"r" vector for tridiagonal solution
      real u(msl,npt)       !solution vector from tridiagonal solution
      real dz(msl,npt)      !dzsoi adjusted for snow
      real ocvts(npt)       !sum (cv*(tsoi[n  ])*dzsoi
      real ncvts(npt)       !sum (cv*(tsoi[n+1])*dzsoi
      real m1               !intermediate variable for calculating r, a, b, c
      real m2               !intermediate variable for calculating r, a, b, c
      real m3               !intermediate variable for calculating r, a, b, c
* -----------------------------------------------------------------

* index soil points

      nsoi = 0
      do k = 1,npt
         if (.not. lakpoi(k)) then
            nsoi = nsoi + 1
            isoi(nsoi) = k
         end if
      end do

* need to use blended snow/soil for first layer

      do ik = 1, nsoi
         k = isoi(ik)
         dz(1,k) = dzsur(k)
      end do
      do j = 2, nsl
         do ik = 1, nsoi
            k = isoi(ik)
            dz(j,k) = dzsoi(j,k)
         end do
      end do

* sum cv*tsoi*dz for energy check

      do ik = 1, nsoi
         k = isoi(ik)
         ocvts(k) = 0.
      end do

      do j = 1, nsl
         do ik = 1, nsoi
            k = isoi(ik)
            ocvts(k) = ocvts(k) + cv(j,k)*tsoi(j,k)*dz(j,k) 
         end do
      end do

* set up vector r and vectors a, b, c that define tridiagonal matrix

      j = 1
      do ik = 1, nsoi
         k = isoi(ik)
         m2 = dz(j,k)/tk(j,k) + dz(j+1,k)/tk(j+1,k)
         m3 = dtim / (dz(j,k)*cv(j,k))
         r(j,ik) = tsoi(j,k) + fgr(k)*m3 - (tsoi(j,k)-tsoi(j+1,k))*m3/m2
         a(j,ik) = 0.
         b(j,ik) = 1. + m3/m2
         c(j,ik) = -m3/m2
      end do

      j = nsl
      do ik = 1, nsoi
         k = isoi(ik)
         m1 = dz(j-1,k)/tk(j-1,k) + dz(j,k)/tk(j,k)
         m3 = dtim / (dz(j,k)*cv(j,k))
         r(j,ik) = tsoi(j,k) + (tsoi(j-1,k)-tsoi(j,k))*m3/m1
         a(j,ik) = -m3/m1
         b(j,ik) = 1. + m3/m1
         c(j,ik) = 0.
      end do

      do j = 2, nsl-1
         do ik = 1, nsoi
            k = isoi(ik)
            m1 = dz(j-1,k)/tk(j-1,k) + dz(j  ,k)/tk(j  ,k)
            m2 = dz(j  ,k)/tk(j  ,k) + dz(j+1,k)/tk(j+1,k)
            m3 = dtim / (dz(j,k)*cv(j,k))
            r(j,ik) = tsoi(j,k) + (tsoi(j-1,k)-tsoi(j,k))*m3/m1 -
     &                (tsoi(j,k)-tsoi(j+1,k))*m3/m2
            a(j,ik) = -m3/m1
            b(j,ik) = 1. + m3/m1 + m3/m2
            c(j,ik) = -m3/m2
         end do
      end do

* solve for tsoi: a, b, c, r, u go from 1 to nsoi. tsoi = 1 to npt

      call tridia (nsoi ,nsl ,a ,b ,c ,r ,u) 
      do j = 1, nsl
         do ik = 1, nsoi
            k = isoi(ik)
            tsoi(j,k) = u(j,ik)
         end do
      end do

* sum cv*tsoi*dz for energy check and check for energy conservation

      do ik = 1, nsoi
         k = isoi(ik)
         ncvts(k) = 0.
      end do

      do j = 1,nsl
         do ik = 1, nsoi
            k = isoi(ik)
            ncvts(k) = ncvts(k) + cv(j,k)*tsoi(j,k)*dz(j,k) 
         end do
      end do

      do ik = 1, nsoi
         k = isoi(ik)
         errsoi(k) = (ncvts(k)-ocvts(k)) / dtim - fgr(k)
      end do

      return
      end
 
