#include <preproc.h>
      subroutine laktem (npt   ,nsl   ,errsoi,dtim ,nstep  ,
     &                   zsoi  ,dzsoi ,tsoi  ,sag  ,lakpoi ,
     &                   fira  ,fsh   ,fgev  ,lat  ,fv     , 
     &                   fsm   ,tg    ,z0mg  ,ist  ,istdlak)
 
#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       laktem.F
* purpose:           lake temperatures
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------
 
* ------------------------ notes ----------------------------------
* calculate lake temperatures from one-dimensional thermal
* stratification model based on eddy diffusion concepts to 
* represent vertical mixing of heat
 
* d ts    d            d ts     1 ds
* ---- = -- [(km + ke) ----] + -- --
*  dt    dz             dz     cw dz   
 
* where: ts = temperature (kelvin)
*         t = time (s)
*         z = depth (m)
*        km = molecular diffusion coefficient (m**2/s)
*        ke = eddy diffusion coefficient (m**2/s)
*        cw = heat capacity (j/m**3/kelvin)
*         s = heat source term (w/m**2)

* there are two types of lakes: 
*    deep lakes are 50 m. shallow lakes are 10 m deep.
*    for unfrozen deep lakes:    ke > 0 and    convective mixing
*    for unfrozen shallow lakes: ke = 0 and no convective mixing

* use crank-nicholson method to set up tridiagonal system of equations to
* solve for ts at time n+1, where the temperature equation for layer i is
* r_i = a_i [ts_i-1] n+1 + b_i [ts_i] n+1 + c_i [ts_i+1] n+1

* the solution conserves energy as

* cw*([ts(  1)] n+1 - [ts(  1)] n)*dz(  1)/dt + ... +
* cw*([ts(nsl)] n+1 - [ts(nsl)] n)*dz(nsl)/dt = fin

* where 
* [ts] n   = old temperature (kelvin)
* [ts] n+1 = new temperature (kelvin)
* fin      = heat flux into lake (w/m**2)
*          = beta*sag-fira-fsh-fgev-fsm + phi(1) + ... + phi(nsl) 
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
* input
      integer npt           !number of points
      integer nsl           !number of soil layers
      integer nstep         !time step
      integer ist(npt)      !"soil texture" type
      integer istdlak       !value of ist for deep lakes

      real dzsoi(msl,npt)   !soil layer thickness (m)
      real zsoi(msl,npt)    !depth (m)

      real sag(npt)         !absorbed solar radiation (w/m**2) 
      real fira(npt)        !net infrared radiation (w/m**2) [+ = to atm]
      real fsh(npt)         !sensible heat flux (w/m**2) [+ = to atm]
      real fgev(npt)        !ground latent heat flux (w/m**2) [+ = to atm]
      real fsm(npt)         !snow melt (w/m**2)
      real tg(npt)          !surface temperature (kelvin)
      real fv(npt)          !friction velocity (m/s)
      real z0mg(npt)        !ground roughness length, momentum (m)
      real lat(npt)         !latitude (radians)

      logical lakpoi(npt)   !true if lake point

      real dtim             !time step (s)
* input/output
      real tsoi(msl,npt)    !lake temperature (kelvin)
      real errsoi(npt)      !lake temperature energy conservation err (w/m**2)
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <phycon.h>
#include <lakcon.h>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer i,j,k,ik      !do loop or array index
      integer nlak          !number of lake points 
      integer ilak(npt)     !index to lake points 
      integer nuf           !number of lakes that are not frozen
      integer ind(npt)      !index to lakes that are not frozen
     
      real a(msl,npt)       !"a" vector for tridiagonal matrix
      real b(msl,npt)       !"b" vector for tridiagonal matrix
      real c(msl,npt)       !"c" vector for tridiagonal matrix
      real r(msl,npt)       !"r" vector for tridiagonal solution
      real u(msl,npt)       !solution vector from tridiagonal solution
      real rhow(msl,npt)    !density of water (kg/m**3)
      real phi(msl,npt)     !solar radiation absorbed by layer (w/m**2)
      real kme(msl,npt)     !molecular + eddy diffusion coefficient (m**2/s)

      real ws(npt)          !surface friction velocity (m/s)
      real ks(npt)          !coefficient
      real in(npt)          !relative flux of solar radiation into layer
      real out(npt)         !relative flux of solar radiation out of layer
      real ri(npt)          !richardson number
      real fin(npt)         !heat flux into lake - flux out of lake (w/m**2)
      real ocvts(npt)       !(cwat*(tsoi[n  ])*dzsoi
      real ncvts(npt)       !(cwat*(tsoi[n+1])*dzsoi

      real m1               !intermediate variable for calculating r, a, b, c
      real m2               !intermediate variable for calculating r, a, b, c
      real m3               !intermediate variable for calculating r, a, b, c
      real ke               !eddy diffusion coefficient (m**2/s)
      real km               !molecular diffusion coefficient (m**2/s)
      real f1               !density of water as a function of temper (kg/m**3)
      real x                !argument for f1
      real zin              !depth at top of layer (m)
      real zout             !depth at bottom of layer (m)
      real drhodz           !d [rhow] /dz (kg/m**4)
      real n2               !brunt-vaisala frequency (/s**2)
      real num              !used in calculating ri
      real den              !used in calculating ri
      real tav              !used in aver temp for convectively mixed layers
      real nav              !used in aver temp for convectively mixed layers
      real phidum           !temporary value of phi
      real u2m              !2 m wind speed (m/s)

* statement function 

      f1(x) = 1000. * ( 1.0 - 1.9549e-05 * (abs(x-277.))**1.68 )
* -----------------------------------------------------------------

* index lake points

      nlak = 0
      do k = 1, npt
         if (lakpoi(k)) then
            nlak = nlak + 1
            ilak(nlak) = k
         end if
      end do

* lake density

      do j = 1, nsl
         do ik = 1, nlak
            k = ilak(ik) 
            rhow(j,k) = f1(tsoi(j,k))
         end do
      end do

* eddy diffusion +  molecular diffusion coefficient:
* eddy diffusion coefficient used for unfrozen deep lakes only

      km = tkwat/cwat
      do ik = 1, nlak
         k = ilak(ik) 
         fin(k) = beta(ist(k))*sag(k) - (fira(k)+fsh(k)+fgev(k)+fsm(k))
         u2m = max(1.0,fv(k)/vkc*log(2./z0mg(k)))
         ws(k) = 1.2e-03 * u2m
         ks(k) = 6.6 * sqrt( abs(sin(lat(k))) ) * (u2m**(-1.84))
      end do

      do j = 1, nsl-1
         do ik = 1, nlak
            k = ilak(ik) 
            drhodz = (rhow(j+1,k)-rhow(j,k)) / (zsoi(j+1,k)-zsoi(j,k))
            n2 = -grav / rhow(j,k) * drhodz
            num = 40. * n2 * (vkc*zsoi(j,k))**2
            den = max( (ws(k)**2) * exp(-2.*ks(k)*zsoi(j,k)), 1.e-10 )
            ri(k) = ( -1. + sqrt( max(1.+num/den, 0.) ) ) / 20.
         end do
         do ik = 1, nlak
            k = ilak(ik) 
            if (ist(k).eq.istdlak .and. tg(k) .gt. tfrz) then
               ke = vkc*ws(k)*zsoi(j,k)/p0 * exp(-ks(k)*zsoi(j,k)) / 
     &              (1.+37.*ri(k)*ri(k))
            else
               ke = 0.
            end if
            kme(j,k) = km + ke 
         end do
      end do

CDIR$ IVDEP
      do ik = 1, nlak
         k = ilak(ik) 
         kme(nsl,k) = kme(nsl-1,k)
      end do

* heat source term: unfrozen lakes only

      do j = 1, nsl
         do ik = 1, nlak
            k = ilak(ik) 
            zin  = zsoi(j,k) - 0.5*dzsoi(j,k)
            zout = zsoi(j,k) + 0.5*dzsoi(j,k)
            in(k)  = exp( -eta(ist(k))*(max(  zin-za(ist(k)),0. )) )
            out(k) = exp( -eta(ist(k))*(max( zout-za(ist(k)),0. )) )
         end do
         do ik = 1, nlak
            k = ilak(ik) 
            if (tg(k) .gt. tfrz) then
               phidum = (in(k)-out(k)) * sag(k) * (1.-beta(ist(k)))
            else if (j.eq.1) then
               phidum= sag(k) * (1.-beta(ist(k)))
            else
               phidum = 0.
            end if
            phi(j,k) = phidum
         end do
      end do

* sum cwat*tsoi*dzsoi for energy check

      do ik = 1, nlak
         k = ilak(ik)
         ocvts(k) = 0.
      end do

      do j = 1, nsl
         do ik = 1, nlak
            k = ilak(ik)
            ocvts(k) = ocvts(k) + cwat*tsoi(j,k)*dzsoi(j,k) 
         end do
      end do

* set up vector r and vectors a, b, c that define tridiagonal matrix

      j = 1
      do ik = 1, nlak
         k = ilak(ik)
         m2 = dzsoi(j,k)/kme(j,k) + dzsoi(j+1,k)/kme(j+1,k)
         m3 = dtim/dzsoi(j,k)
         r(j,ik) = tsoi(j,k) + (fin(k)+phi(j,k))*m3/cwat -
     &             (tsoi(j,k)-tsoi(j+1,k))*m3/m2
         a(j,ik) = 0.
         b(j,ik) = 1. + m3/m2
         c(j,ik) = -m3/m2
      end do

      j = nsl
      do ik = 1, nlak
         k = ilak(ik)
         m1 = dzsoi(j-1,k)/kme(j-1,k) + dzsoi(j,k)/kme(j,k)
         m3 = dtim/dzsoi(j,k)
         r(j,ik) = tsoi(j,k) + phi(j,k)*m3/cwat +
     &             (tsoi(j-1,k)-tsoi(j,k))*m3/m1
         a(j,ik) = -m3/m1
         b(j,ik) = 1. + m3/m1
         c(j,ik) = 0.
      end do

      do j = 2, nsl-1
         do ik = 1, nlak
            k = ilak(ik)
            m1 = dzsoi(j-1,k)/kme(j-1,k) + dzsoi(j  ,k)/kme(j  ,k)
            m2 = dzsoi(j  ,k)/kme(j  ,k) + dzsoi(j+1,k)/kme(j+1,k)
            m3 = dtim/dzsoi(j,k)
            r(j,ik) = tsoi(j,k) + phi(j,k)*m3/cwat +
     &                (tsoi(j-1,k)-tsoi(j,k))*m3/m1 -
     &                (tsoi(j,k)-tsoi(j+1,k))*m3/m2
            a(j,ik) = -m3/m1
            b(j,ik) = 1. + m3/m1 + m3/m2
            c(j,ik) = -m3/m2
         end do
      end do

* solve for tsoi: a, b, c, r, u go from 1 to nsoi. tsoi = 1 to npt

      call tridia (nlak ,nsl ,a ,b ,c ,r ,u) 
      do j = 1, nsl
         do ik = 1, nlak
            k = ilak(ik)
            tsoi(j,k) = u(j,ik)
         end do
      end do

* convective mixing: make sure cwat*dzsoi*ts is conserved. mixing
* is only allowed for unfrozen deep lakes. mix every 3 time steps

      if (mod(nstep,3) .eq. 0) then
         nuf = 0
         do ik = 1, nlak
            k = ilak(ik)
            if (ist(k).eq.istdlak .and. tg(k) .gt. tfrz) then
               nuf = nuf + 1
               ind(nuf) = k
            end if
         end do

         do j = 1, nsl-1
            do ik = 1, nuf
               k = ind(ik)
c
c Put an epsilon on density comparison to achieve smooth solution separation
c behavior
c
               if (rhow(j,k) .gt. rhow(j+1,k)+1.e-4) then
                  tav = 0.
                  nav = 0.
                  do i = 1, j+1
                     tav = tav + tsoi(i,k)*dzsoi(i,k)
                     nav = nav + dzsoi(i,k)
                  end do
                  tav = tav/nav

                  do i = 1, j+1
                     tsoi(i,k) = tav
                     rhow(i,k) = f1( tsoi(i,k) )
                  end do

               end if
            end do
         end do
      end if

* sum cwat*tsoi*dzsoi and total energy into lake for energy check

      do ik = 1, nlak
         k = ilak(ik)
         ncvts(k) = 0.
      end do

      do j = 1, nsl
         do ik = 1, nlak
            k = ilak(ik)
            ncvts(k) = ncvts(k) + cwat*tsoi(j,k)*dzsoi(j,k) 
            fin(k) = fin(k) + phi(j,k)
         end do
      end do

      do ik = 1, nlak
         k = ilak(ik)
         errsoi(k) = (ncvts(k)-ocvts(k)) / dtim - fin(k)
      end do

      return
      end
 
