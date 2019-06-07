#include <preproc.h>
      subroutine snoalb (npt,nband,h2osno,ta,coszen,ind,soot,alb)

#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       snoalb.F
* purpose:           snow albedos
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
* input
      integer npt         !number of points 
      integer ind         !0=direct beam, 1=diffuse radiation
      integer nband       !number of waveband classes

      real ta(npt)        !air temperature (kelvin)
      real coszen(npt)    !cosine solar zenith angle
      real soot(npt)      !soot concentration snow (mass fraction)
      real h2osno(npt)    !water equivalent snow (mm h2o)

* output
      real alb(mband,npt) !snow albedo
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer i,j         !loop index
      integer ib          !waveband class
      integer ncal        !number of points with coszen > 0 and h2osno > 0
      integer ipt(npt)    !index for points with h2osno > 0 and coszen > 0

      real b(mband)       !constants for zenith angle correction to grain size
      real e(mband)       !constants for ln(s2)
      real f(mband)       !constants for ln(s2)
      real c1(mband)      !constants for ah
      real c2(mband)      !constants for ah
      real c3(mband)      !constants for ah
      real al(mband)      !minimum albedo (dirty snow)
      real kh1(mband)     !constants for a1
      real kh2(mband)     !constants for a1
      real kl(mband)      !constants for a2

      real udif     !cosine zenith angle for diffuse radiation
      real u(npt)   !cosine zenith angle
      real r(npt)   !snow grain radius (um)
      real rc       !grain radius for cold snow (um)
      real rm       !grain radius for melting snow (um)
      real er(npt)  !effective snow grain radius, with zenith angle effect (um)
      real lsr      !ln(sqrt(er))
      real ls2      !ln(s2), where s2 = empirical soot concentration parameter
      real s(npt)   !soot concentration snow (mass fraction)
      real ls2os    !ln(s2/s)
      real ah       !maximum albedo (clean snow)
      real am       !1/2 (ah+al)
      real a1       !albedo if s < s2
      real a2       !albedo if s > s2
      real sa       !snow albedo
      real sgr      !temporary snow grain radius

      save b,e,f,c1,c2,c3,al,kh1,kh2,kl,rc,rm,udif
      data udif /0.65/
      data rc,rm /80.,1000./
      data b /1.2812,1.2642/
      data e /-6.3550,-5.9328/
      data f /-1.9752,-1.8102/
      data c1 /0.9988,1.0667/
      data c2 /-0.0020,0./
      data c3 /0.,-0.1811/
      data al /0.070,0.059/
      data kh1 /0.450,0.580/
      data kh2 /-0.040,-0.040/
      data kl /0.454,0.567/
* -----------------------------------------------------------------

* zero albedos for all points

      do ib = 1, nband
         do j = 1, npt
            alb(ib,j) = 0.
         end do
      end do

* index points with snow and with coszen > 0

      ncal = 0
      do j = 1, npt
         if (coszen(j).gt.0. .and. h2osno(j).gt.0.) then
            ncal = ncal + 1
            ipt(ncal) = j
         end if
      end do
      if (ncal .eq. 0) then
         return
      endif

      do j = 1, ncal
         i = ipt(j)

* minimum soot

         s(i) = max(soot(i), 1.e-09)

* snow grain radius

         sgr = exp( 4.38+0.23*(ta(i)-263.15) )
         if (ta(i) .lt. 263.15) sgr = rc 
         if (ta(i) .gt. 274.15) sgr = rm       
         r(i) = sgr

* adjust for zenith angle: direct (ind=0) = coszen; diffuse(ind=1) = udif

         u(i) = (udif-coszen(i))*ind + coszen(i)
      end do

* albedo

      do ib = 1, nband
         do j = 1, ncal
            i = ipt(j)
            er(i) = r(i) * (1.+(u(i)-udif)/b(ib))**2
            lsr = log( sqrt(er(i)) )
            ls2 = e(ib) + f(ib)*lsr
            ls2os = ls2 - log(s(i))
            ah = c1(ib) + c2(ib)*sqrt(er(i)) + c3(ib)*lsr
            am = 0.5 * (ah+al(ib))
            a1 = min( ah, 
     &        ah-(ah-am)*exp(-kh1(ib)*ls2os+kh2(ib)*ls2os*ls2os) )
            a2 = max( al(ib), 
     &        al(ib)+(am-al(ib))*exp(-kl(ib)*(log(s(i))-ls2)) )
            if (s(i) .lt. exp(ls2)) then
               sa = a1
            else
               sa = a2
            end if
            alb(ib,i) = sa
         end do
      end do

      return
      end
 
