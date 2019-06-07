#include <preproc.h>
      subroutine stomata (npt  ,ivt    ,tfrz ,mpe    ,apar  ,
     &                    tv   ,ei     ,ea   ,tgcm   ,pgcm  ,
     &                    o2   ,co2    ,igs  ,btran  ,rb    ,
     &                    rs   ,psn    ,qe25 ,aqe    ,kc25  ,
     &                    ko25 ,vcmx25 ,akc  ,ako    ,avcmx ,
     &                    bp   ,mp     ,foln ,folnmx ,c3psn )

#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       stomata.F
* purpose:           leaf stomatal resistance and leaf photosynthesis
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
* input
      integer npt      !number of points
      integer ivt(npt) !vegetation physiology type

      real tfrz        !freezing point (kelvin)
      real mpe         !prevents division by zero errors

      real tv(npt)     !foliage temperature (kelvin)
      real ei(npt)     !vapor pressure inside leaf (sat vapor press at tv) (pa)
      real ea(npt)     !vapor pressure of canopy air (pa)
      real apar(npt)   !par absorbed per unit lai (w/m**2)
      real o2(npt)     !atmospheric o2 concentration (pa)
      real co2(npt)    !atmospheric co2 concentration (pa)
      real pgcm(npt)   !air pressure at agcm reference height (pa)
      real tgcm(npt)   !air temperature at agcm reference height (kelvin)
      real igs(npt)    !growing season index (0=off, 1=on)
      real btran(npt)  !soil water transpiration factor (0 to 1)
      real foln(npt)   !foliage nitrogen concentration (%)
      real rb(npt)     !boundary layer resistance (s/m)

      real qe25(mvt)   !quantum efficiency at 25c (umol co2 / umol photon)
      real ko25(mvt)   !o2 michaelis-menten constant at 25c (pa)
      real kc25(mvt)   !co2 michaelis-menten constant at 25c (pa)
      real vcmx25(mvt) !maximum rate of carboxylation at 25c (umol co2/m**2/s)
      real aqe(mvt)    !q10 for qe25
      real ako(mvt)    !q10 for ko25
      real akc(mvt)    !q10 for kc25
      real avcmx(mvt)  !q10 for vcmx25
      real bp(mvt)     !minimum leaf conductance (umol/m**2/s)
      real mp(mvt)     !slope for conductance-to-photosynthesis relationship 
      real folnmx(mvt) !foliage nitrogen concentration when f(n)=1 (%)
      real c3psn(mvt)  !photosynthetic pathway: 0. = c4, 1. = c3

* output
      real rs(npt)     !leaf stomatal resistance (s/m)
      real psn(npt)    !foliage photosynthesis (umol co2 /m**2/ s) [always +]
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------

      integer i,k,ik   !array/loop index
      integer ncal     !number of points to calculate rs
      integer ipt(npt) !pointer
      integer iter     !iteration index
      integer niter    !number of iterations

      data niter /3/
      save niter

      real ab          !used in statement functions
      real bc          !used in statement functions
      real f1          !generic temperature response (statement function)
      real f2          !generic temperature inhibition (statement function)
      real tc          !foliage temperature (degree celsius)
      real cs          !co2 concentration at leaf surface (pa)
      real kc          !co2 michaelis-menten constant (pa)
      real ko          !o2 michaelis-menten constant (pa)
      real a,b,c,q     !intermediate calculations for rs
      real r1,r2       !roots for rs
      real fnf         !foliage nitrogen adjustment factor (0 to 1)
      real ppf         !absorb photosynthetic photon flux (umol photons/m**2/s)
      real wc          !rubisco limited photosynthesis (umol co2/m**2/s)
      real wj          !light limited photosynthesis (umol co2/m**2/s)
      real we          !export limited photosynthesis (umol co2/m**2/s)
      real cp(npt)     !co2 compensation point (pa)
      real ci(npt)     !internal co2 (pa)
      real awc(npt)    !intermediate calcuation for wc
      real vcmx(npt)   !maximum rate of carboxylation (umol co2/m**2/s)
      real j(npt)      !electron transport (umol co2/m**2/s)
      real cea(npt)    !constrain ea or else model blows up
      real cf(npt)     !s m**2/umol -> s/m

      f1(ab,bc) = ab**((bc-25.)/10.)
      f2(ab) = 1. + exp((-2.2e05+710.*(ab+273.16))/(8.314*(ab+273.16)))
* -----------------------------------------------------------------

* initialize rs=rsmax and psn=0 because will only do calculations
* for apar > 0, in which case rs <= rsmax and psn >= 0

      ncal = 0
      do k = 1, npt
         cf(k) = pgcm(k)/(8.314*tgcm(k))*1.e06 
         rs(k) = 1./bp(ivt(k)) * cf(k)
         psn(k) = 0.
         if (apar(k) .gt. 0.) then
            ncal = ncal + 1
            ipt(ncal) = k
         end if
      end do
      if (ncal.eq.0) then
         return
      end if	
C
      do ik = 1, ncal
         k = ipt(ik)
         i = ivt(k)
         fnf = min( foln(k)/max(mpe,folnmx(i)), 1.0 )
         tc = tv(k)-tfrz                            
         ppf = 4.6*apar(k)                  
         j(k) = ppf*qe25(i)
         kc = kc25(i) * f1(akc(i),tc)       
         ko = ko25(i) * f1(ako(i),tc)
         awc(k) = kc * (1.+o2(k)/ko)
         cp(k) = 0.5*kc/ko*o2(k)*0.21
         vcmx(k) = vcmx25(i) * f1(avcmx(i),tc) / f2(tc) * fnf * btran(k)

* first guess ci

         ci(k) = 0.7*co2(k)*c3psn(i) + 0.4*co2(k)*(1.-c3psn(i))  

* rb: s/m -> s m**2 / umol

         rb(k) = rb(k)/cf(k) 

* constrain ea

         cea(k) = max(0.25*ei(k)*c3psn(i)+0.40*ei(k)*(1.-c3psn(i)), 
     &                  min(ea(k),ei(k)) ) 
      end do

* ci iteration

      do iter = 1, niter
         do ik = 1, ncal
            k = ipt(ik)
            i = ivt(k)
            wj = max(ci(k)-cp(k),0.)*j(k)/(ci(k)+2.*cp(k))*c3psn(i) +
     &           j(k)*(1.-c3psn(i))
            wc = max(ci(k)-cp(k),0.)*vcmx(k)/(ci(k)+awc(k))*c3psn(i) +
     &           vcmx(k)*(1.-c3psn(i))
            we = 0.5*vcmx(k)*c3psn(i) + 
     &           4000.*vcmx(k)*ci(k)/pgcm(k)*(1.-c3psn(i))
            psn(k) = min(wj,wc,we) * igs(k)
         end do
         do ik = 1,ncal
            k = ipt(ik)
            i = ivt(k)
            cs = max( co2(k)-1.37*rb(k)*pgcm(k)*psn(k), mpe )
            a = mp(i)*psn(k)*pgcm(k)*cea(k) / (cs*ei(k)) + bp(i)
            b = ( mp(i)*psn(k)*pgcm(k)/cs + bp(i) ) * rb(k) - 1.
            c = -rb(k)
            if (b .ge. 0.) then
               q = -0.5*( b + sqrt(b*b-4.*a*c) )
            else
               q = -0.5*( b - sqrt(b*b-4.*a*c) )
            end if
            r1 = q/a
            r2 = c/q
            rs(k) = max(r1,r2)
            ci(k) = max( cs-psn(k)*pgcm(k)*1.65*rs(k), 0. )
         end do
      end do

* rs, rb:  s m**2 / umol -> s/m 

      do ik = 1, ncal
         k = ipt(ik)
         rs(k) = rs(k)*cf(k)
         rb(k) = rb(k)*cf(k) 
      end do

      return
      end
 
