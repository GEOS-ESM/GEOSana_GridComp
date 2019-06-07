#include <preproc.h>
#ifdef linux
#  define REAL16 double precision
#else
#  define REAL16 real*16
#endif
      subroutine gaulat(a, w, k)

#include <vartyp.h>

* ------------------------ code history ------------------------------
* source file:       gaulat.F
* purpose:           gaussian latitudes: does not give ccm latitudes exactly
*                    so can not use if coupled to ccm
* date last revised: July 1995 - lsm version 1
* author:            Gordon Bonan, derived from ccm2
* --------------------------------------------------------------------

* ------------------------ input/output variables --------------------
* input
      integer k            ! number of latitudes pole to pole

* output
      real a(k)            ! sine of latitudes
      real w(k)            ! gaussian weights
* --------------------------------------------------------------------

* ------------------------ local variables ---------------------------
      pointer (psinlat,sinlat)
      pointer (pwgt,wgt)
      REAL16 sinlat(k)    ! sine of latitudes
      REAL16 wgt(k)       ! gaussian weights
      REAL16 one          ! 1. in REAL16.  Needed by atan
      REAL16 eps          ! convergence criterion
      parameter (one = 1.)
#if defined( CRAY_T3E )
#include <machine.h>
#else
#  if (defined linux)
      parameter (eps = 2.2204460492503131D-016)
#  else
      parameter (eps = 1.D-27)
#  endif
#endif
      REAL16 pi           ! value of pi
      REAL16 c            ! constant combination
      REAL16 fk           ! real k
      REAL16 xz           ! abscissa estimate
      REAL16 pkm1         ! |
      REAL16 pkm2         ! |-polynomials
      REAL16 pkmrk        ! |
      REAL16 pk           ! |
      REAL16 sp           ! current iteration latitude increment
      REAL16 avsp         ! |sp|
      REAL16 fn           ! real n

      integer kk           ! k/2 (number of latitudes in hemisphere)
      integer is           ! latitude index
      integer iter         ! iteration counter
      integer n,l          ! indices
*
*------------------------------Externals--------------------------------
*
      external bsslzr      ! provides zeroes of Bessel function or 
*                          !  estimates thereof
* --------------------------------------------------------------------
      call memget('GAULAT',psinlat,2*k)
      call memget('GAULAT',pwgt,2*k)
#if defined( CRAY_T3E )
      eps = 10.*epsmac  ! Convergence criteria ten times machine epsilon
#endif
      pi = 4.*atan(one)
 
* The value eps, used for convergence tests in the iterations, 
* can be changed.  Newton iteration is used to find the abscissas.
 
      c = (1.-(2./pi)**2)*0.25
      fk = k
      kk = k/2
      call bsslzr(sinlat,kk)

      do is=1,kk
        xz = cos(sinlat(is)/sqrt((fk+0.5)**2+c))
 
* This is the first approximation to xz
 
        iter = 0
   10   pkm2 = 1.
        pkm1 = xz
        iter = iter + 1
        if (iter.gt.10) then
 
* Error exit
 
          write (6,*) 'gaulat: no convergence in 10 iterations'
          call stoprun
        end if
 
* Computation of the legendre polynomial
 
        do n=2,k
          fn = n
          pk = ((2.*fn-1.)*xz*pkm1-(fn-1.)*pkm2)/fn
          pkm2 = pkm1
          pkm1 = pk
        end do
        pkm1 = pkm2
        pkmrk = (fk*(pkm1-xz*pk))/(1.-xz**2)
        sp = pk/pkmrk
        xz = xz - sp
        avsp = abs(sp)
        if (avsp.gt.eps) go to 10
        sinlat(is) = xz
        wgt(is) = (2.*(1.-xz**2))/(fk*pkm1)**2
      end do
      if (k.ne.kk*2) then
 
* For odd k computation of weight at the equator
 
         sinlat(kk+1) = 0.
         pk = 2./fk**2
         do n=2,k,2
           fn = n
           pk = pk*fn**2/(fn-1.)**2
         end do
         wgt(kk+1) = pk
      end if
 
* Complete the sets of abscissas and weights, using the symmetry.
 
      do n=1,kk
        l = k + 1 - n
        a(n) = sinlat(n)
        a(l) = -sinlat(n)

        w(n) = wgt(n)
        w(l) = wgt(n)
      end do

      call memfree(psinlat)
      call memfree(pwgt)

      return
      end
 
