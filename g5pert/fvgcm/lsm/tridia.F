#include <preproc.h>
      subroutine tridia (npt, n, a, b, c, r, u)
 
#include <vartyp.h>      
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       tridia.F
* purpose:           solve tridiagonal system of equations
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* -----------------------------------------------------------------
   
* ------------------------ notes ----------------------------------
* solve for u given the set of equations f * u = r, where u is a 
* vector of length n, r is a vector of length n, and f is an n x n 
* tridiagonal matrix defined by the vectors a, b, c [each of length n]. 
* a(1) and c(n) are undefined and are not referenced by the subroutine. 
 
*    |b(1) c(1)   0  ...                     |   |u(1  )|   |r(1  )|
*    |a(2) b(2) c(2) ...                     |   |u(2  )|   |r(2  )|
*    |               ...                     | * | ...  | = | ...  |
*    |               ... a(n-1) b(n-1) c(n-1)|   |u(n-1)|   |r(n-1)|
*    |               ...   0    a(n  ) b(n  )|   |u(n  )|   |r(n  )|
     
* -----------------------------------------------------------------

* ------------------------ input/output variables -----------------
      integer n                             !number of soil layers
      integer npt                           !number of points
      real a(msl,npt)                       !input vector
      real b(msl,npt)                       !input vector
      real c(msl,npt)                       !input vector
      real r(msl,npt)                       !input vector
      real u(msl,npt)                       !solution vector
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      real gam(msl,npt)
      real bet(npt)
      integer k,j                           ! loop index
* -----------------------------------------------------------------

      if (npt .eq. 0) return

      do k = 1, npt
         bet(k) = b(1,k)
         u(1,k) = r(1,k) / bet(k)
      end do
      do j = 2, n
         do k = 1, npt
            gam(j,k) = c(j-1,k) / bet(k)
            bet(k) = b(j,k) - a(j,k) * gam(j,k)
            u(j,k) = (r(j,k) - a(j,k)*u(j-1,k)) / bet(k)
         end do
      end do
      do j = n-1, 1, -1
         do k = 1, npt
            u(j,k) = u(j,k) - gam(j+1,k) * u(j+1,k)
         end do
      end do
 
      return
      end
 
