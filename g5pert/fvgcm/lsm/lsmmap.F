#include <preproc.h>
      subroutine lsmmap (wt)

#include <vartyp.h>
#include <lsmpar.h>
#include <lsmspmd.h>

* ------------------------ code history ------------------------------
* source file:       lsmmap.F
* purpose:           build subgrid <-> grid mapping indices and weights
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      J. Truesdale, Feb. 1996
* reviewed:          G. Bonan, Feb. 1996
* --------------------------------------------------------------------

* ------------------------ notes -------------------------------------
* build mapping indices and weights: [lon] x [lat] grid <->
* [lpt] vector of land points <-> [kpt] vector of subgrid points. 
* allow for variable longitudinal resolution: [nlon] <= [lon]

* example: [kpt] vector of subgrid points -> [lsmlon] x [lsmlat] grid
* for 1 <= l <= lpt: 
* t_xy( ixy(l), jxy(l) ) =  t_vec( kvec(l,   1) ) * wsg2g(l,   1) + ... +
*                           t_vec( kvec(l,msub) ) * wsg2g(l,msub) 

* if the grid cell has less than [msub] subgrid points, an arbitrary 
* value of [kvec] is used to index the subgrid vector. the weight is zero

* example: [lsmlon] x [lsmlat] grid -> [kpt] vector of subgrid points
* for 1 <= k <= kpt: t_vec(k) = t_xy( ixy(klnd(k)), jxy(klnd(k)) )
* --------------------------------------------------------------------

* ------------------------ input/output variables --------------------
* input
      real wt(lsmlon,lsmlat,msub)  !subgrid weights
* --------------------------------------------------------------------

* ------------------------ common block variables ----------------------
#include <dynmem.h>
#include <lsmtc.h>
#include <lsmtcdyn.h>
* --------------------------------------------------------------------

* ------------------------ local variables ---------------------------
      integer i,j,k,l,m            !indices
      real sumwv(lsmlon,lsmlat)    !sum of wt (for error checking)
* --------------------------------------------------------------------

* initialize [kvec] to arbitrary value with weight = 0
* for [msub] subgrid types for [lpt] land cells

      do m = 1, msub
         do l = 1, lpt
            kvec(l,m)  = l
            wsg2g(l,m) = 0.
         end do
      end do

* determine subgrid vector index and weight for land points

      k = 0
      l = 0
      do j = 1, lsmlat
         begkptspmd(j) = k+1
         beglptspmd(j) = l+1
         do i = 1, numlon(j)
            if (surf2d(i,j) .gt. 0) then           !land point
               l = l+1                             !land index
               ixy(l) = i                          !longitude index
               jxy(l) = j                          !latitude index
               do m = 1, msub
                  if (wt(i,j,m) .gt. 0.) then      !valid subgrid point
                     k = k+1                       !subgrid index
                     kvec(l,m) = k                 !subgrid index for land point
                     wsg2g(l,m) = wt(i,j,m)        !subgrid weights
                     klnd(k) = l                   !land index for subgrid point
                  end if
               end do
            end if
         end do
         endkptspmd(j)=k
         endlptspmd(j)=l
      end do
* error check: make sure weights sum to one for each land cell

!$omp parallel do default(shared) private(i,j)
      do j = 1, lsmlat
         do i = 1, numlon(j)
            if (surf2d(i,j) .gt. 0) then
               sumwv(i,j) = 0.
            else
               sumwv(i,j) = 1.
            end if
         end do
      end do

      do m = 1, msub
         do l = 1, lpt
            i = ixy(l)
            j = jxy(l)
            sumwv(i,j) = sumwv(i,j) + wsg2g(l,m)
         end do
      end do

      do j = 1, lsmlat
         do i = 1, numlon(j)
            if (abs(sumwv(i,j)-1.) .gt. 1.e-05) then
               write (*,*) 'LSMMAP fatal error: weights do not sum to 1'
               write (*,*) 'lon = ',i,' lat = ',j,' :sum = ',sumwv(i,j)
               write (*,*) 'surf2d(i,j)= ', surf2d(i,j)   ! added by SJL
               call stoprun
            end if
         end do
      end do

      return
      end
 
