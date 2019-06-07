#include <preproc.h>
      subroutine rankn (ndim, ai, nmax, indx)
      implicit none

* ------------------------ code history ---------------------------
* source file:       rank.F
* purpose:           return indices of the nmax largest numbers in the 
*                    array arrin
* date last revised: March 1997
* author:            Mariana Vertenstein
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ input variables ------------------------
      integer ndim              !array length
      real ai(ndim)             !input array
      integer nmax              !dimension of output index array 
* -----------------------------------------------------------------

* ------------------------ output variables -----------------------
      integer indx(nmax)        !indices into sorted array
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer indsort(ndim)     !temporary index array
      real arrsort(ndim)        !temporary array
      real item                 !temporary sort variable
      integer i,ii,n,istop      !indices
* -----------------------------------------------------------------

* sort input array in desending order

      do i=1,ndim
         indsort(i) = i
         arrsort(i) = ai(i)
      end do

      do ii=1,ndim-1
         istop=ndim-ii
         do i=1,istop
            if (arrsort(i+1).gt.arrsort(i)) then
               item         = arrsort(i)
               arrsort(i)   = arrsort(i+1)
               arrsort(i+1) = item
               item         = indsort(i)
               indsort(i)   = indsort(i+1)
               indsort(i+1) = item
            endif
         end do
      end do

* keep index of only top nmax elements of sorted array

      do i=1,nmax
         indx(i) = indsort(i)
      end do

      return
      end

 
