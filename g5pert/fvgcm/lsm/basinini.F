#include <preproc.h>
      subroutine basinini 

#include <vartyp.h>
#include <lsmpar.h>
#include <basinpar.h>

* ------------------------ code history ------------------------------
* source file:       basini.F
* purpose:           determine basin matrix over land points
* date last revised: April 1997
* author:            Mariana Vertenstein
* standardized:      
* reviewed:          
* --------------------------------------------------------------------

* ------------------------ common block variables --------------------
#include <basin.h>
#include <dynmem.h>
#include <basindyn.h>
#include <lsmtc.h>
#include <lsmctl.h>
* --------------------------------------------------------------------

* ------------------------ local variables ---------------------------
      integer i,j,l,n    !indices
* --------------------------------------------------------------------

* determine drainage basin matrix over land points

      l = 0                                
      do j = 1, lsmlat
         do i = 1, numlon(j)
            if (surf2d(i,j) .gt. 0) then                 !land point
               l = l+1
               do n=1,ndrnmax
                  drnfrac(n,l) = drngfrac(n,i,j)
                  drnbasn(n,l) = drngbasn(n,i,j)
                  drnarea  (l) = drngarea  (i,j)
               end do
            endif
         end do
      end do

      if (l .ne. lpt) then
         write (6,*) 'BASININI fatal error: discrepancy in land points'
         write (6,*) 'expecting ',lpt,' points but found ',l
         call stoprun
      end if
               
* set flag to write basin history field

      ncbasin = .true.

* initialize history counter

      if (nsrest .eq. 0 .or. nsrest .eq. 3) then
         do i = 1, ndrn
            bascnt = 0
         end do
      endif

      return
      end
 
