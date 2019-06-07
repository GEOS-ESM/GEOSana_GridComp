#include <preproc.h>
      subroutine mlf2xy (fldvec, fldxy)

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ notes ----------------------------------
* sub-grid to grid mapping: average a subgrid input vector [fldvec]
* of length [kpt] to a 2-d [lsmlon] x [lsmlat] output array [fldxy]
* setting non-land points to zero
* ---------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <dynmem.h>
#include <lsmtc.h>
#include <lsmtcdyn.h>
#include <ncd.h>
* -----------------------------------------------------------------

* ------------------------ input/output variables ---------------
      real fldvec(kpt,msl)          !input  vector field
      real fldxy(lsmlon,lsmlat,msl) !output xy field
* ---------------------------------------------------------------

* ------------------------ local variables ----------------------
      integer i,j,k,l,m,n
* ---------------------------------------------------------------

* initialize field on land points to zero and to the 'special value'
* at all other grid points

      do n = 1,msl
         do j = 1,lsmlat
            do i = 1,numlon(j)
               if (surf2d(i,j) .gt. 0) then    
                  fldxy(i,j,n) = 0.    !land point
               else
                  fldxy(i,j,n) = spval !ocean point
               endif
            end do
            do i = numlon(j)+1,lsmlon
               fldxy(i,j,n) = spval    !outside of reduced grid domain 
            end do
         end do
      end do

* perform sub-grid to grid mapping

      do n = 1,msl              !subsoil layer
         do m = 1,msub          !subgrid points for each land point
            do l = 1,lpt        !land point index for [lsmlon] x [lsmlat] grid
               i = ixy(l)       !longitude index for land point
               j = jxy(l)       !latitude index for land point
               k = kvec(l,m)    !subgrid vector index: [1] to [kpt]
               fldxy(i,j,n) = fldxy(i,j,n) + fldvec(k,n)*wsg2g(l,m)
            end do
         end do
      end do

      return
      end
 
