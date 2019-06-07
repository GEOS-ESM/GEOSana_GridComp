#include <preproc.h>
#ifdef COUP_CCM
      subroutine atmini (londim ,latdim ,beglatspmd,endlatspmd,lwupxy ,
     &	                 tsxy   ,asdirxy,asdifxy   ,aldirxy   ,aldifxy,
     &                   snowxy, oroatm ) 
#else
      subroutine atmini (londim ,latdim ,beglatspmd,endlatspmd,lwupxy ,
     &	                 tsxy   ,asdirxy,asdifxy   ,aldirxy   ,aldifxy,
     &                   snowxy) 
#endif

#include <vartyp.h>
#include <lsmpar.h>
#include <lsmspmd.h>
* ------------------------ code history --------------------------------
* source file:       lsmini.F
* purpose:           lsm initialization. return land surface variables to
*                    atm model if initial run
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:      
* reviewed:         
* ----------------------------------------------------------------------

* ------------------------ notes ---------------------------------------
* Initialize land surface variables for atmosphere model.
* Only invoked upon an initial run.
* ----------------------------------------------------------------------

* ------------------------ input/output variables ----------------------
* input from atmospheric model

      integer londim                !number of atm longitude points
      integer latdim                !number of atm latitude points
      integer beglatspmd    !beg lat on proc, degenerates to 1 for nonspmd
      integer endlatspmd    !end lat on proc, degenerates to atmlat for nonspmd

* output to atmospheric model
 
      real asdirxy(londim,latdim)   !albedo - visible waveband, direct
      real asdifxy(londim,latdim)   !albedo - visible waveband, diffuse
      real aldirxy(londim,latdim)   !albedo - near infrared waveband, direct
      real aldifxy(londim,latdim)   !albedo - near infrared waveband, diffuse
      real tsxy   (londim,latdim)   !surface (skin) temperature
      real lwupxy (londim,latdim)   !emitted longwave radiation (w/m**2)
      real snowxy (londim,latdim)   !water equivalent snow (m)
#ifdef COUP_CCM
      real oroatm (londim,latdim)   !0 = non-land. 1 = land
#endif
* ----------------------------------------------------------------------

* ------------------------ common block variables ----------------------
#include <dynmem.h>
#include <lsmtc.h>
#include <lsmtcdyn.h>
#include <lsmtvdyn.h>
#include <phycon.h>
* ----------------------------------------------------------------------

* ------------------------ local variables -----------------------------
      integer i,j,k,l,m               !indices
* ----------------------------------------------------------------------

* ----------------------------------------------------------------------
* return surface variables for atmospheric model
* ----------------------------------------------------------------------

* initialize to zero for land only, other than snowxy which is 
* initialized to zero everywhere
* beglatspmd degenerates to 1 for nonspmd case
* endlatspmd degenerates to lsmlat for nonspmd case
*
      do j = beglatspmd,endlatspmd
         do i = 1, numlon(j)
            if (surf2d(i,j) .gt. 0) then 
               asdirxy(i,j) = 0.
               asdifxy(i,j) = 0.
               aldirxy(i,j) = 0.
               aldifxy(i,j) = 0.
               tsxy(i,j)    = 0.
               lwupxy(i,j)  = 0.
            endif
            snowxy(i,j)  = 0.
         end do
      end do

* [kpt] vector of subgrid points -> [lpt] vector of land points -> 
* [lsmlon] x [lsmlat] grid 

      do m = 1, msub           !subgrid points for each land point
         do l = beglptspmd(beglatspmd),endlptspmd(endlatspmd)
            i = ixy(l)         !longitude index
            j = jxy(l)         !latitude index
            k = kvec(l,m)      !lsm subgrid vector index
            asdirxy(i,j) = asdirxy(i,j) + albd(1,k) * wsg2g(l,m)
            asdifxy(i,j) = asdifxy(i,j) + albi(1,k) * wsg2g(l,m)
            aldirxy(i,j) = aldirxy(i,j) + albd(2,k) * wsg2g(l,m)
            aldifxy(i,j) = aldifxy(i,j) + albi(2,k) * wsg2g(l,m)
            tsxy(i,j)    = tsxy(i,j)    + tg(k)     * wsg2g(l,m)
            lwupxy(i,j)  = lwupxy(i,j)  + sb*(tg(k)**4)  * wsg2g(l,m)
            snowxy(i,j)  = snowxy(i,j)  + h2osno(k)/1000.* wsg2g(l,m)
         end do
      end do

#ifdef COUP_CCM
* return 2-d orography based on surf2d to atm model. 
* NOTE: oroatm received from the atm may be modified 
* (for example, if antartica is set to true in subroutine fsurdat).
* NOTE: if a restart or branch run, the atm obtains oroatm from its 
* own restart dataset.

      do j = beglatspmd,endlatspmd
         do i = 1, numlon(j)
            if (surf2d(i,j) .gt. 0) then
               oroatm(i,j) = 1.
            else
               oroatm(i,j) = 0.
            end if
         end do
      end do
#endif

      return
      end
 
