#include <preproc.h>
      subroutine basindrv(beglatspmd, endlatspmd, qover, qdrai)

#include <vartyp.h>
#include <lsmpar.h>
#include <lsmspmd.h>
#include <basinpar.h>

* ------------------------ code history ---------------------------
* source file:       basindrv.F
* purpose:           
* date last revised: April 1997
* author:            Mariana Vertenstein
* standardized:      
* reviewed:         
* -----------------------------------------------------------------

* ------------------- common block variables ----------------------
#include <basin.h>
#include <dynmem.h>
#include <basindyn.h>
#include <lsmtcdyn.h>
#include <lsmhis.h>
* -----------------------------------------------------------------

* ------------------- input variables -----------------------------
      integer beglatspmd    !beg lat on proc, degenerates to 1 for nonspmd
      integer endlatspmd    !end lat on proc, degenerates to atmlat for nonspmd
      real qover(begkptspmd(beglatspmd):endkptspmd(endlatspmd))	!surface runoff (mm h2o /s)
      real qdrai(begkptspmd(beglatspmd):endkptspmd(endlatspmd)) !sub-surface runoff (mm h2o /s) 
* -----------------------------------------------------------------

* ---------------------- local variables --------------------------
      integer i,k,l,m,n        !loop/array indices
      real wt                  !subgrid weight
* -----------------------------------------------------------------

* determine total runoff over land points (sum over subgrid points)

      do l = beglptspmd(beglatspmd),endlptspmd(endlatspmd)
         runoff(l) = 0.
      end do

      do m=1,msub
         do l = beglptspmd(beglatspmd),endlptspmd(endlatspmd)
            k = kvec(l,m)       !lsm subgrid vector index
            wt = wsg2g(l,m)     !lsm subgrid weight
            if (wt .ne. 0.) then 
               runoff(l) = runoff(l) + wt*(qover(k)+qdrai(k))
            end if
         end do
      end do

* determine drainage vector. note that drnfrac[n,l] is the 
* percent drainage at land point 'l' to the basin referenced by 
* drnbasn[n,l]. For each basin, runoff os summed over all land 
* points that drain into that basin. Since runoff is to be in 
* units of kg/sec, multiply runoff by the area of the corresponding
* grid cell since conversion from mm/sec to kg/sec follows from:
* kg/sec = mm/sec * m/1000mm * 1000kg/m^3 * [area of grid cell]
* The factor of 1.e6 is inserted since drnarea is in units of km^2.

      do n=1,nbasmax            ! initialize drainage over all basins
         drnvec(n) = 0.
      end do

      do n=1,ndrnmax            ! loop over top ndrnmax basins
         do l = beglptspmd(beglatspmd),endlptspmd(endlatspmd)
            i = drnbasn(n,l)
            drnvec(i) = drnvec(i) +
     $           drnfrac(n,l)*runoff(l)*drnarea(l)*1.e6
         end do
      end do

* initialize accumulated basin drainage if necessary

      if (bascnt .eq. 0) then
         do n=1,ndrn
            bashist(n) = 0.  
         end do
      endif

* accumulate field 

      do n=1,ndrn
         bashist(n) = bashist(n) + drnvec(n)
      end do

* end of history interval: normalize accumulated values by number 
* of counts. counter is off by one because it is updated later

      if (.not. ehi) return

      do n=1,ndrn
         bashist(n) = bashist(n) / float(bascnt+1)
      end do

      return
      end
 
