#include <preproc.h>
      subroutine histslf (name, npt, field, begkpt, endkpt)

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       histslf.F
* purpose:           accumulate single-level field over history time interval
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* [field] is a "small" vector of length [npt], which corresponds
* to locations [begkpt] to [endkpt] in the "big" vector of length
* [kpt]. save accumulated field values in the array [slfval], which
* is dimensioned for [mslflds] as (kpt,mslflds)
* -----------------------------------------------------------------

* ------------------------ input variables ------------------------
      character*8 name   !field name
      integer npt        !number of points
      real field(npt)    !field values for current time step
      integer begkpt     !beginning location in vector of length kpt
      integer endkpt     !ending    location in vector of length kpt
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <dynmem.h>
#include <lsmhis.h>
#include <lsmhisdyn.h>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer i,k        !do loop index
      integer n          !field index
* -----------------------------------------------------------------

* find field index. return if "name" is not on active list

      n = 0
      do i = 1, nslflds
         if (name .eq. slfnam(i)) n = i
      end do
      if (n .eq. 0) return

* initialize accumulated field value if necessary

      if (slfcnt(n) .eq. 0) then
         if      (slftyp(n) .eq. naver) then   !time average field
            do k = begkpt, endkpt
               slfval(k,n) = 0.  
            end do
         else if (slftyp(n) .eq. ninst) then   !instantaneous field value
            do k = begkpt, endkpt
               slfval(k,n) = 1.e36
            end do
         else if (slftyp(n) .eq. nmaxi) then   !maximum field value
            do k = begkpt, endkpt
               slfval(k,n) = -1.e36
            end do
         else if (slftyp(n) .eq. nmini) then   !minimum field value
            do k = begkpt, endkpt
               slfval(k,n) = 1.e36
            end do
         else
            do k = begkpt, endkpt
               slfval(k,n) = 1.e36
            end do
         end if
      end if

* accumulate field 

      if      (slftyp(n) .eq. naver) then      !time average field
         do k = begkpt, endkpt
            slfval(k,n) = slfval(k,n) + field(k-begkpt+1)
         end do

      else if (slftyp(n) .eq. ninst) then      !instantaneous field value
         do k = begkpt, endkpt
            slfval(k,n) = field(k-begkpt+1)
         end do

      else if (slftyp(n) .eq. nmaxi) then      !maximum field value
         do k = begkpt, endkpt
            slfval(k,n) = max( slfval(k,n), field(k-begkpt+1) )
         end do

      else if (slftyp(n) .eq. nmini) then      !minimum field value
         do k = begkpt, endkpt
            slfval(k,n) = min( slfval(k,n), field(k-begkpt+1) )
         end do

      end if

      if (.not. ehi) return

* end of history interval: normalize accumulated values by number 
* of counts. counter is off by one because it is updated later

      if (slftyp(n) .eq. naver) then  
         do k = begkpt, endkpt
            slfval(k,n) = slfval(k,n) / float(slfcnt(n)+1)
         end do
      end if

      return
      end
 
