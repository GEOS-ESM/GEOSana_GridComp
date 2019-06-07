#include <preproc.h>
      subroutine histmlf (name, npt, field, begkpt, endkpt)

#include <vartyp.h>
#include <lsmpar.h>

* ------------------------ code history ---------------------------
* source file:       histmlf.F
* purpose:           accumulate multi-level field over history time interval
* date last revised: March 1996 - lsm version 1
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ notes ----------------------------------
* [field] is a "small" multi-level vector of length [npt], corresponding
* to locations [begkpt] to [endkpt] in the "big" vector of length 
* [kpt]. save accumulated field values in the array [mlfval], which
* is dimensioned for [mmlflds] and [msl] levels as (kpt,msl,mmlflds)
* -----------------------------------------------------------------

* ------------------------ input variables ------------------------
      character*8 name    !field name
      integer npt         !number of points
      real field(msl,npt) !field values for current time step
      integer begkpt      !beginning location in vector of length kpt
      integer endkpt      !ending    location in vector of length kpt
* -----------------------------------------------------------------

* ------------------------ common block variables -----------------
#include <dynmem.h>
#include <lsmhis.h>
#include <lsmhisdyn.h>
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      integer i,j,k      !do loop index
      integer n          !field index
* -----------------------------------------------------------------

* find field index. return if "name" is not on active list

      n = 0
      do i = 1, nmlflds
         if (name .eq. mlfnam(i)) n = i
      end do
      if (n .eq. 0) return

* initialize accumulated field value if necessary

      do j = 1, msl
      if (mlfcnt(n) .eq. 0) then
         if      (mlftyp(n) .eq. naver) then   !time average field
            do k = begkpt, endkpt
               mlfval(k,j,n) = 0.  
            end do
         else if (mlftyp(n) .eq. ninst) then   !instantaneous field value
            do k = begkpt, endkpt
               mlfval(k,j,n) = 1.e36
            end do
         else if (mlftyp(n) .eq. nmaxi) then   !maximum field value
            do k = begkpt, endkpt
               mlfval(k,j,n) = -1.e36
            end do
         else if (mlftyp(n) .eq. nmini) then   !minimum field value
            do k = begkpt, endkpt
               mlfval(k,j,n) = 1.e36
            end do
         else
            do k = begkpt, endkpt
               mlfval(k,j,n) = 1.e36
            end do
         end if
      end if
      end do

* accumulate field 

      do j = 1, msl
      if      (mlftyp(n) .eq. naver) then      !time average field
         do k = begkpt, endkpt
            mlfval(k,j,n) = mlfval(k,j,n) + field(j,k-begkpt+1)
         end do

      else if (mlftyp(n) .eq. ninst) then      !instantaneous field value
         do k = begkpt, endkpt
            mlfval(k,j,n) = field(j,k-begkpt+1)
         end do

      else if (mlftyp(n) .eq. nmaxi) then      !maximum field value
         do k = begkpt, endkpt
            mlfval(k,j,n) = max( mlfval(k,j,n), field(j,k-begkpt+1) )
         end do

      else if (mlftyp(n) .eq. nmini) then      !minimum field value
         do k = begkpt, endkpt
            mlfval(k,j,n) = min( mlfval(k,j,n), field(j,k-begkpt+1) )
         end do

      end if
      end do

      if (.not. ehi) return

* end of history interval: normalize accumulated values by number 
* of counts. counter is off by one because it is updated later

      if (mlftyp(n) .eq. naver) then   
         do j = 1, msl
            do k = begkpt, endkpt
               mlfval(k,j,n) = mlfval(k,j,n) / float(mlfcnt(n)+1)
            end do
         end do
      end if

      return
      end
 
