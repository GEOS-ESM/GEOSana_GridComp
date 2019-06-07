#include <preproc.h>
      subroutine rank (n, a, miss, i_1, i_2)
      implicit none

* ------------------------ code history ---------------------------
* source file:       rank.F
* purpose:           return index of the largest number in the array arrin
* date last revised: April 1996
* author:            Gordon Bonan
* standardized:
* reviewed:
* -----------------------------------------------------------------

* ------------------------ input variables ------------------------
      integer n            !array length
      real a(0:n)          !array
      integer miss         !missing data value
* -----------------------------------------------------------------

* ------------------------ output variables -----------------------
      integer i_1          !index to largest value in array
      integer i_2          !index to second largest value in array
* -----------------------------------------------------------------

* ------------------------ local variables ------------------------
      real a_max           !maximum value in array
      integer i            !array index
      real del_max,del_tst !delta for maximum test
* -----------------------------------------------------------------

* find index of largest non-zero number.
* note that the following test includes del_max so only 
* 6 significant figures are retained in the comparison.

      i_1   = miss
      a_max = -9999.
      del_max = a_max/1.e6

      do i = 0, n
         if (a(i).gt.0. .and. (a(i)-del_max).gt.a_max) then
            a_max = a(i)
            i_1  = i
            del_max = a_max/1.e6
         else
            del_tst = a_max/1.e7
            if (a(i).gt. 0. .and. (a(i)-del_tst).gt.a_max) then
               write(6,*)'rank: WARNING maximum test failed',
     $              ' at 10^6 and passed at 10^7'
            endif
         end if
      end do

* i_1 = miss indicates no values > 0. this is an error

      if (i_1 .eq. miss) then
         write (6,*) 'rank error: i_1 = missing'
         stop
      end if

* find index of second largest non-zero number. 
* i_2 = miss if there are no more values > 0

      i_2   = miss
      a_max = -9999.

      do i = 0, n
         if (i .ne. i_1) then
            if (a(i).gt.0. .and. a(i).gt.a_max) then
               a_max = a(i)
               i_2  = i
            end if
         end if
      end do

      return
      end
 
