      function leap_year(ny)
!
! Determine if year ny is a leap year
!
! Author: S.-J. Lin
      implicit none
      logical leap_year
      integer ny
      integer ny00

!
! No leap years prior to 1900
!
      parameter ( ny00 = 1900 )   ! The threshold for starting leap-year 

      if( mod(ny,4) .eq. 0 .and. ny .ge. ny00 ) then
          leap_year = .true.
      else
          leap_year = .false.
      endif

      return 
      end
