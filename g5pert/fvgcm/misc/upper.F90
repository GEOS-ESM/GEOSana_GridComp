      subroutine upper(string)

!***********************************************************************
!
!     upper.f - change lower case letter to upper case letter          *
!                                                                      *
!     George Lai Tue Jun 28 16:37:00 1994                              *
!                                                                      *
!***********************************************************************

      implicit         none

      character*(*) string
      integer i, n
      integer a, z, dist
#if ! (defined LINUX)
      a = ichar('a')
      z = ichar('z')
      n = len(string)
      dist = ichar('A') - a

      do i = 1,n
        if (ichar(string(i:i)) .ge. a .and.       &
            ichar(string(i:i)) .le. z) then
          string(i:i) = char(ichar(string(i:i))+dist)
        endif
      end do
#endif

      return
      end
