
      subroutine trimleft (string)

!***********************************************************************
!                                                                      *
!     trimleft.f - trim blanks and tabs from left of the string        *
!                                                                      *
!     Note:    This is a 'in-place' routine, the input string is       *
!              also the output string.                                 *
!                                                                      *
!     Lats Modified:  Mon Aug 25 13:45:52 EDT 1997                     *
!                                                                      *
!***********************************************************************

      implicit         none

      integer          tab, blank
      parameter        (tab   =  9)
      parameter        (blank = 32)

      character*(*)    string
      logical          found
      integer          n, code, first, i, j

      n = len(string)
      if (n .eq. 0) return

      i = 1
      found = .false.

      do while (i .le. n .and. .not. found)
        code = ichar(string(i:i))
        if (code .eq. blank .or. code .eq. tab) then
          i = i + 1
        else
          found = .true.
          first = i
        end if
      end do

      if (found .and. first .gt. 1) then
        do i = first, n
          j = i - first + 1
          string(j:j) = string(i:i)
          string(i:i) = ' '
        end do
      endif

      return
      end

