      integer function len_trim (string)

!***********************************************************************
!                                                                      * 
!     len_trim.f - return the length of string without the trailing    * 
!                  blanks and tabs                                     * 
!                                                                      * 
!     George Lai Sat May 29 10:39:00 1993                              * 
!                                                                      * 
!***********************************************************************

      implicit         none

      integer          tab, blank
      parameter        (tab   =  9)
      parameter        (blank = 32)

      character*(*)    string
      integer          n, code

      len_trim = 0
      n = len(string)
      if (n .eq. 0) return

      do while (n .gt. 0)
        code = ichar(string(n:n))
        if (code .eq. blank  .or. code .eq. tab) then
          len_trim = n - 1
	  n = n - 1
        else
	  len_trim = n
          n = 0
        end if
      end do

      return 
      end
