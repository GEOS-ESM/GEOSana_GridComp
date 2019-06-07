
      subroutine linint(im,km,pin,q1,undef,pout,q2)
 
! Linear interpolation of 2d field (length,level) on eta to
! a specified pressure level (length).  Assumes input field 
! is top-down (k=1 is model top).


! Input
      real pin(im,km)       ! coordinate at which q1 is defined
      real q1(im,km)       ! original data
      real undef           ! undefined value assigned if out of bound

! Output
      real q2(im)          ! output interpolated data
      real pout            ! output coordinate 
 
! local variables
      real s
      logical notfound	   ! true = continue to search for bounding model levels

      do i=1,im

! Start searching
      q2(i) = undef
      if(pout .ge. pin(i,1) .and. pout .le. pin(i,km) ) then
        notfound = .true.
        L = 0
        do while (L.lt.(km-1) .and. notfound) 
          L = L + 1
          if (pout .ge. pin(i,L) .and. pout .le. pin(i,L+1)) then
            s  = (pout-pin(i,L)) / (pin(i,L+1) -pin(i,L))
            q2(i) = q1(i,L) + s*(q1(i,L+1) - q1(i,L))
            notfound = .false.
          endif
        enddo
      endif
      enddo

      return
      end



