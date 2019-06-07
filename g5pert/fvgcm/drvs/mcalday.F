      subroutine mcalday(nymd, nhms, calday)
      implicit none

! input:
      integer nymd
      integer nhms
! Output:
      real calday                    ! Julian day (1 to 366 for non-leap year)
                                     ! Julian day (-1 to -367 for   leap year)
! Local:
      logical leapyr                 ! leap_year?
      logical leap_year
      real tsec
      integer n, nsecf, m, mm
      integer dd, ds
      integer days(12)
      integer ny

      data days /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      nsecf(n)  = n/10000*3600 + mod(n,10000)/100* 60 + mod(n,100)

      ny = nymd / 10000
      mm = mod(nymd, 10000) / 100 
      dd = mod(nymd,   100)

      ds = dd -1

      if( mm .ne. 1) then
      do m=1, mm-1
         if( m.eq.2  .and. leap_year(ny) ) then 
             ds = ds + 29
         else
             ds = ds + days(m)
         endif
      enddo
      endif

      tsec = ds * 86400 + nsecf(nhms)

      calday = tsec / 86400.  + 1.
      if( leap_year(ny) ) calday = -calday

      return
      end
