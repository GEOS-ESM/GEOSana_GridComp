	subroutine time_shift(idate1, itime1, idate2, itime2,    &
                              odate, otime)
!     +                        op)             !optional


        use precision      
        implicit none



	integer(i4), intent(in)  ::  idate1
	integer(i4), intent(in)  ::  itime1
	integer(i4), intent(in)  ::  idate2
	integer(i4), intent(in)  ::  itime2
	integer(i4), intent(out) ::  odate
	integer(i4), intent(out) ::  otime

!        character*1,  intent(in), optional :: op
!        real,   intent(in), optional :: op

! Local Variables

        logical present
        character*256   line 
        integer(i4)  yy1, mm1, dd1, hh1, mn1, ss1
        integer(i4)  yy2, mm2, dd2, hh2, mn2, ss2


        write (line, '(i8.8)')  idate1
        read  (line, '(i4.4,2i2.2)') yy1, mm1, dd1
        write (line, '(i6.6)')  itime1
        read  (line, '(3i2.2)') hh1, mn1, ss1

        write (line, '(i8.8)')  idate2
        read  (line, '(i4.4,2i2.2)') yy2, mm2, dd2
        write (line, '(i6.6)')  itime2
        read  (line, '(3i2.2)') hh2, mn2, ss2
!        write(*,*) 'hh2 mn2 ss2', hh2, mn2, ss2


!        if ( present(op) ) then
!            write(*,*)'test'
!        else
            call p_shift(yy1, mm1, dd1, hh1, mn1, ss1,        &
                         yy2, mm2, dd2, hh2, mn2, ss2)
!        endif


        write (line, '(i4.4,2i2.2)') yy1, mm1, dd1
        read  (line, '(i8.8)')  odate
        write (line, '(3i2.2)') hh1, mn1, ss1
        read  (line, '(i6.6)')  otime
         


	return

	end


!
!
!



	subroutine p_shift(yy, mm, dd, hh, mn, ss,             &
                           yy2, mm2, dd2, hh2, mn2, ss2)

        use precision
        implicit none

                   
        integer(i4), intent(inout)  ::  yy,  mm,  dd,  hh,  mn,  ss
        integer(i4), intent(in)     ::  yy2, mm2, dd2, hh2, mn2, ss2

! Local Variables


        logical leap_year

        integer(i4) days(12)

        data days /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

        ss=ss+ss2
        mn=mn+mn2
        hh=hh+hh2

        dd=dd+dd2
        mm=mm+mm2
        yy=yy+yy2

        if (ss >=60) then
           ss = ss - 60
           mn = mn +1
        endif

	if (mn >= 60) then
           mn = mn - 60
           hh = hh + 1
        endif

        if (hh >= 24) then
           hh = hh - 24
           dd = dd + 1
        endif

        if (mm==2 .and. leap_year(yy)) then
          days(mm)=29
        endif

        if (dd > days(mm)) then
           dd = dd - days(mm)
           mm = mm + 1
        endif

        if (mm > 12) then
           mm = mm - 12
           yy = yy + 1
        endif

!	write(*,*) 'hh mn ss', hh, mn, ss


	return
	end
