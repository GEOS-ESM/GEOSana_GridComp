  program read_pblh_lidar_cats_cnn

! !USES:

      use kinds, only: r_kind,r_double,i_kind
      use constants, only: zero,one_tenth,one,deg2rad,rad2deg,three
      use gsi_4dvar, only: l4dvar,l4densvar,time_4dvar,winlen
      use obsmod, only: iadate,offtime_data,bmiss
      use deter_sfc_mod, only: deter_sfc2
      use mpimod, only: npe
      use netcdf

      implicit none

!     Declare local parameters
      integer(i_kind),parameter:: ntime = 236
      integer(i_kind),parameter:: maxnum = 400000 ! 150000
      real(r_double),parameter:: r360 = 360.0_r_double
      real(r_double),parameter:: missing = -99999.9_r_double

      character(len=225):: filename(32767),filename2(ntime),obstype
      character(len=225):: filename_each
      character(20):: sis

      character(10) date
      character*255 argv
      logical first,outside,inflate_error,lexist,more_data
      integer(i_kind) nfile, iarg, argc, iargc
      integer(i_kind) ana_yy1,ana_mm1,ana_dd1,ana_hh1,ana_yy0,ana_mm0,ana_dd0,ana_hh0
      integer(i_kind) iret,istat,i,j,k,kk,ii
      integer(i_kind) ikx,nkx,kx,nreal,ilat,ilon,iout
      character(19) site
      character(8) instrument
      character(10) ymd
      character(4) yyyy
      character(2) mm
      character(2) dd
      character(1) nightday
      integer(i_kind) idx
      integer(i_kind) iyyyy,imm,idd
      integer(i_kind) mm_new, dd_new
      integer(i_kind) unit_number, nr
      integer(i_kind) iostat

      integer(i_kind) :: ncid,ncid2(ntime),ierr,ierr2,dimid1,dimid2,dimid3,dimid4,dimid5,nobs
      integer(i_kind) :: ntime_1440
      integer(i_kind) :: varid1,varid2,varid3,varid4,varid5,varid6
      integer(i_kind) :: varid7,varid8,varid9,varid10,varid11
!     integer(i_kind) :: iyear, imonth, idate, ihour, iminute
      real(r_double), allocatable,dimension(:) :: time
      real(kind=4), allocatable,dimension(:)   :: lat, lon, PBL_z
      real(kind=8), allocatable,dimension(:) :: jdy
      integer(i_kind), allocatable,dimension(:) :: ryear, rmonth, rdate, rhour, rmin
      real(r_double), allocatable,dimension(:) :: rsec
      !real(r_double), allocatable,dimension(:,:) :: ATB

      integer(i_kind) :: total_nobs(ntime)
      integer(i_kind) ana_time(ntime), win_time(ntime), obs_time
      real(r_kind), dimension(maxnum,ntime) :: slat, slon, spblh
      integer(i_kind), dimension(maxnum,ntime) :: sflag_nightday
      integer(i_kind), dimension(maxnum,ntime) :: syear, smonth, sdate, shour, smin
      real(r_double), dimension(maxnum,ntime) :: ssec

      obstype = 'pblh'
      sis = 'catscnn'

!     set analysis window and output file name
      ana_yy0 = 2015
      ana_yy1 = ana_yy0
      ana_mm0 = 8
      ana_dd0 = 3
      ana_hh0 = 9
      win_time(1) = ana_yy0*1000000+ana_mm0*10000+ana_dd0*100+ana_hh0
      print*, win_time(1)
      call add_interval(ana_mm1, ana_dd1, ana_hh1, 3, ana_mm0, ana_dd0, ana_hh0)
      ana_time(1) = ana_yy1*1000000+ana_mm1*10000+ana_dd1*100+ana_hh1
      do i = 2, ntime
         call add_interval(ana_mm1, ana_dd1, ana_hh1, 6, ana_mm0, ana_dd0, ana_hh0)
         ana_mm0 = ana_mm1
         ana_dd0 = ana_dd1
         ana_hh0 = ana_hh1
         win_time(i) = ana_yy0*1000000+ana_mm0*10000+ana_dd0*100+ana_hh0
         call add_interval(ana_mm1, ana_dd1, ana_hh1, 3, ana_mm0, ana_dd0, ana_hh0)
         ana_time(i) = ana_yy1*1000000+ana_mm1*10000+ana_dd1*100+ana_hh1
         print*, i-1, win_time(i-1), ana_time(i-1), win_time(i)
      end do

!     Read in original input data file names
      argc =  iargc()
      if ( argc .lt. 1 ) stop
      iarg = 0
      nfile = 0
      do i = 1, 32767
         iarg = iarg + 1
         if ( iarg .gt. argc ) exit
         call GetArg ( iArg, argv )
         nfile = nfile + 1
         filename(nfile) = argv
         print*, 'filename = ', trim(filename(nfile))
      end do
      print*, 'nfile = ', nfile ! 775

!     Read data from each file and Sort data into analysis windows
      total_nobs = 0 
      unit_number = 10
      do k = 1, nfile

         print*, "k = ", k
!        Check if pblh file exists
         print*, "input filename = ", filename(k)
         inquire(file=trim(filename(k)),exist=lexist)
         if (.not.lexist) cycle

!        Read data

         ! Open the file
         open(unit=unit_number, file=filename(k), status="old", form="unformatted", access="stream")

         ! Read the first 4 bytes (nr: data length)
         read(unit_number) nr
         print *, "Number of elements (nr):", nr

         ! Allocate memory
         allocate(PBL_z(nr), lat(nr), lon(nr))
         allocate(jdy(nr))

         ! Read CNN_PBL_z (float32)
         read(unit_number) PBL_z
         print *, "CNN_PBL_z read successfully."
         print *, "CNN_PBL_z(1)=",PBL_z(1)

         ! Read CNN_lat (float32)
         read(unit_number) lat
         print *, "CNN_lat read successfully."
         print *, "CNN_lat(1)=",lat(1)

         ! Read CNN_lon (float32)
         read(unit_number) lon
         print *, "CNN_lon read successfully."
         print *, "CNN_lon(1)=",lon(1)

         ! Read CNN_jdy (float64)
         read(unit_number) jdy
         print *, "CNN_jdy read successfully."
         print *, "CNN_jdy(1)=",jdy(1)

         ! Deallocate memory and close the file
         !deallocate(PBL_z, lat, lon, jdy)
         !close(unit_number)

         ! Read date
         filename_each = filename(k)
         ymd = filename_each(27:36)
         yyyy = filename_each(27:30)
         mm = filename_each(32:33)
         dd = filename_each(35:36)
         nightday = filename_each(14:15) ! day or night from file name: D or N

         read(yyyy, '(I4)', iostat=iostat) iyyyy
         read(mm, '(I2)', iostat=iostat) imm
         read(dd, '(I2)', iostat=iostat) idd
         !print*, 'yyyy-mm-dd=', ymd
         !print*, 'yyyy,mm,dd=', yyyy,mm,dd
         !print*, 'iyyyy,imm,idd=', iyyyy,imm,idd
         !print*, 'day or night=', nightday

         allocate(ryear(nr),rmonth(nr),rdate(nr),rhour(nr),rmin(nr),rsec(nr))
!        Sorting based on analysis window 
         do i = 1, nr

            ! Obs time info
            !   a) yyyy-mm-dd from file name 
            if (i==1) then
               ryear(i)  = iyyyy
               rmonth(i) = imm
               rdate(i)  = idd
            else
               ryear(i)  = ryear(i-1)
               rmonth(i) = rmonth(i-1)
               rdate(i)  = rdate(i-1)
            end if

            !   b) hh:min:sec from variable in .dat file
            !   Convert Julian date (fractional day-of-year) -> hour,min,sec
            call djul_day_inv(jdy(i),rhour(i),rmin(i),rsec(i))
            !   if yyyy-mm-dd is changed in one file, we need to add one day (yyyy-mm-dd)
            if ( i>1 .and. jdy(i) < jdy(i-1) ) then 
               call add_oneday(mm_new, dd_new, imm, idd)
               rmonth(i) = mm_new
               rdate(i)  = dd_new
            end if
            !print*, 'ryear(i),rmonth(i),rdate(i)=',ryear(i),rmonth(i),rdate(i)

            obs_time = ryear(i)*1000000+rmonth(i)*10000+rdate(i)*100+rhour(i)
            !print*, i, "obs_time=", obs_time
            kk = 0
            do j = 2, ntime
               if (obs_time>=win_time(j-1) .and. obs_time<win_time(j)) then
                  kk = j-1
                  !print*, "obs_time=", obs_time, " kk=", kk
                  !print*, "win_time=", win_time(j-1), win_time(j)
                  !print*, "ana_time=", ana_time(j-1)
                  exit
               end if
            end do

            if (kk == 0) cycle
            total_nobs(kk) = total_nobs(kk) + 1
            ii = total_nobs(kk)
            slat(ii,kk) = lat(i) 
            slon(ii,kk) = lon(i) 
            syear(ii,kk) = ryear(i) 
            smonth(ii,kk) = rmonth(i) 
            sdate(ii,kk) = rdate(i) 
            shour(ii,kk) = rhour(i) 
            smin(ii,kk) = rmin(i) 
            ssec(ii,kk) = rsec(i) 
            spblh(ii,kk) = PBL_z(i)
            if (nightday=='N') then
               sflag_nightday(ii,kk) = 0 ! night
            else
               sflag_nightday(ii,kk) = 1 ! day
            end if

         end do

         deallocate(lat, lon, PBL_z, jdy)
         deallocate(ryear, rmonth, rdate, rhour, rmin, rsec)

      end do ! end of nfile

!     Write into netcdf files for each analysis window
      do i = 1, ntime
         if (total_nobs(i) < 1) cycle

         write(filename2(i), '(A14,I10,A4)') 'cats_pblh_cnn.', ana_time(i), '.nc4'
!        ierr2 =  NF90_CREATE(trim(filename2(i)), NF90_CLOBBER, ncid2(i)) ! overwrite this file if it already exists
         ierr2 =  NF90_CREATE(trim(filename2(i)), NF90_NETCDF4, ncid2(i)) 
         if (ierr2 /= nf90_noerr) call handle_err(ierr2,"create")
!        ierr2 = nf90_open(trim(filename2(i)),nf90_write, ncid2(i))
!        if (ierr2 /= nf90_noerr) call handle_err(ierr2,"open")
         print*, "filename2 = ", filename2(i)
         print*, "total_nobs =", total_nobs(i)

         ierr2 = nf90_def_dim(ncid2(i), 'Time', 1, dimid3)
         ierr2 = nf90_def_dim(ncid2(i), 'nobs', total_nobs(i), dimid1)

         ierr2 = nf90_def_var(ncid2(i), 'Ana_Time', NF90_INT64, dimid3, varid9)
         ierr2 = nf90_def_var(ncid2(i), 'lat', NF90_FLOAT, dimid1, varid1)
         ierr2 = nf90_def_var(ncid2(i), 'lon', NF90_FLOAT, dimid1, varid2)
         ierr2 = nf90_def_var(ncid2(i), 'Year', NF90_INT64, dimid1, varid3)
         ierr2 = nf90_def_var(ncid2(i), 'Month', NF90_INT64, dimid1, varid4)
         ierr2 = nf90_def_var(ncid2(i), 'Date', NF90_INT64, dimid1, varid5)
         ierr2 = nf90_def_var(ncid2(i), 'Hour', NF90_INT64, dimid1, varid6)
         ierr2 = nf90_def_var(ncid2(i), 'Minute', NF90_INT64, dimid1, varid7)
         ierr2 = nf90_def_var(ncid2(i), 'Second', NF90_DOUBLE, dimid1, varid8)
         ierr2 = nf90_def_var(ncid2(i), 'PBL_Height', NF90_FLOAT, dimid1, varid10)
         ierr2 = nf90_def_var(ncid2(i), 'Flag_nightday', NF90_INT64, dimid1, varid11)

         ierr2 = nf90_put_att(ncid2(i), varid9, 'units', 'YYYYMMDDHH')
         ierr2 = nf90_put_att(ncid2(i), varid1, 'units', 'degrees_north')
         ierr2 = nf90_put_att(ncid2(i), varid2, 'units', 'degrees_east')
         ierr2 = nf90_put_att(ncid2(i), varid10, 'units', 'meters')
         ierr2 = nf90_put_att(ncid2(i), varid11, 'units', '0_night_1_day')

         ierr2 = nf90_enddef(ncid2(i))

         kk = total_nobs(i)
         ierr2 = nf90_put_var(ncid2(i), varid9, ana_time(i))  
         ierr2 = nf90_put_var(ncid2(i), varid1, slat(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid2, slon(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid3, syear(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid4, smonth(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid5, sdate(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid6, shour(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid7, smin(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid8, ssec(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid10, spblh(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid11, sflag_nightday(1:kk,i))  

         ierr2 = NF90_CLOSE(ncid2(i))
      end do
  end program read_pblh_lidar_cats_cnn

  subroutine add_interval(ana_mm1, ana_dd1, ana_hh1, interval, ana_mm0, ana_dd0, ana_hh0)
     use kinds, only: i_kind
     integer(i_kind) :: ana_mm1, ana_dd1, ana_hh1, interval, ana_mm0, ana_dd0, ana_hh0
     ana_mm1 = ana_mm0
     ana_dd1 = ana_dd0
     ana_hh1 = ana_hh0 + interval
     if (ana_hh1 >= 24) then
        ana_hh1 = ana_hh1 - 24
        ana_dd1 = ana_dd0 + 1
        if (ana_mm0 == 8 .and. ana_dd1 > 31) then
           ana_mm1 = 9
           ana_dd1 = 1
        end if
        if (ana_mm0 == 9 .and. ana_dd1 > 30) then
           ana_mm1 = 10
           ana_dd1 = 1
        end if
     end if
  end subroutine add_interval

  subroutine add_oneday(mm1, dd1, mm0, dd0)
     use kinds, only: i_kind
     integer(i_kind) :: mm1, dd1, mm0, dd0
     mm1 = mm0
     dd1 = dd0+1
     if (mm0 == 8 .and. dd1 > 31) then
        mm1 = 9
        dd1 = 1
     end if
     if (mm0 == 9 .and. dd1 > 30) then
        mm1 = 10
        dd1 = 1
     end if
  end subroutine add_oneday

  subroutine djul_day_inv(jday,hour,min,sec)
    use kinds, only: r_kind,r_double,i_kind
    implicit none
    ! Convert Julian date (fractional day-of-year) -> hour,min,sec

    !Transcription of Dennis Hlavka's IDL code. Original description follows.
    !    This subroutine calculates the hour, minute, and second of the fractional
    !    day-of-year (pc) and returns the three parms to the calling program.
    !    The 'second' parameter includs a fractional part.

    ! Variable declarations
    real(kind=8),intent(in) :: jday     ! Fractional day-of-year input
    real(r_double) :: secs                ! Total seconds 
    real(r_double),intent(out) :: sec     ! Remaining seconds
    integer,intent(out) :: hour, min      ! Hours and minutes

    ! Perform calculations
    secs = jday * 86400.0      ! Convert jday to total seconds
    hour = int(secs / 3600)    ! Calculate hours
    min = int((secs - hour * 3600) / 60)  ! Calculate minutes
    sec = secs - (hour * 3600 + min * 60) ! Calculate remaining seconds (fractional)

    ! Display the results
    !print *, "Fractional day-of-year, hr, min, sec: ", jday, hour, min, sec
  end subroutine djul_day_inv

