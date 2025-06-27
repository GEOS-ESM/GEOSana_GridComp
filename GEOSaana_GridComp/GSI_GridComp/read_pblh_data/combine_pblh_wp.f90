  program read_pblh_radar_wp

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
      integer(i_kind),parameter:: ntime = 161
      integer(i_kind),parameter:: maxnum = 2500 ! 150000
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
      integer(i_kind) iret,im,iy,idd,ihh,istat,i,j,k,kk,ii
      integer(i_kind) ikx,nkx,kx,nreal,ilat,ilon,iout
      character(5) sid
      integer(i_kind) idx

      integer(i_kind) :: ncid,ncid2(ntime),ierr,ierr2,dimid1,dimid2,dimid3,dimid4,dimid5,nobs
      integer(i_kind) :: ntime_stn,nlat,nlon,nlev
      integer(i_kind) :: varid1,varid2,varid3,varid4,varid5,varid6
      integer(i_kind) :: varid7,varid8,varid9,varid10,varid11,varid12,varid13,varid14,varid15
!     integer(i_kind) :: iyear, imonth, idate, ihour, iminute
      !real(r_kind), allocatable,dimension(:) :: time, local_time
      !real(r_kind), allocatable,dimension(:) :: pblh, cloud_frac, error, negative_snr
      real(r_kind), allocatable,dimension(:,:,:,:) :: time, local_time
      real(r_kind), allocatable,dimension(:,:,:,:) :: pblh, cloud_frac, error, negative_snr
      real(r_kind), allocatable,dimension(:) :: lat, lon, lev
      integer(i_kind), allocatable,dimension(:) :: ryear, rmonth, rdate, rhour, rmin
      real(r_double), allocatable,dimension(:) :: rsec

      integer(i_kind) :: total_nobs(ntime)
      integer(i_kind) ana_time(ntime), win_time(ntime), obs_time
      real(r_kind), dimension(maxnum,ntime) :: slat, slon, slev, spblh, scloud_frac, serror,snegative_snr
      integer(i_kind), dimension(maxnum,ntime) :: syear, smonth, sdate, shour, smin
      real(r_double), dimension(maxnum,ntime) :: ssec
      character(5), dimension(maxnum,ntime) :: ssid
      !real(r_double), dimension(nheights,ntime) :: sliadr_data_alt
      !real(r_double), dimension(maxnum,nheights,ntime) :: sATB

      obstype = 'pblh'
      sis = 'wp'

!     set analysis window and output file name
      ana_yy0 = 2015
      ana_yy1 = ana_yy0
      ana_mm0 = 8
      ana_dd0 = 20
      ana_hh0 = 21
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
      print*, 'nfile = ', nfile

!     Read data from each file and Sort data into analysis windows
      total_nobs = 0 
      do k = 1, nfile ! 91 stations

!        Check if pblh file exists
         print*, "input filename = ", filename(k)
         inquire(file=trim(filename(k)),exist=lexist)
         if (.not.lexist) cycle

!        Read data
         ierr =  NF90_OPEN(trim(filename(k)),0,ncid)
         if (ierr /= nf90_noerr) call handle_err(ierr,"open")

         ierr = NF90_INQ_DIMID(ncid,'time',dimid1)
         if (ierr /= nf90_noerr) call handle_err(ierr,"ntime_stn")
         ierr = NF90_INQ_DIMID(ncid,'lev',dimid2)
         if (ierr /= nf90_noerr) call handle_err(ierr,"nlev")
         ierr = NF90_INQ_DIMID(ncid,'lat',dimid3)
         if (ierr /= nf90_noerr) call handle_err(ierr,"nlat")
         ierr = NF90_INQ_DIMID(ncid,'lon',dimid4)
         if (ierr /= nf90_noerr) call handle_err(ierr,"nlon")

         ierr = nf90_inquire_dimension(ncid, dimid1, len = ntime_stn)
         ierr = nf90_inquire_dimension(ncid, dimid2, len = nlev)
         ierr = nf90_inquire_dimension(ncid, dimid3, len = nlat)
         ierr = nf90_inquire_dimension(ncid, dimid4, len = nlon)
         print*, 'read_pblh: ntime_stn=', ntime_stn

         ! Reverse order of dimensions as stated in ncdump:
         allocate(lat(1), lon(1), lev(1), time(1,1,1,ntime_stn), pblh(1,1,1,ntime_stn) )
         allocate(local_time(1,1,1,ntime_stn))
         allocate(cloud_frac(1,1,1,ntime_stn))
         allocate(error(1,1,1,ntime_stn))
         allocate(negative_snr(1,1,1,ntime_stn))
         allocate(ryear(ntime_stn), rmonth(ntime_stn), rdate(ntime_stn), rhour(ntime_stn), rmin(ntime_stn))
         allocate(rsec(ntime_stn))

!        Latitude/Longitude and observation time
!        LATITUDE: degrees_north
         ierr = NF90_INQ_VARID(ncid,'lat',varid1)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid1,lat)
!        LONGITUDE: degrees_east
         ierr = NF90_INQ_VARID(ncid,'lon',varid2)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid2,lon)
!        lev: Altitude [meters]
         ierr = NF90_INQ_VARID(ncid,'lev',varid3)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid3,lev)
!        TIME: hours since 2009-05-17 00:00:00 (UTC)"
         ierr = NF90_INQ_VARID(ncid,'time_grads',varid4)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid4,time)
!        local_time: hours since 2009-05-17 00:00:00 (local)"
         ierr = NF90_INQ_VARID(ncid,'time_local',varid5)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid5,local_time)
!        PBL Height estimate: [meters]
!        Above ground level
         ierr = NF90_INQ_VARID(ncid,'pblh',varid6)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid6,pblh)
         !Total cloud area fraction (obtained from CERES_SYN1deg_Ed4.1 3-hourly)
         ierr = NF90_INQ_VARID(ncid,'cloud_fraction',varid7)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid7,cloud_frac)
         !2x smoothed vertical resolution of profile (error estimate)
         ierr = NF90_INQ_VARID(ncid,'vert_resolution_2x',varid8)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid8,error)
         !negative indicator (0: no negative SNR value below absolute max of SNR, 1: has negative SNR value below absolute max of SNR
         ierr = NF90_INQ_VARID(ncid,'negativeSNR_indicator',varid9)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid9,negative_snr)

!        Read the global attribute 
!        Here, string global attribute doesn't work.
         !ierr = NF90_INQUIRE_ATTRIBUTE(ncid, nf90_global, 'SITE',xtype, len, attnum)
         !if (ierr == nf90_noerr) ierr = NF90_GET_ATT_TEXT(ncid,nf90_global,'SITE',site)
         !print*, "site=",site
         !ierr = NF90_INQUIRE_ATTRIBUTE(ncid, nf90_global, 'INSTRUMENT')
         !if (ierr == nf90_noerr) ierr = NF90_GET_ATT_TEXT(ncid,nf90_global,'INSTRUMENT',instrument)
         !print*, "instrument=",instrument

         ! Read station id from file name
         filename_each = filename(k)
         idx = scan(filename_each,'_')
         filename_each = filename_each(idx+1:)
         ! station number (e.g., 3203 or 95759)
         idx = scan(filename_each,'.')
         sid= filename_each(:idx-1) 

         print*, 'sid = ', sid
         print*, 'pblh = ', pblh
         print*, 'lat,lon,lev= ', lat,lon,lev

!        Sorting based on analysis window 
         do i = 1, ntime_stn ! ntime for each station

            ! Convert hours to calendar date 
            ! hours since 2009-05-17 00:00:00 UTC
            call calendar_date(time(1,1,1,i),ryear(i),rmonth(i),rdate(i),rhour(i),rmin(i),rsec(i))
            obs_time = ryear(i)*1000000+rmonth(i)*10000+rdate(i)*100+rhour(i)
            print*, i, "obs_time=", obs_time
            kk = 0
            do j = 2, ntime ! analysis time loop
               if (obs_time>=win_time(j-1) .and. obs_time<win_time(j)) then
                  kk = j-1
                  print*, "obs_time=", obs_time, " kk=", kk
                  exit
               end if
            end do

            if (kk == 0) cycle
            total_nobs(kk) = total_nobs(kk) + 1
            ii = total_nobs(kk)
            ! slat (number, time)
            slat(ii,kk) = lat(1) 
            slon(ii,kk) = lon(1) 
            slev(ii,kk) = lev(1) 
            syear(ii,kk) = ryear(i) 
            smonth(ii,kk) = rmonth(i) 
            sdate(ii,kk) = rdate(i) 
            shour(ii,kk) = rhour(i) 
            smin(ii,kk) = rmin(i) 
            ssec(ii,kk) = rsec(i) 
            if (pblh(1,1,1,i) > 10000) then
               print*,'if pblh>10000 -> pblh=',pblh(1,1,1,i)
               spblh(ii,kk) = -9999.0
               scloud_frac(ii,kk) = -9999.0
               serror(ii,kk) = -9999.0
               snegative_snr(ii,kk) = -9999.0
            else
               spblh(ii,kk) = pblh(1,1,1,i)
               scloud_frac(ii,kk) = cloud_frac(1,1,1,i)
               serror(ii,kk) = error(1,1,1,i)
               snegative_snr(ii,kk) = negative_snr(1,1,1,i)
            end if
            ssid(ii,kk) = sid

         end do

         ierr = NF90_CLOSE(ncid)
         if (ierr /= nf90_noerr) call handle_err(ierr,"close")

         deallocate(lat, lon, lev, time, local_time, pblh, cloud_frac, error, negative_snr)
         deallocate(ryear, rmonth, rdate, rhour, rmin, rsec)

      end do ! end of nfile

!     Write into netcdf files for each analysis window
      do i = 1, ntime
         if (total_nobs(i) < 1) cycle

         write(filename2(i), '(A13,I10,A4)') 'wp_pblh.', ana_time(i), '.nc4'
!        ierr2 =  NF90_CREATE(trim(filename2(i)), NF90_CLOBBER, ncid2(i)) ! overwrite this file if it already exists
         ierr2 =  NF90_CREATE(trim(filename2(i)), NF90_NETCDF4, ncid2(i)) 
         if (ierr2 /= nf90_noerr) call handle_err(ierr2,"create")
!        ierr2 = nf90_open(trim(filename2(i)),nf90_write, ncid2(i))
!        if (ierr2 /= nf90_noerr) call handle_err(ierr2,"open")
         print*, "filename2 = ", filename2(i)
         print*, "total_nobs =", total_nobs(i)

         ierr2 = nf90_def_dim(ncid2(i), 'Time', 1, dimid3)
         ierr2 = nf90_def_dim(ncid2(i), 'nobs', total_nobs(i), dimid1)
         ierr2 = nf90_def_dim(ncid2(i), 'Sid_maxstrlen', 5, dimid4)

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
         ierr2 = nf90_def_var(ncid2(i), 'Altitude', NF90_FLOAT, dimid1, varid11)
         ierr2 = nf90_def_var(ncid2(i), 'Cloud_Fraction', NF90_FLOAT, dimid1, varid12)
         ierr2 = nf90_def_var(ncid2(i), 'Vert_Resolution_2x', NF90_FLOAT, dimid1, varid13)
         ierr2 = nf90_def_var(ncid2(i), 'Negative_SNR_Indicator', NF90_FLOAT, dimid1, varid14)
         ierr2 = nf90_def_var(ncid2(i), 'Station_ID', NF90_CHAR,(/dimid4, dimid1/), varid15)

         ierr2 = nf90_put_att(ncid2(i), varid9, 'units', 'YYYYMMDDHH')
         ierr2 = nf90_put_att(ncid2(i), varid1, 'units', 'degrees_north')
         ierr2 = nf90_put_att(ncid2(i), varid2, 'units', 'degrees_east')
         ierr2 = nf90_put_att(ncid2(i), varid10, 'units', 'meters')
         ierr2 = nf90_put_att(ncid2(i), varid11, 'units', 'meters')
         !ierr2 = nf90_put_att(ncid2(i), varid10, 'missing_value', 1.e+15)
         ierr2 = nf90_put_att(ncid2(i), varid10, 'missing_value', -9999.0)
         ierr2 = nf90_put_att(ncid2(i), varid12, 'long_name', 'total cloud area fraction (obtained from CERES_SYN1deg_Ed4.1 3-hourly)')
         ierr2 = nf90_put_att(ncid2(i), varid13, 'long_name', '2x smoothed vertical resolution of profile (error estimate)')
         ierr2 = nf90_put_att(ncid2(i), varid14, 'long_name', 'negative indicator (0: no negative SNR value below absolute max of SNR, 1: has negative SNR value below absolute max of SNR')

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
         ierr2 = nf90_put_var(ncid2(i), varid11, slev(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid12, scloud_frac(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid13, serror(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid14, snegative_snr(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid15, ssid(1:kk,i))

         ierr2 = NF90_CLOSE(ncid2(i))
      end do
  end program read_pblh_radar_wp

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

  subroutine calendar_date(org_hr,year,month,day,hrs,min,sec)

    ! Returns the year, month, day, hr, min, sec for the specified date.

    use kinds, only: r_kind,r_double,i_kind

    implicit none

    real(r_kind),intent(in)  :: org_hr ! hours since 2009-05-17 00:00:00 UTC
    real(r_double) :: julian_date ! 
    real(r_double) :: julian_date_2009_05_17 ! 2454968.5 = 2009-05-17 00 UTC 
    integer,intent(out)  :: year
    integer,intent(out)  :: month
    integer,intent(out)  :: day
    integer,intent(out)  :: hrs
    integer,intent(out)  :: min
    logical :: is_leap_year

    real(r_double),intent(out) :: sec

    integer :: i,j,k,l,n,jd
    real(r_double) :: frac_day

    print*, 'org_hr [hours since 2009-05-17 00 UTC]=',org_hr
    ! The calendar date for 2454968.5 is 2009-May-17, 00 UTC.
    julian_date_2009_05_17 = 2454968.5
    julian_date = real(org_hr,kind=r_double)/24.0 + julian_date_2009_05_17 
    jd = int(julian_date)

    l = jd+68569
    n = 4*l/146097
    l = l-(146097*n+3)/4
    i = 4000*(l+1)/1461001
    l = l-1461*i/4+31
    j = 80*l/2447
    k = l-2447*j/80
    l = j/11
    j = j+2-12*l
    i = 100*(n-49)+i+l

    year  = i
    month = j
    day   = k

    frac_day = julian_date - real(jd,kind=r_double) + 0.5
    hrs = int(frac_day*24.0)
    min = int((frac_day - hrs/24.0) * 1440.0)
    sec = (frac_day - hrs/24.0 - min/1440.0) * 86400.0

    if (sec >= 60.0) then
       sec = sec - 60.0
       min = min + 1
    end if

    if (min >= 60) then
       min = min - 60
       hrs = hrs + 1
    end if

    if (hrs >= 24.0) then
       hrs = hrs - 24
       day = day + 1
    end if

    ! check if day is later than the last day of the month depending on the month
    if (month==1 .or. month==3 .or. month==5 .or. month==7 .or. month==8 .or. month==10 .or. month==12 ) then
       if (day > 31) then
          month = month + 1
          day = 1
          if (month > 12) then
             month = 1
             year = year+1
          end if
       end if
    end if

    if (month==4 .or. month==6 .or. month==9 .or. month==11 ) then
       if (day > 30) then
          month = month + 1
          day = 1
       end if
    end if

    ! check if this is a leap year for Feb
    if (month==2) then 
       ! leap year 
       is_leap_year = .false.
       if ( mod(year,4)==0) then
          if ( mod(year,100)==0) then
             if ( mod(year,400)==0) then
                is_leap_year = .true.
             end if
          else
             is_leap_year = .true.
          end if
       end if

       if (is_leap_year) then
          if (day > 29) then
             month = 3
          end if 
       else 
          if (day > 28) then
             month = 3
          end if 
       end if 

    end if

  end subroutine calendar_date
