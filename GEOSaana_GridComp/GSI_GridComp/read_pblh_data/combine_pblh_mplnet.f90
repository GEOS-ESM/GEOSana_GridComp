  program read_pblh_lidar_mplnet

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
      character(19) site
      character(8) instrument
      integer(i_kind) idx

      integer(i_kind) :: ncid,ncid2(ntime),ierr,ierr2,dimid1,dimid2,dimid3,dimid4,dimid5,nobs
      integer(i_kind) :: ntime_1440
      integer(i_kind) :: varid1,varid2,varid3,varid4,varid5,varid6
      integer(i_kind) :: varid7,varid8,varid9,varid10,varid11,varid12,varid13,varid14
!     integer(i_kind) :: iyear, imonth, idate, ihour, iminute
      real(r_double), allocatable,dimension(:) :: time
      real(r_kind), allocatable,dimension(:) :: lat, lon, pblh, sfc_alti, flag_cloud_screen
      integer(i_kind), allocatable,dimension(:) :: ryear, rmonth, rdate, rhour, rmin
      real(r_double), allocatable,dimension(:) :: rsec
      !real(r_double), allocatable,dimension(:,:) :: ATB

      integer(i_kind) :: total_nobs(ntime)
      integer(i_kind) ana_time(ntime), win_time(ntime), obs_time
      real(r_kind), dimension(maxnum,ntime) :: slat, slon, spblh, ssfc_alti, sflag_cloud_screen
      integer(i_kind), dimension(maxnum,ntime) :: syear, smonth, sdate, shour, smin
      real(r_double), dimension(maxnum,ntime) :: ssec
      character(19), dimension(maxnum,ntime) :: ssite
      character(8), dimension(maxnum,ntime) :: sinstrument
      !real(r_double), dimension(nheights,ntime) :: sliadr_data_alt
      !real(r_double), dimension(maxnum,nheights,ntime) :: sATB

      obstype = 'pblh'
      sis = 'mplnet'

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
      do k = 1, nfile

!        Check if pblh file exists
         print*, "input filename = ", filename(k)
         inquire(file=trim(filename(k)),exist=lexist)
         if (.not.lexist) cycle

!        Read data
         ierr =  NF90_OPEN(trim(filename(k)),0,ncid)
         if (ierr /= nf90_noerr) call handle_err(ierr,"open")

         ierr = NF90_INQ_DIMID(ncid,'time',dimid1)
         if (ierr /= nf90_noerr) call handle_err(ierr,"ntime_1440")

         ierr = nf90_inquire_dimension(ncid, dimid1, len = ntime_1440)
         print*, 'read_pblh: ntime_1440=', ntime_1440

         allocate(lat(ntime_1440), lon(ntime_1440), time(ntime_1440), pblh(ntime_1440), sfc_alti(ntime_1440))
         allocate(flag_cloud_screen(ntime_1440))
         !allocate(lat(time), lon(time), pblh(time), sfc_elev(norbits))
         !allocate(sfc_mask(norbits), liadr_data_alt(nheights), ATB(norbits, nheights))
         allocate(ryear(ntime_1440), rmonth(ntime_1440), rdate(ntime_1440), rhour(ntime_1440), rmin(ntime_1440))
         allocate(rsec(ntime_1440))

!        Latitude/Longitude and observation time
!        LATITUDE: degrees_north
         ierr = NF90_INQ_VARID(ncid,'LATITUDE',varid1)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid1,lat)
!        LONGITUDE: degrees_east
         ierr = NF90_INQ_VARID(ncid,'LONGITUDE',varid2)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid2,lon)
!        TIME: JULIAN (UTC) days since -4713-01-01 12:00:00 UTC
         ierr = NF90_INQ_VARID(ncid,'TIME',varid3)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid3,time)
         !if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid3,ryear)
         !if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid4,rmonth)
         !if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid5,rdate)
         !if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid6,rhour)
         !if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid7,rminute)
!        PBL Height: [km]
!        height above surface at the top of mixed layer (planetary boundary layer or mixing layer)
         ierr = NF90_INQ_VARID(ncid,'MIXED_LAYER_HEIGHT',varid4)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid4,pblh)
!        Surface Altitude: transceiver height above sea level [km]
         ierr = NF90_INQ_VARID(ncid,'SURFACE_ALTITUDE',varid5)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid5,sfc_alti)
!        FLAG_CLOUD_SCREEN:
!        1b, 2b, 4b, 8b ;
!        FLAG_MEANINGS = "cloud_free", "cloud_fraction_>_0%", "cloud_detection_fail", "no_cloud_product" ;
         ierr = NF90_INQ_VARID(ncid,'FLAG_CLOUD_SCREEN',varid6)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid6,flag_cloud_screen)

!        Read the global attribute 
!        Here, string global attribute doesn't work.
         !ierr = NF90_INQUIRE_ATTRIBUTE(ncid, nf90_global, 'SITE',xtype, len, attnum)
         !if (ierr == nf90_noerr) ierr = NF90_GET_ATT_TEXT(ncid,nf90_global,'SITE',site)
         !print*, "site=",site
         !ierr = NF90_INQUIRE_ATTRIBUTE(ncid, nf90_global, 'INSTRUMENT')
         !if (ierr == nf90_noerr) ierr = NF90_GET_ATT_TEXT(ncid,nf90_global,'INSTRUMENT',instrument)
         !print*, "instrument=",instrument

         ! Read site and instrument from file name
         filename_each = filename(k)
         idx = scan(filename_each,'_')
         filename_each = filename_each(idx+1:)
         do while (filename_each(1:3)/='MPL')
            idx = scan(filename_each,'_')
            filename_each= filename_each(idx+1:)
         end do
         ! instrument (e.g., MPL44104)
         idx = scan(filename_each,'_')
         instrument=filename_each(:idx-1)
         ! site (e.g. GSFC, Singapore)
         filename_each = filename_each(idx+1:)
         idx = scan(filename_each,'.')
         site = filename_each(:idx-1) 

         print*, 'site, instrument, pblh = ', site, instrument, pblh

!        Sorting based on analysis window 
         do i = 1, ntime_1440

            ! Convert Julian date 
            call calendar_date_realsec(time(i),ryear(i),rmonth(i),rdate(i),rhour(i),rmin(i),rsec(i))
            obs_time = ryear(i)*1000000+rmonth(i)*10000+rdate(i)*100+rhour(i)
            print*, i, "obs_time=", obs_time
            kk = 0
            do j = 2, ntime
               if (obs_time>=win_time(j-1) .and. obs_time<win_time(j)) then
                  kk = j-1
                  print*, "obs_time=", obs_time, " kk=", kk
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
            spblh(ii,kk) = pblh(i)*1000.0 
            ssfc_alti(ii,kk) = sfc_alti(i)*1000.0 
            sflag_cloud_screen(ii,kk) = flag_cloud_screen(i)
            ssite(ii,kk) = site
            sinstrument(ii,kk) = instrument

            !do j = 1, nheights 
            !   sliadr_data_alt(j,kk) = liadr_data_alt(j)
            !   sATB(ii,j,kk) = ATB(i,j)
            !end do
         end do

         ierr = NF90_CLOSE(ncid)
         if (ierr /= nf90_noerr) call handle_err(ierr,"close")

         deallocate(lat, lon, time, pblh, sfc_alti, flag_cloud_screen)
         deallocate(ryear, rmonth, rdate, rhour, rmin, rsec)

      end do ! end of nfile

!     Write into netcdf files for each analysis window
      do i = 1, ntime
         if (total_nobs(i) < 1) cycle

         write(filename2(i), '(A13,I10,A4)') 'mplnet_pblh.', ana_time(i), '.nc4'
!        ierr2 =  NF90_CREATE(trim(filename2(i)), NF90_CLOBBER, ncid2(i)) ! overwrite this file if it already exists
         ierr2 =  NF90_CREATE(trim(filename2(i)), NF90_NETCDF4, ncid2(i)) 
         if (ierr2 /= nf90_noerr) call handle_err(ierr2,"create")
!        ierr2 = nf90_open(trim(filename2(i)),nf90_write, ncid2(i))
!        if (ierr2 /= nf90_noerr) call handle_err(ierr2,"open")
         print*, "filename2 = ", filename2(i)
         print*, "total_nobs =", total_nobs(i)

         ierr2 = nf90_def_dim(ncid2(i), 'Time', 1, dimid3)
         ierr2 = nf90_def_dim(ncid2(i), 'nobs', total_nobs(i), dimid1)
         ierr2 = nf90_def_dim(ncid2(i), 'Site_maxstrlen', 19, dimid4)
         ierr2 = nf90_def_dim(ncid2(i), 'Instrument_maxstrlen', 8, dimid5)

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
         ierr2 = nf90_def_var(ncid2(i), 'Surface_Altitude', NF90_FLOAT, dimid1, varid11)
         ierr2 = nf90_def_var(ncid2(i), 'Flag_Cloud_Screen', NF90_FLOAT, dimid1, varid12)
         ierr2 = nf90_def_var(ncid2(i), 'Site', NF90_CHAR,(/dimid4, dimid1/), varid13)
         ierr2 = nf90_def_var(ncid2(i), 'Instrument', NF90_CHAR,(/dimid5, dimid1/), varid14)
         !ierr2 = nf90_def_var(ncid2(i), 'Total_Attenuated_Backscatter_532', NF90_DOUBLE, (/ dimid1,dimid2 /), varid12)

         ierr2 = nf90_put_att(ncid2(i), varid9, 'units', 'YYYYMMDDHH')
         ierr2 = nf90_put_att(ncid2(i), varid1, 'units', 'degrees_north')
         ierr2 = nf90_put_att(ncid2(i), varid2, 'units', 'degrees_east')
         ierr2 = nf90_put_att(ncid2(i), varid10, 'units', 'meters')
         ierr2 = nf90_put_att(ncid2(i), varid11, 'units', 'meters')

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
         ierr2 = nf90_put_var(ncid2(i), varid11, ssfc_alti(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid12, sflag_cloud_screen(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid13, ssite(1:kk,i))
         ierr2 = nf90_put_var(ncid2(i), varid14, sinstrument(1:kk,i))
         !ierr2 = nf90_put_var(ncid2(i), varid12, sATB(1:kk,1:nheights,i))  

         ierr2 = NF90_CLOSE(ncid2(i))
      end do
  end program read_pblh_lidar_mplnet

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

  subroutine calendar_date_realsec(julian_date,year,month,day,hrs,min,sec)

    ! Returns the year, month, day, hr, min, sec for the specified Julian date.
    ! input unit: days since -4713-01-01 12:00:00 UTC
    !
    ! https://jacobwilliams.github.io/Fortran-Astrodynamics-Toolkit/proc/calendar_date_realsec.html

    use kinds, only: r_kind,r_double,i_kind

    implicit none

    real(r_double),intent(in)  :: julian_date !! julian date
    ! unit: days since -4713-01-01 12:00:00 UTC
    integer,intent(out)  :: year
    integer,intent(out)  :: month
    integer,intent(out)  :: day
    integer,intent(out)  :: hrs
    integer,intent(out)  :: min

    real(r_double),intent(out) :: sec

    integer :: i,j,k,l,n,jd
    real(r_double) :: frac_day

    print*, julian_date
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

    if (month==8 .and. day > 31) then
       month = 9
       day = 1
    end if

    if (month==9 .and. day > 30) then
       month = 10
       day = 1
    end if

  end subroutine calendar_date_realsec
