  program read_pblh_lidar_calipso

! !USES:

      use kinds, only: r_double,r_double,i_kind
      use constants, only: zero,one_tenth,one,deg2rad,rad2deg,three
      use gridmod, only: diagnostic_reg,regional,nlon,nlat,&
           tll2xy,txy2ll,rotate_wind_ll2xy,rotate_wind_xy2ll,&
           rlats,rlons
      use convinfo, only: nconvtype,ctwind, &
           icuse,ictype,ioctype
      use gsi_4dvar, only: l4dvar,l4densvar,time_4dvar,winlen
      use obsmod, only: iadate,offtime_data,bmiss
      use deter_sfc_mod, only: deter_sfc2
      use mpimod, only: npe
      use netcdf

      implicit none

!     Declare local parameters
      integer(i_kind),parameter:: ntime = 40
      integer(i_kind),parameter:: maxnum = 150000
      integer(i_kind),parameter:: nheights = 583
      real(r_double),parameter:: r360 = 360.0_r_double
      real(r_double),parameter:: missing = -99999.9_r_double

      character(len=225):: filename(32767),filename2(ntime),obstype
      character(20):: sis

      real(r_double),allocatable,dimension(:,:):: cdata_all

      character(10) date
      character*255 argv
      logical first,outside,inflate_error,lexist,more_data
      integer(i_kind) nfile, iarg, argc, iargc
      integer(i_kind) ana_yy1,ana_mm1,ana_dd1,ana_hh1,ana_yy0,ana_mm0,ana_dd0,ana_hh0
      integer(i_kind) iret,im,iy,idd,ihh,istat,i,j,k,kk,ii
      integer(i_kind) ikx,nkx,kx,nreal,ilat,ilon,iout
      integer(i_kind) ntest,nchanl
      integer(i_kind) pblhqm,maxobs,idomsfc
!     integer(i_kind),dimension(8):: obs_time,anal_time,lobs_time
!     real(r_double) ltime,cenlon_tmp
      real(r_double) usage
      real(r_double) pblhob,pblhoe,pblhelev
      real(r_double) dlat,dlon,dlat_earth,dlon_earth
      real(r_double) dlat_earth_deg,dlon_earth_deg
      real(r_double) cdist,disterr,disterrmax,rlon00,rlat00
      real(r_double) :: tsavg,ff10,sfcr,zz
      real(r_double) dx,dy,dx1,dy1,w00,w10,w01,w11

      integer(i_kind) idate5(5),minobs,minan
      real(r_double) time_correction,timeobs,time,toff,t4dv,zeps

!     real(r_single) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblh_calipso
      real(r_double) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblh_calipso

      integer(i_kind) :: ncid,ncid2(ntime),ierr,ierr2,dimid1,dimid2,dimid3,norbits
      integer(i_kind) :: varid1,varid2,varid3,varid4,varid5,varid6
      integer(i_kind) :: varid7,varid8,varid9,varid10,varid11,varid12,varid13
!     integer(i_kind) :: iyear, imonth, idate, ihour, iminute
      real(r_double), allocatable,dimension(:) :: lat, lon, pblh, sfc_elev, sfc_mask, liadr_data_alt
      real(r_double), allocatable,dimension(:) :: ryear, rmonth, rdate, rhour, rminute
      real(r_double), allocatable,dimension(:,:) :: ATB

      integer(i_kind) :: total_norbits(ntime)
      integer(i_kind) ana_time(ntime), win_time(ntime), obs_time
      real(r_double), dimension(maxnum,ntime) :: slat, slon, spblh, ssfc_elev, ssfc_mask
      real(r_double), dimension(maxnum,ntime) :: syear, smonth, sdate, shour, sminute
      real(r_double), dimension(nheights,ntime) :: sliadr_data_alt
      real(r_double), dimension(maxnum,nheights,ntime) :: sATB

      obstype = 'pblh'
      sis = 'calipso'

!     set analysis window and output file name
      ana_yy0 = 2015
      ana_yy1 = ana_yy0
      ana_mm0 = 9
      ana_dd0 = 29
      ana_hh0 = 3
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
      total_norbits = 0 
      do k = 1, nfile

!        Check if pblh file exists
         print*, "input filename = ", filename(k)
         inquire(file=trim(filename(k)),exist=lexist)
         if (.not.lexist) cycle

!        Read data
         ierr =  NF90_OPEN(trim(filename(k)),0,ncid)
         if (ierr /= nf90_noerr) call handle_err(ierr,"open")

         ierr = NF90_INQ_DIMID(ncid,'norbits',dimid1)
         if (ierr /= nf90_noerr) call handle_err(ierr,"norbits")
         ierr = NF90_INQ_DIMID(ncid,'nheights',dimid2)
         if (ierr /= nf90_noerr) call handle_err(ierr,"nheights")

         ierr = nf90_inquire_dimension(ncid, dimid1, len = norbits)
!        ierr = nf90_inquire_dimension(ncid, dimid2, len = nheights)
         print*, 'read_pblh: norbits=', norbits, ' nheights=', nheights

         allocate(lat(norbits), lon(norbits), pblh(norbits), sfc_elev(norbits))
         allocate(sfc_mask(norbits), liadr_data_alt(nheights), ATB(norbits, nheights))
         allocate(ryear(norbits), rmonth(norbits), rdate(norbits), rhour(norbits), rminute(norbits))

!        Latitude/Longitude and observation time
         ierr = NF90_INQ_VARID(ncid,'lat',varid1)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid1,lat)
         ierr = NF90_INQ_VARID(ncid,'lon',varid2)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid2,lon)
         ierr = NF90_INQ_VARID(ncid,'Year',varid3)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid3,ryear)
         ierr = NF90_INQ_VARID(ncid,'Month',varid4)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid4,rmonth)
         ierr = NF90_INQ_VARID(ncid,'Date',varid5)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid5,rdate)
         ierr = NF90_INQ_VARID(ncid,'Hour',varid6)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid6,rhour)
         ierr = NF90_INQ_VARID(ncid,'Minute',varid7)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid7,rminute)
!        PBL_Height: meters
         ierr = NF90_INQ_VARID(ncid,'PBL_Height',varid8)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid8,pblh)
!        SurfaceElevation: meters
         ierr = NF90_INQ_VARID(ncid,'SurfaceElevation',varid9)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid9,sfc_elev)
!        0=shallow ocean, 1=land, 2=coastlines, 3=shallow inland water, 4=intermittent water, 
!        5=deep inland water, 6=continental ocean, 7=deep ocean
         ierr = NF90_INQ_VARID(ncid,'Land_Water_Mask',varid10)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid10,sfc_mask)
!        Lidar_Data_Altitude: km
         ierr = NF90_INQ_VARID(ncid,'Lidar_Data_Altitude',varid11)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid11,liadr_data_alt)
!        Total_Attenuated_Backscatter_532: km^-1 sr^-1
         ierr = NF90_INQ_VARID(ncid,'Total_Attenuated_Backscatter_532',varid12)
         if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid12,ATB)

!        print*, 'pblh = ', pblh

!        Sorting based on analysis window 
         do i = 1, norbits
            obs_time = ryear(i)*1000000+rmonth(i)*10000+rdate(i)*100+rhour(i)
!           print*, i, "obs_time=", obs_time
            kk = 0
            do j = 2, ntime
               if (obs_time>=win_time(j-1) .and. obs_time<win_time(j)) then
                  kk = j-1
                  print*, "obs_time=", obs_time, " kk=", kk
                  exit
               end if
            end do

            if (kk == 0) cycle
            total_norbits(kk) = total_norbits(kk) + 1
            ii = total_norbits(kk)
            slat(ii,kk) = lat(i) 
            slon(ii,kk) = lon(i) 
            syear(ii,kk) = ryear(i) 
            smonth(ii,kk) = rmonth(i) 
            sdate(ii,kk) = rdate(i) 
            shour(ii,kk) = rhour(i) 
            sminute(ii,kk) = rminute(i) 
            spblh(ii,kk) = pblh(i) 
            ssfc_elev(ii,kk) = sfc_elev(i) 
            ssfc_mask(ii,kk) = sfc_mask(i)
            do j = 1, nheights 
               sliadr_data_alt(j,kk) = liadr_data_alt(j)
               sATB(ii,j,kk) = ATB(i,j)
            end do
         end do

         ierr = NF90_CLOSE(ncid)
         if (ierr /= nf90_noerr) call handle_err(ierr,"close")

         deallocate(lat, lon, pblh, sfc_elev)
         deallocate(sfc_mask, liadr_data_alt, ATB)
         deallocate(ryear, rmonth, rdate, rhour, rminute)

      end do ! end of nfile

!     Write into netcdf files for each analysis window
      do i = 1, ntime
         if (total_norbits(i) < 1) cycle

         write(filename2(i), '(A13,I10,A4)') 'calipso_pblh.', ana_time(i), '.nc4'
!        ierr2 =  NF90_CREATE(trim(filename2(i)), NF90_CLOBBER, ncid2(i)) ! overwrite this file if it already exists
         ierr2 =  NF90_CREATE(trim(filename2(i)), NF90_NETCDF4, ncid2(i)) 
         if (ierr2 /= nf90_noerr) call handle_err(ierr2,"create")
!        ierr2 = nf90_open(trim(filename2(i)),nf90_write, ncid2(i))
!        if (ierr2 /= nf90_noerr) call handle_err(ierr2,"open")
         print*, "filename2 = ", filename2(i)
         print*, "total_norbits =", total_norbits(i)

         ierr2 = nf90_def_dim(ncid2(i), 'Time', 1, dimid3)
         ierr2 = nf90_def_dim(ncid2(i), 'norbits', total_norbits(i), dimid1)
         ierr2 = nf90_def_dim(ncid2(i), 'nheights', nheights, dimid2)

         ierr2 = nf90_def_var(ncid2(i), 'Ana_Time', NF90_INT64, dimid3, varid13)
         ierr2 = nf90_def_var(ncid2(i), 'lat', NF90_DOUBLE, dimid1, varid1)
         ierr2 = nf90_def_var(ncid2(i), 'lon', NF90_DOUBLE, dimid1, varid2)
         ierr2 = nf90_def_var(ncid2(i), 'Year', NF90_DOUBLE, dimid1, varid3)
         ierr2 = nf90_def_var(ncid2(i), 'Month', NF90_DOUBLE, dimid1, varid4)
         ierr2 = nf90_def_var(ncid2(i), 'Date', NF90_DOUBLE, dimid1, varid5)
         ierr2 = nf90_def_var(ncid2(i), 'Hour', NF90_DOUBLE, dimid1, varid6)
         ierr2 = nf90_def_var(ncid2(i), 'Minute', NF90_DOUBLE, dimid1, varid7)
         ierr2 = nf90_def_var(ncid2(i), 'PBL_Height', NF90_DOUBLE, dimid1, varid8)
         ierr2 = nf90_def_var(ncid2(i), 'SurfaceElevation', NF90_DOUBLE, dimid1, varid9)
         ierr2 = nf90_def_var(ncid2(i), 'Land_Water_Mask', NF90_DOUBLE, dimid1, varid10)
         ierr2 = nf90_def_var(ncid2(i), 'Lidar_Data_Altitude', NF90_DOUBLE, dimid2, varid11)
         ierr2 = nf90_def_var(ncid2(i), 'Total_Attenuated_Backscatter_532', NF90_DOUBLE, (/ dimid1,dimid2 /), varid12)

         ierr2 = nf90_put_att(ncid2(i), varid13, 'units', 'YYYYMMDDHH')
         ierr2 = nf90_put_att(ncid2(i), varid1, 'units', 'degrees_north')
         ierr2 = nf90_put_att(ncid2(i), varid2, 'units', 'degrees_east')
         ierr2 = nf90_put_att(ncid2(i), varid8, 'units', 'meters')
         ierr2 = nf90_put_att(ncid2(i), varid9, 'units', 'meters')
!        ierr2 = nf90_put_att(ncid2(i), varid10, 'interp', '0=shallow ocean, 1=land, 2=coastlines, 3=shallow inland water, 4=intermittent water, 5=deep inland water, 6=continental ocean, 7=deep ocean')
         ierr2 = nf90_put_att(ncid2(i), varid11, 'units', 'km')
         ierr2 = nf90_put_att(ncid2(i), varid12, 'units', 'km^-1 sr^-1')

         ierr2 = nf90_enddef(ncid2(i))

         kk = total_norbits(i)
         ierr2 = nf90_put_var(ncid2(i), varid13, ana_time(i))  
         ierr2 = nf90_put_var(ncid2(i), varid1, slat(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid2, slon(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid3, syear(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid4, smonth(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid5, sdate(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid6, shour(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid7, sminute(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid8, spblh(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid9, ssfc_elev(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid10, ssfc_mask(1:kk,i))  
         ierr2 = nf90_put_var(ncid2(i), varid11, sliadr_data_alt(1:nheights,i))  
         ierr2 = nf90_put_var(ncid2(i), varid12, sATB(1:kk,1:nheights,i))  

         ierr2 = NF90_CLOSE(ncid2(i))
      end do
  end program read_pblh_lidar_calipso

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

