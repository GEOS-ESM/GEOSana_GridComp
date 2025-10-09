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

!     Declare passed variables
      character(len=225):: infile,obstype
      character(20):: sis

!     Declare local parameters
      real(r_double),parameter:: r360 = 360.0_r_double

      real(r_double),allocatable,dimension(:,:):: cdata_all

      character(10) date
      logical first,outside,inflate_error,lexist,more_data
      integer(i_kind) iret,im,iy,idd,ihh,istat,i,j,k
      integer(i_kind) ikx,nkx,kx,nreal,ilat,ilon,iout
      integer(i_kind) kk,klon1,klat1,klonp1,klatp1
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

      integer(i_kind) :: ncid,ierr,dimid1,dimid2,norbits,nheights
      integer(i_kind) :: varid1,varid2,varid3,varid4,varid5,varid6
      integer(i_kind) :: varid7,varid8,varid9,varid10,varid11,varid12
      integer(i_kind) :: iyear, imonth, idate, ihour, iminute
      real(r_double), allocatable,dimension(:) :: lat, lon, pblh, sfc_elev, sfc_mask, liadr_data_alt
      real(r_double), allocatable,dimension(:) :: ryear, rmonth, rdate, rhour, rminute
      real(r_double), allocatable,dimension(:,:) :: ATB

      obstype = 'pblh'
      sis = 'calipso'
      infile = 'CAL_LID_L1-ValStage1-V3-30.2015-07-31T23-55-36ZD.nc4'

!     Check if pblh file exists
      inquire(file=trim(infile),exist=lexist)
      if (.not.lexist) stop

!     Read data
      ierr =  NF90_OPEN(trim(infile),0,ncid)
      if (ierr /= nf90_noerr) call handle_err(ierr,"open")

      ierr = NF90_INQ_DIMID(ncid,"norbits",dimid1)
      if (ierr /= nf90_noerr) call handle_err(ierr,"norbits")
      ierr = NF90_INQ_DIMID(ncid,"nheights",dimid2)
      if (ierr /= nf90_noerr) call handle_err(ierr,"nheights")

      ierr = nf90_inquire_dimension(ncid, dimid1, len = norbits)
      ierr = nf90_inquire_dimension(ncid, dimid2, len = nheights)
      print*, 'read_pblh: norbits=', norbits, ' nheights=', nheights

      allocate(lat(norbits), lon(norbits), pblh(norbits), sfc_elev(norbits))
      allocate(sfc_mask(norbits), liadr_data_alt(nheights), ATB(norbits, nheights))
      allocate(ryear(norbits), rmonth(norbits), rdate(norbits), rhour(norbits), rminute(norbits))

!     Latitude: degrees
      ierr = NF90_INQ_VARID(ncid,"lat",varid1)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid1,lat)
      ierr = NF90_INQ_VARID(ncid,"lon",varid2)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid2,lon)
      ierr = NF90_INQ_VARID(ncid,"Year",varid3)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid3,ryear)
      ierr = NF90_INQ_VARID(ncid,"Month",varid4)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid4,rmonth)
      ierr = NF90_INQ_VARID(ncid,"Date",varid5)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid5,rdate)
      ierr = NF90_INQ_VARID(ncid,"Hour",varid6)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid6,rhour)
      ierr = NF90_INQ_VARID(ncid,"Minute",varid7)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid7,rminute)
!     PBL_Height: meters
      ierr = NF90_INQ_VARID(ncid,"PBL_Height",varid8)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid8,pblh)
!     SurfaceElevation: meters
      ierr = NF90_INQ_VARID(ncid,"SurfaceElevation",varid9)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid9,sfc_elev)
!     0=shallow ocean, 1=land, 2=coastlines, 3=shallow inland water, 4=intermittent water, 
!     5=deep inland water, 6=continental ocean, 7=deep ocean
      ierr = NF90_INQ_VARID(ncid,"Land_Water_Mask",varid10)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid10,sfc_mask)
!     Lidar_Data_Altitude: km
      ierr = NF90_INQ_VARID(ncid,"Lidar_Data_Altitude",varid11)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid11,liadr_data_alt)
!     Total_Attenuated_Backscatter_532: km^-1 sr^-1
      ierr = NF90_INQ_VARID(ncid,"Total_Attenuated_Backscatter_532",varid12)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid12,ATB)

      ierr = NF90_CLOSE(ncid)
      if (ierr /= nf90_noerr) call handle_err(ierr,"close")

      print*, 'pblh=', pblh

      first = .true.
      maxobs=0
      do i = 1, 20

!        Time offset         
         iyear = ryear(i)
         imonth = rmonth(i)
         idate = rdate(i)
         ihour = rhour(i)
         write (date,'(i4,3i2)') iyear,imonth,idate,ihour
         read (date,'( i10)') idate
         print*, 'idate=', idate
         if (first) then
            call time_4dvar(idate,toff)
            first=.false.
         end if

         print*, iyear,imonth,idate,ihour
         print*, 'read_pblh: idate, iadate= ', idate, iadate
         print*, 'read_pblh:', lat(i),lon(i),pblh(i),sfc_elev(i)

      end do

  end program read_pblh_lidar_calipso
