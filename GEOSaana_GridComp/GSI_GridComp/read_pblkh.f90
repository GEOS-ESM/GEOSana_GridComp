!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: gsimod  ---

!
! !INTERFACE:
!
  module read_pblkh

! !USES:

      use kinds, only: r_kind,r_double,i_kind
      use constants, only: zero,one_tenth,one,deg2rad,rad2deg,three,r_missing
      use gridmod, only: diagnostic_reg,regional,nlon,nlat,&
           tll2xy,txy2ll,rotate_wind_ll2xy,rotate_wind_xy2ll,&
           rlats,rlons
      use convinfo, only: nconvtype,ctwind, &
           icuse,ictype,ioctype
      use gsi_4dvar, only: l4dvar,l4densvar,time_4dvar,winlen
      use obsmod, only: iadate,offtime_data,bmiss
      use deter_sfc_mod, only: deter_sfc2, deter_sfc_type_coast
      use deter_sfc_mod, only: deter_hsstdv_model
      use mpimod, only: npe

      implicit none

      private

! !PUBLIC ROUTINES:
!     public read_pblkh

!  interface read_pblkh
!    module procedure read_pblkh_wp_text
!    module procedure read_pblkh_calipso
!  end interface
      public read_pblkh_wp_text 
      public read_pblkh_calipso
      public read_pblkh_gnssro ! to assimilate GNSSRO_HIGH

!---------------------------------------------------------------------------

   CONTAINS

     subroutine read_pblkh_gnssro(nread,ndata,nodata,infile,obstype,lunout,twindin,&
         sis,nobs)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  read_pblkh_gnssro      read obs from msgs in NETCDF files
!
! program history log:
!   2022-10     Y. Zhu     - adapted from read_pblh
!   2023-07     E.-G. Yang - subtracted surface height from gnss ro pblh,
!                            in order to match model pblh definition,
!                            - gnss ro pblh is height above mean sea level.
!                            - model pblh is height above surface height.
!   2023-08     E.-G. Yang - add quality control regarding topography height stdv (meter)
!   2023-08     E.-G. Yang - add quality control regarding coastline
!   2023-10     E.-G. Yang - add quality control regarding data which do not penetrate below 500 m
!   2023-11     E.-G. Yang - add quality control regarding model mixed sfc type
!   2024-02     Y. Zhu     - add observation error model

      use netcdf
      implicit none

!     Declare passed variables
      character(len=*),intent(in):: infile,obstype
      character(20),intent(in):: sis
      integer(i_kind),intent(in):: lunout
      integer(i_kind),intent(inout):: nread,ndata,nodata
      integer(i_kind),dimension(npe),intent(inout):: nobs
      real(r_kind),intent(in):: twindin
!     real(r_kind),dimension(nlat,nlon,nsig),intent(in):: hgtl_full

!     Declare local parameters
      real(r_kind),parameter:: r360 = 360.0_r_kind

      real(r_kind),allocatable,dimension(:,:):: cdata_all

      character(10) date
      logical first,outside,inflate_error,lexist,more_data
      integer(i_kind) iret,im,iy,ihh,idd,i,j,k
      integer(i_kind) ikx,nkx,kx,nreal,ilat,ilon,iout
      integer(i_kind) kk,klon1,klat1,klonp1,klatp1
      integer(i_kind) ntest,nchanl
      integer(i_kind) pblkhqm,maxobs,idomsfc
!     integer(i_kind),dimension(8):: obs_time,anal_time,lobs_time
!     real(r_kind) ltime,cenlon_tmp
      real(r_kind) usage
      real(r_kind) pblkhob,pblkhoe,pblkhelev
      real(r_kind) dlat,dlon,dlat_earth,dlon_earth
      real(r_kind) dlat_earth_deg,dlon_earth_deg
      real(r_kind) cdist,disterr,disterrmax,rlon00,rlat00
      real(r_kind) :: tsavg,ff10,sfcr,zz
!     real(r_kind),dimension(5):: tmp_time
!     real(r_kind) sin2,termg,termr,termrg
!     real(r_kind),dimension(nsig):: zges,hges
      real(r_kind) dx,dy,dx1,dy1,w00,w10,w01,w11

      integer(i_kind) idate5(5),minobs,minan
      real(r_kind) time_correction,timeobs,time,toff,t4dv,zeps

!     real(r_single) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblkh_calipso
      real(r_kind) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblkh_gnssro

      integer(i_kind) :: ncid,ierr,dimid1,dimid2,dimid3,dimid4
      integer(i_kind) :: profiles,levels,ndate,ndatelocal
      integer(i_kind) :: varid1,varid2,varid3,varid4,varid5,varid6
      integer(i_kind) :: varid7,varid8,varid9,varid10,varid11,varid12
      integer(i_kind) :: varid13,varid14,varid15,varid16,varid17,varid18
      integer(i_kind) :: iyear, imonth, iday, ihour, iminute
      integer(i_kind) :: idate
      integer(i_kind) :: ana_time
      integer(i_kind) :: hr, oltime, ltime_mm, ltime_dd, ltime_hh, ltime_min
      real(r_kind) :: hr_min
      real(r_kind), allocatable, dimension(:)    :: lat_pblh, lon_pblh, latstart, lonstart
      real(r_kind), allocatable, dimension(:)    :: pblkh
      real(r_kind), allocatable, dimension(:,:)  :: lat, lon
      real(r_kind), allocatable, dimension(:,:)  :: alt, ref, ref_grad
      integer(i_kind), allocatable, dimension(:) :: gnss_satid, ref_satid, occ_satid,qc_same0
      !gnss_satid: GNSS satellite classification, e.g., 401=GPS, 402=GLONASS
      !ref_satid: GNSS satellite transmitter identifier (1-32)
      !occ_satid: Low Earth Orbit satellite identifier, e.g., COSMIC2=750-755
      character(20), allocatable, dimension(:,:) :: time_start ! UTC time at first level
      character(20), allocatable, dimension(:,:) :: localtime_start ! Local time at first level
      character(20), allocatable, dimension(:,:) :: time_pblh ! UTC time at PBL level

      real(r_kind), allocatable, dimension(:) :: ryear, rmonth, rdate, rhour, rminute
      character(25) :: sfile

      real(r_kind) :: mesh, nobs_within_distance
      real(r_kind) :: hsstdv
      logical :: coast
      integer(i_kind) :: isflg

      ! Observation error model
      real(r_kind) :: cpblh1, cpblh2, cerr1, cerr2

!     Initialize obs err parameters
      cpblh1 = 1500.0_r_kind
      cpblh2 = 4000.0_r_kind
      cerr1 = 250.0_r_kind
      cerr2 = 800.0_r_kind

!     Check if pblkh file exists
      inquire(file=trim(infile),exist=lexist)
      if (.not.lexist) return

!     Read data
      ierr =  NF90_OPEN(trim(infile),0,ncid)
      if (ierr /= nf90_noerr) call handle_err(ierr,"open")

      kx=870  ! GNSSRO PBLH
      nkx=870 ! GNSSRO PBLH

      ierr = NF90_INQ_DIMID(ncid,'profiles',dimid1)
      if (ierr /= nf90_noerr) call handle_err(ierr,"profiles")
      ierr = NF90_INQ_DIMID(ncid,'levels',dimid2)
      if (ierr /= nf90_noerr) call handle_err(ierr,"levels")
      ierr = NF90_INQ_DIMID(ncid,'ndate',dimid3)
      if (ierr /= nf90_noerr) call handle_err(ierr,"ndate")
      ierr = NF90_INQ_DIMID(ncid,'ndatelocal',dimid4)
      if (ierr /= nf90_noerr) call handle_err(ierr,"ndatelocal")

      ierr = nf90_inquire_dimension(ncid, dimid1, len = profiles)
      ierr = nf90_inquire_dimension(ncid, dimid2, len = levels)
      ierr = nf90_inquire_dimension(ncid, dimid3, len = ndate)
      ierr = nf90_inquire_dimension(ncid, dimid4, len = ndatelocal)
      print*, 'read_pblkh,',trim(infile),':: profiles=',profiles,' levels=',levels,' ndate=',ndate,' ndatelocal=',ndatelocal

!     Allocate
      allocate(lat_pblh(profiles), lon_pblh(profiles), latstart(profiles), lonstart(profiles))
      allocate(pblkh(profiles))
      !allocate(lat(profiles,levels), lon(profiles,levels))
      !allocate(alt(profiles,levels), ref(profiles,levels), ref_grad(profiles,levels))
      ! transpose
      allocate(lat(levels,profiles), lon(levels,profiles))
      allocate(alt(levels,profiles), ref(levels,profiles), ref_grad(levels,profiles))
      allocate(gnss_satid(profiles),ref_satid(profiles),occ_satid(profiles),qc_same0(profiles))
      ! not transpose
      allocate(time_start(profiles,ndate), localtime_start(profiles,ndatelocal), time_pblh(profiles,ndate))
      !allocate(time_start(ndate,profiles), localtime_start(ndatelocal,profiles), time_pblh(ndate,profiles))
      allocate(ryear(profiles), rmonth(profiles), rdate(profiles), rhour(profiles), rminute(profiles))

!     Read Variables

!     Analysis Time
      !gnssro_pbl_obs_2015093018.nc4
      !source_file=gnssro_obs_2015093018.nc4
      ierr = NF90_GET_ATT(ncid,nf90_global,'source_file',sfile)
      if (ierr /= nf90_noerr) call handle_err(ierr,"source_file")
      print*, "yeg_read_pblkh_ana_time: sfile=", sfile
      read(sfile(12:21),'(i10)') ana_time
      print*, "yeg_read_pblkh_ana_time: ana_time=", ana_time
      !read(infile(16:25),'(i10)') ana_time
      print*, "yeg_read_pblkh: kx=",kx

!     Latitude at PBLH: degrees
      ierr = NF90_INQ_VARID(ncid,'latitude',varid1)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid1,lat_pblh)
!     Longiude at PBLH: degrees
      ierr = NF90_INQ_VARID(ncid,'longitude',varid2)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid2,lon_pblh)
      ierr = NF90_INQ_VARID(ncid,'time_pblh',varid3)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid3,time_pblh)
!     Split time_pblh
      call split_time(profiles,ndate,time_pblh,ryear,rmonth,rdate,rhour,rminute)

!     PBLH: meters
      ierr = NF90_INQ_VARID(ncid,'PBLH',varid4)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid4,pblkh)
!     Atmospheric refractivity [N]
      ierr = NF90_INQ_VARID(ncid,'ref',varid5)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid5,ref)
!     Geometric altitude [Meters]
      ierr = NF90_INQ_VARID(ncid,'alt',varid6)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid6,alt)
!     Refractivity Gradient [N-unit/km]
      ierr = NF90_INQ_VARID(ncid,'ref_gradient',varid7)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid7,ref_grad)
!     GNSS satellite classification, e.g., 401=GPS, 402=GLONASS
      ierr = NF90_INQ_VARID(ncid,'gnss_satid',varid8)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid8,gnss_satid)
!     GNSS satellite transmitter identifier (1-32)
      ierr = NF90_INQ_VARID(ncid,'ref_satid',varid9)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid9,ref_satid)
!     Low Earth Orbit satellite identifier, e.g., COSMIC2=750-755
      ierr = NF90_INQ_VARID(ncid,'occ_satid',varid10)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid10,occ_satid)


!     Latitude at first level
      ierr = NF90_INQ_VARID(ncid,'latstart',varid11)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid11,latstart)
!     Longitude at first level
      ierr = NF90_INQ_VARID(ncid,'lonstart',varid12)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid12,lonstart)
!     UTC time at first level
      ierr = NF90_INQ_VARID(ncid,'time_start',varid13)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid13,time_start)
!     LOCAL time at first level ("yyyy-mm-dd_hh:mimi")
      ierr = NF90_INQ_VARID(ncid,'localtime_start',varid14)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid14,localtime_start)
!     QC = 0 if same, 1 if different low and high GNSSRO PBLRF
      ierr = NF90_INQ_VARID(ncid,'QC',varid15)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid15,qc_same0)

      ierr = NF90_CLOSE(ncid)
      if (ierr /= nf90_noerr) call handle_err(ierr,"close")

      nreal=18 ! Temporary setting

      maxobs=count(ryear/=0)
      print*, "yeg_read_pblkh: maxobs =",maxobs
      allocate(cdata_all(nreal,maxobs))
      nchanl=0
      nread=0

      first = .true.
      !maxobs=0

      do i = 1, profiles
         print*, " -------------------------"
         print*, "YEGGG:i=",i
!        Time offset
         iyear = ryear(i)
         imonth = rmonth(i)
         iday = rdate(i)
         ihour = rhour(i)
         iminute = rminute(i)

         write (date,'(i4,3i2.2)') iyear,imonth,iday,ihour
         read (date,'( i10)') idate
         print*, 'idate=', idate
         if (first) then
            !call time_4dvar(idate,toff)
            call time_4dvar(ana_time,toff)
            !toff: time since start of 4D-Var window (hours)
            !????????????
            !?toff: here we don't have ana_time. how can we get toff
            print*, 'yeg_pblkh: idate,ana_time,toff=',idate,ana_time,toff
            first=.false.
         end if

         nread=nread+1

!        Is pblkh in convinfo file
         ikx=0
         do j=1,nconvtype
            if(kx == ictype(j)) then
               ikx=j
               exit
            end if
         end do
         print*, "yeg_pblkh: ikx=",ikx
         if(ikx == 0) cycle

         lon_deg = lon_pblh(i)
         lat_deg = lat_pblh(i)
         if(lon_deg>= r360)lon_deg=lon_deg-r360
         if(lon_deg < zero)lon_deg=lon_deg+r360
         dlon_earth_deg=lon_deg
         dlat_earth_deg=lat_deg
         dlon_earth=lon_deg*deg2rad
         dlat_earth=lat_deg*deg2rad
         if(regional)then
            call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)    ! convert to rotated coordinate
            if(diagnostic_reg) then
               call txy2ll(dlon,dlat,rlon00,rlat00)
               ntest=ntest+1
               cdist=sin(dlat_earth)*sin(rlat00)+cos(dlat_earth)*cos(rlat00)* &
                    (sin(dlon_earth)*sin(rlon00)+cos(dlon_earth)*cos(rlon00))
               cdist=max(-one,min(cdist,one))
               disterr=acos(cdist)*rad2deg
               disterrmax=max(disterrmax,disterr)
            end if
            if(outside) cycle   ! check to see if outside regional domain
         else
            dlat = dlat_earth
            dlon = dlon_earth
            call grdcrd1(dlat,rlats,nlat,1)
            call grdcrd1(dlon,rlons,nlon,1)
         endif

         if(offtime_data) then

!          in time correction for observations to account for analysis
!                time being different from obs file time.
           write(date,'( i10)') idate
           read (date,'(i4,3i2)') iy,im,idd,ihh
           idate5(1)=iyear
           idate5(2)=imonth
           idate5(3)=iday
           idate5(4)=ihour
           idate5(5)=0
           call w3fs21(idate5,minobs)    !  obs ref time in minutes relative to historic date
           idate5(1)=iadate(1)
           idate5(2)=iadate(2)
           idate5(3)=iadate(3)
           idate5(4)=iadate(4)
           idate5(5)=0
           call w3fs21(idate5,minan)    !  analysis ref time in minutes relative to historic date

!          Add obs reference time, then subtract analysis time to get obs time relative to analysis

           time_correction=float(minobs-minan)/60._r_kind

         else
           time_correction=zero
         end if

!        Time check
!-----------------------------------------------------------------------
!        timeobs (obs - analysis time [hrs]) should be calculated here.

!        in time correction for observations to account for analysis
!              time being different from obs file time.
         ! obs time

         idate5(1)=iyear
         idate5(2)=imonth
         idate5(3)=iday
         idate5(4)=ihour
         idate5(5)=iminute
         call w3fs21(idate5,minobs)    !  obs ref time in minutes relative to historic date
         idate5(1)=iadate(1)
         idate5(2)=iadate(2)
         idate5(3)=iadate(3)
         idate5(4)=iadate(4)
         idate5(5)=0
         call w3fs21(idate5,minan)    !  analysis ref time in minutes relative to historic date

!        Add obs reference time, then subtract analysis time to get obs time relative to analysis
         timeobs = float(minobs-minan)/60._r_kind
!-----------------------------------------------------------------------
!        timeobs=real(rhour(i), r_double)
         time=timeobs + time_correction
         t4dv=timeobs + toff
         zeps=1.0e-8_r_kind
         if (t4dv<zero  .and.t4dv>      -zeps) t4dv=zero
         if (t4dv>winlen.and.t4dv<winlen+zeps) t4dv=winlen
         t4dv=t4dv + time_correction
         if (l4dvar.or.l4densvar) then
           if (t4dv<zero.OR.t4dv>winlen) cycle
         else
           if((real(abs(time)) > real(ctwind(ikx)) .or. real(abs(time)) > real(twindin))) cycle
         end if
         print*, 'read_pblkh: idate, iadate, timeobs, toff, t4dv=', idate, iadate, timeobs, toff, t4dv


         !------------------------------------------
         ! Local Time (only available for 2015/Aug)
         !------------------------------------------
         ! Converting UTC to Local Time (UTC + hr)
         ! longitude (degree East) -> East: lon>0, West: lon<0
         ! 1) hourly obs (MPLNET and Radar wind profiler)
         ! 2) non-hourly obs (CALIPSO, CATS, raob, GNSSRO)
         if (dlon_earth_deg>180) then
            hr_min = (dlon_earth_deg-360.0)/15.
         else
            hr_min = dlon_earth_deg/15.
         end if

         print*, "yeg_read_pblkh 870 idate=",idate,", dlon_earth_deg=",dlon_earth_deg,", local[hr_min]=",hr_min
         call convert_localtime_min(ltime_mm, ltime_dd, ltime_hh, ltime_min, hr_min, imonth, iday, ihour,iminute)
         oltime = iyear*1000000 + ltime_mm*10000 + ltime_dd*100 + ltime_hh
         print*, "yeg_read_pblkh 870 Local time, oltime=",oltime

!--Outputs

         ndata=ndata+1
         nodata=nodata+1
         iout=ndata

         if(ndata > maxobs) then
            write(6,*)'READ_PBLH:  ***WARNING*** ndata > maxobs for ',obstype
            ndata = maxobs
         end if

!--------------------------------
         !pblkhelev=sfc_elev(i)
         pblkh_gnssro=pblkh(i)

!---------------------------------------------------
!        1) Subtract surface height (zz) from pblkh
!           to match model pblkh (above surface height)

!        Get information from surface file necessary for conventional data here
         call deter_sfc2(dlat_earth,dlon_earth,t4dv,idomsfc,tsavg,ff10,sfcr,zz)
         print*, "read_pblkh L451: pblkh_gnssro=",pblkh_gnssro
         print*, "read_pblkh L452: zz=",zz
         pblkh_gnssro = pblkh_gnssro - zz
         print*, "read_pblkh L454: pblkh_gnssro-zz=",pblkh_gnssro

!        2) Quality mark for pblkh

!        A) Get model topography height stdv at gnss ro pblkh location
         call deter_hsstdv_model(dlat,dlon,hsstdv)
         print*, 'read_pblkh L460: dlat,dlon,hsstdv=',dlat,dlon,hsstdv
         print*, 'read_pblkh L461: topo_stdv=',hsstdv

!        B) Check for coastline and sfc type from model at gnss ro pblkh location
         call deter_sfc_type_coast(dlat_earth,dlon_earth,t4dv,isflg,coast)
         print*, 'read_pblkh coast L466: dlat_earth,dlon_earth=',dlat_earth,dlon_earth
         print*, 'read_pblkh coast L467: isflg,coast=',isflg,coast

!        C) Exclude data that do not penetrate below 500 m
         print*, 'read_pblkh L469: alt(1,i), alt(2,i)=',alt(1,i),alt(2,i)
         print*, 'read_pblkh L470: alt(1,i)-zz, alt(2,i)-zz=',alt(1,i)-zz,alt(2,i)-zz

         pblKhqm=0
         if (pblkh_gnssro .le. 0.0 .or. pblkh_gnssro .gt. 6000.0_r_kind) then
            pblkhqm=15
            !pblkhob=-999.0_r_kind
            print*, 'read_pblkh L472: qc=15, pblkh_gnssro=',pblkh_gnssro
         else if (qc_same0(i).eq.0.0) then ! if qc from obs is 0 -> low and high gnssro pblrf are the same, so exclude
            print*, 'read_pblkh L473: qc=0, meaning low and high pblrf obs are the same'
            pblkhqm=11
         else if (idomsfc>=3) then ! dominate model sfc type >=3 <-- mixed if any surrounding grid has different sfc type (water(0),land(1),ice(2))
            print*, 'read_pblkh L474: idomsfc>=3, idomsfc=',idomsfc
            pblkhqm=20
         else if (hsstdv > 200.0) then ! topography stdv > 200 (steep slope)
            print*, 'read_pblkh L475: qc=9, hsstdv=',hsstdv
            pblkhqm=9
!         else if (coast) then ! if coastline
!            print*, 'read_pblkh L478: qc=10, isflg,coast=',isflg,coast
!            pblkhqm=10
         else if (alt(1,i)-zz>500) then ! if lowest observed altitude - zz [agl] > 500 m
            print*, 'read_pblkh L479: alt(1,i)-zz=',alt(1,i)-zz
            pblkhqm=10
         end if

         if (pblkh_gnssro <= zero) then
            pblkhob=r_missing
         else
            pblkhob=pblkh_gnssro
         end if
!        3) Take log for pblkh
         !if (pblkh_gnssro <= zero) then
         !   pblkhob=r_missing
         !else if (pblkh_gnssro > zero .and. pblkh_gnssro <= one) then
         !   pblkh_gnssro=one
         !   pblkhob=log(pblkh_gnssro)
         !else
         !   pblkhob=log(pblkh_gnssro)
         !end if

!        Set usage variable
         usage = 0.
         if(icuse(ikx) <= 0) usage=150.
         if(pblkhqm == 9 .or. pblkhqm == 10 .or. pblkhqm == 11 .or. pblkhqm == 15 .or. pblkhqm == 20) usage=150.

!        Set inflate_error logical
         inflate_error=.false.
         if (pblkhqm == 3) inflate_error=.true.

!        setup for sort of averaged obs
         !pblkhoe=0.1_r_kind  ! temporarily 100/1400
         !pblkhoe=1.2_r_kind  ! 0.1 -> 1.0 -> 1.5 -> 1.2 (latest)
         !pblkhoe=0.05_r_kind*pblkhob !  0.05 (Jo/n=1000) --> 0.2 (Jo/n=66) --> 0.3 (Jo/n=30)
         if (pblkhob == r_missing) then
            pblkhoe = r_missing
         else
            if(pblkhob <= cpblh1) then
               pblkhoe = cerr1
            else if(pblkhob > cpblh1 .and. pblkhob < cpblh2) then
               pblkhoe = cerr1 + (pblkhob-cpblh1)* &
                           (cerr2-cerr1)/(cpblh2-cpblh1)
            else
               pblkhoe = cerr2
            endif
         endif

         ! Inflate obs error if adjacent obs exists within the distance (mesh), depending on the number of obs.
         mesh=125 ! 250 ! [km]
         call count_dist_obs(profiles,dlat_earth_deg,dlon_earth_deg,lat_pblh,lon_pblh,ryear,mesh,nobs_within_distance)
         ! Multiply sqrt of the number of obs if adjacent obs exists: * sqrt(2) or sqrt(3)
         if (nobs_within_distance > 1.0) pblkhoe=pblkhoe*sqrt(nobs_within_distance)
         print*, "YEG: read_pblkh_ober: pblkhoe=",pblkhoe
         ! if (nkx==120) pblkhoe=0.07_r_kind
         ! if (nkx==181) pblkhoe=0.10_r_kind
         if (inflate_error) pblkhoe=pblkhoe*1.05_r_kind

         cdata_all(1,iout)=pblkhoe                 ! pblkh error (cb)
         cdata_all(2,iout)=dlon                    ! grid relative longitude
         cdata_all(3,iout)=dlat                    ! grid relative latitude
         cdata_all(4,iout)=r_missing !pblkhelev                ! pblkh obs elevation
         cdata_all(5,iout)=pblkhob                 ! pblkh obs
!        cdata_all(6,iout)=r_missing
         cdata_all(6,iout)=kx
         cdata_all(7,iout)=t4dv                    ! time
         cdata_all(8,iout)=ikx                     ! type
         print*, 'YEG_read_pblkh L513: ikx=',ikx
         cdata_all(9,iout)=oltime                  ! local time (YYYYMMDDHH) for all obs types!!!!!!!!!!!!!!!!
         cdata_all(10,iout)=pblkhqm                ! quality mark
         print*, 'YEG_read_pblkh L516: pblkhqm=',pblkhqm,',usage=',usage
         cdata_all(11,iout)=usage                  ! usage parameter
         cdata_all(12,iout)=dlon_earth_deg         ! earth relative longitude (degrees)
         cdata_all(13,iout)=dlat_earth_deg         ! earth relative latitude (degrees)

         !sfc type
         cdata_all(14,iout)=r_missing ! pblriosfc
         cdata_all(15,iout)=idomsfc

         !SAT info
!        GNSS satellite classification, e.g., 401=GPS, 402=GLONASS
         !gnss_satid
         cdata_all(16,iout)=gnss_satid(i) ! gnss sat class
!        GNSS satellite transmitter identifier (1-32)
         !ref_satid
         cdata_all(17,iout)=ref_satid(i) ! gnss sat transmit id
!        Low Earth Orbit satellite identifier, e.g., COSMIC2=750-755
         !occ_satid
         cdata_all(18,iout)=occ_satid(i)

      end do
!   Normal exit

     ilat=3
     ilon=2
     write(6,*) 'YEG3_read_pblkh: ndata, nodata, iout =',ndata,nodata,iout
     write(6,*) 'YEG4_read_pblkh: nobs=',nobs
!   Write observation to scratch file
     call count_obs(ndata,nreal,ilat,ilon,cdata_all,nobs)
     write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
     write(lunout) ((cdata_all(j,i),j=1,nreal),i=1,ndata)
     deallocate(cdata_all)

     deallocate(lat_pblh, lon_pblh, latstart, lonstart, pblkh)
     deallocate(lat, lon, alt, ref, ref_grad)
     deallocate(gnss_satid, ref_satid, occ_satid)
     deallocate(time_start, localtime_start, time_pblh)
     deallocate(ryear, rmonth, rdate, rhour, rminute)

     end subroutine read_pblkh_gnssro


     subroutine split_time(nprofiles,ndate,time,ryear,rmonth,rdate,rhour,rminute)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  split_time      split time variable to year, month, date, hour, min
!
! program history log:
!   2023-04     E.-G. Yang - split YYYY-MM-DDTHH:MM:00Z to have separate variables.

     implicit none

!    Declare passed variables
     integer(i_kind), intent(in) :: nprofiles, ndate
     character(len=*), dimension(nprofiles,ndate), intent(in) :: time
     real(r_kind), dimension(:), allocatable, intent(out) :: ryear,rmonth,rdate,rhour,rminute
     integer :: i

!    Declare local parameters

     allocate(ryear(nprofiles), rmonth(nprofiles), rdate(nprofiles), rhour(nprofiles), rminute(nprofiles))
     do i = 1, nprofiles
        print *, time(i,1)
        read(time(i,1)(1:4),'(f4.0)') ryear(i)
        read(time(i,1)(6:7),'(f2.0)') rmonth(i)
        read(time(i,1)(9:10),'(f2.0)') rdate(i)
        read(time(i,1)(12:13),'(f2.0)') rhour(i)
        read(time(i,1)(15:16),'(f2.0)') rminute(i)
        print *, ryear(i),rmonth(i),rdate(i),rhour(i),rminute(i)
     end do

     end subroutine split_time


     subroutine convert_localtime_hr(ltime_mm, ltime_dd, ltime_hh, interval, ana_mm, ana_dd, ana_hh)
     use kinds, only: i_kind
     integer(i_kind), intent(out) :: ltime_mm, ltime_dd, ltime_hh
     integer(i_kind), intent(in ) :: interval, ana_mm, ana_dd, ana_hh
     ltime_mm = ana_mm
     ltime_dd = ana_dd
     ltime_hh = ana_hh + interval

     if (ltime_hh >= 24) then

        ltime_hh = ltime_hh - 24
        ltime_dd = ana_dd + 1
        if (ana_mm == 8 .and. ltime_dd > 31) then
           ltime_mm = 9
           ltime_dd = 1
        end if
        if (ana_mm == 9 .and. ltime_dd > 30) then
           ltime_mm = 10
           ltime_dd = 1
        end if

     elseif (ltime_hh < 0) then

        ltime_hh = 24 + ltime_hh
        ltime_dd = ana_dd - 1
        if (ana_mm == 9 .and. ana_dd == 1) then
           ltime_mm = 8
           ltime_dd = 31
        end if

     end if
     end subroutine convert_localtime_hr

     subroutine convert_localtime_min(ltime_mm, ltime_dd, ltime_hh, ltime_min, interval_hr_min, ana_mm, ana_dd, ana_hh,ana_min)
     use kinds, only: i_kind,r_kind
     integer(i_kind), intent(out) :: ltime_mm, ltime_dd, ltime_hh, ltime_min
     integer(i_kind), intent(in ) :: ana_mm, ana_dd, ana_hh, ana_min
     real(r_kind), intent(in ) :: interval_hr_min
     real(r_kind) :: ana_hh_real, ltime_hh_real

     ! localtime = UTC + interval_hr_min
     ltime_mm = ana_mm
     ltime_dd = ana_dd
     ! 1) Convert HH:MM to HH.xx
     !ana_hh_real: ana_hh + float(ana_min) to account for minutes
     ana_hh_real = real(ana_hh) + real(ana_min)/60.0
     ! 2) UTC + interval => Local time [hr]
     !ltime_hh = ana_hh + interval_hr_min
     ltime_hh_real = ana_hh_real + interval_hr_min
     ! 3) local time [hr] -> integer HH:MM
     ! 0.9 -> 0, 1.9 -> 1
     ! for negative, -0.9 -> -1, -1.9 -> -2
     ltime_hh = floor(ltime_hh_real) !
     ! 4) local time [minutes] -> integer HH:MM
     ltime_min = nint((ltime_hh_real - ltime_hh)*60.0) ! Round (1.6 -> 2)

     ! only for Aug and Sep
     if (ltime_hh >= 24) then

        ltime_hh = ltime_hh - 24
        ltime_dd = ana_dd + 1
        if (ana_mm == 8 .and. ltime_dd > 31) then
           ltime_mm = 9
           ltime_dd = 1
        end if
        if (ana_mm == 9 .and. ltime_dd > 30) then
           ltime_mm = 10
           ltime_dd = 1
        end if

     elseif (ltime_hh < 0) then

        ltime_hh = 24 + ltime_hh
        ltime_dd = ana_dd - 1
        if (ana_mm == 9 .and. ana_dd == 1) then
           ltime_mm = 8
           ltime_dd = 31
        end if

     end if
     end subroutine convert_localtime_min

     subroutine count_dist_obs(nobs,plat,plon,lat,lon,ryear,mesh,nobs_within_distance)
     use kinds, only: i_kind,r_kind
     use buddycheck_mod, only: gc_dist
     implicit none
     integer(i_kind), intent(in)  :: nobs ! the total number of obs
     real(r_kind), dimension(nobs), intent(in)  :: lat, lon, ryear
     real(r_kind), intent(in)  :: mesh ! (km)
     real(r_kind), intent(in)  :: plat, plon ! (km)
     real(r_kind), intent(out) :: nobs_within_distance ! the number of obs within distance

     integer(i_kind) :: i
     real(r_kind) :: o1_lat, o1_lon, o2_lat, o2_lon, dist

     ! Count the number of obs within distance (mesh)

     o1_lat = plat
     o1_lon = plon

     print*, "plat,plon=",plat,plon
     nobs_within_distance=0
     do i = 1, nobs

        if (ryear(i)/=0.0) then
           o2_lat = lat(i)
           o2_lon = lon(i)

           !function gc_dist(inlat1,inlon1,inlat2,inlon2) [meters]
           !lat1,lon1,lat2,lon2 in degrees
           ! *0.001_r_kind
           dist = gc_dist(o1_lat,o1_lon,o2_lat,o2_lon)*0.001_r_kind

           if (dist <= mesh) then

              nobs_within_distance = nobs_within_distance+1
              if (dist/=0.0) then
                 print*,"yeg_pblkh_distance, dist[km]=",dist,",lat,lon=",lat(i),lon(i)
              end if

           end if
        end if

     end do
     if (nobs_within_distance>1) then
        print*,"yeg_pblkh_distance, nobs_within_distance=",nobs_within_distance
     end if

     end subroutine count_dist_obs


















     subroutine read_pblkh_calipso(nread,ndata,nodata,infile,obstype,lunout,twindin,&
         sis,nobs)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  read_pblkh_calipso     read obs from msgs in NETCDF files
!
! program history log:
!   2022-10     Y. Zhu - adapted from read_pblkh

      use netcdf
      implicit none

!     Declare passed variables
      character(len=*),intent(in):: infile,obstype
      character(20),intent(in):: sis
      integer(i_kind),intent(in):: lunout
      integer(i_kind),intent(inout):: nread,ndata,nodata
      integer(i_kind),dimension(npe),intent(inout):: nobs
      real(r_kind),intent(in):: twindin
!     real(r_kind),dimension(nlat,nlon,nsig),intent(in):: hgtl_full

!     Declare local parameters
      real(r_kind),parameter:: r360 = 360.0_r_kind

      real(r_kind),allocatable,dimension(:,:):: cdata_all

      character(10) date
      logical first,outside,inflate_error,lexist,more_data
      integer(i_kind) iret,im,iy,ihh,idd,i,j,k
      integer(i_kind) ikx,nkx,kx,nreal,ilat,ilon,iout
      integer(i_kind) kk,klon1,klat1,klonp1,klatp1
      integer(i_kind) ntest,nchanl
      integer(i_kind) pblkhqm,maxobs,idomsfc
!     integer(i_kind),dimension(8):: obs_time,anal_time,lobs_time
!     real(r_kind) ltime,cenlon_tmp
      real(r_kind) usage
      real(r_kind) pblkhob,pblkhoe,pblkhelev
      real(r_kind) dlat,dlon,dlat_earth,dlon_earth
      real(r_kind) dlat_earth_deg,dlon_earth_deg
      real(r_kind) cdist,disterr,disterrmax,rlon00,rlat00
      real(r_kind) :: tsavg,ff10,sfcr,zz
!     real(r_kind),dimension(5):: tmp_time
!     real(r_kind) sin2,termg,termr,termrg
!     real(r_kind),dimension(nsig):: zges,hges
      real(r_kind) dx,dy,dx1,dy1,w00,w10,w01,w11

      integer(i_kind) idate5(5),minobs,minan
      real(r_kind) time_correction,timeobs,time,toff,t4dv,zeps

!     real(r_single) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblkh_calipso
      real(r_kind) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblkh_calipso

      integer(i_kind) :: ncid,ierr,dimid1,dimid2,norbits,nheights
      integer(i_kind) :: varid1,varid2,varid3,varid4,varid5,varid6
      integer(i_kind) :: varid7,varid8,varid9,varid10,varid11,varid12
      integer(i_kind) :: iyear, imonth, idate, ihour, iminute
      real(r_kind), allocatable, dimension(:) :: lat, lon, pblkh, sfc_elev, sfc_mask, liadr_data_alt
      real(r_kind), allocatable, dimension(:) :: ryear, rmonth, rdate, rhour, rminute
      real(r_kind), allocatable, dimension(:,:) :: ATB

!     Check if pblkh file exists
      inquire(file=trim(infile),exist=lexist)
      if (.not.lexist) return

!     Read data
      ierr =  NF90_OPEN(trim(infile),0,ncid)
      if (ierr /= nf90_noerr) call handle_err(ierr,"open")

      ierr = NF90_INQ_DIMID(ncid,'norbits',dimid1)
      if (ierr /= nf90_noerr) call handle_err(ierr,"norbits")
      ierr = NF90_INQ_DIMID(ncid,'nheights',dimid2)
      if (ierr /= nf90_noerr) call handle_err(ierr,"nheights")

      ierr = nf90_inquire_dimension(ncid, dimid1, len = norbits)
      ierr = nf90_inquire_dimension(ncid, dimid2, len = nheights)
      print*, 'read_pblkh: norbits=', norbits, ' nheights=', nheights

      allocate(lat(norbits), lon(norbits), pblkh(norbits), sfc_elev(norbits))
      allocate(sfc_mask(norbits), liadr_data_alt(nheights), ATB(norbits, nheights))
      allocate(ryear(norbits), rmonth(norbits), rdate(norbits), rhour(norbits), rminute(norbits))

!     Latitude: degrees
      ierr = NF90_INQ_VARID(ncid,'Latitude',varid1)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid1,lat)
      ierr = NF90_INQ_VARID(ncid,'Longitude',varid2)
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
!     PBL_Height: meters
      ierr = NF90_INQ_VARID(ncid,'PBL_Height',varid8)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid8,pblkh)
!     SurfaceElevation: meters
      ierr = NF90_INQ_VARID(ncid,'SurfaceElevation',varid9)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid9,sfc_elev)
!     0=shallow ocean, 1=land, 2=coastlines, 3=shallow inland water, 4=intermittent water, 
!     5=deep inland water, 6=continental ocean, 7=deep ocean
      ierr = NF90_INQ_VARID(ncid,'Land_Water_Mask',varid10)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid10,sfc_mask)
!     Lidar_Data_Altitude: km
      ierr = NF90_INQ_VARID(ncid,'Lidar_Data_Altitude',varid11)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid11,liadr_data_alt)
!     Total_Attenuated_Backscatter_532: km^-1 sr^-1
      ierr = NF90_INQ_VARID(ncid,'Total_Attenuated_Backscatter_532',varid12)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid12,ATB)

      ierr = NF90_CLOSE(ncid)
      if (ierr /= nf90_noerr) call handle_err(ierr,"close")

      allocate(cdata_all(nreal,maxobs))
      nchanl=0
      nread=0

      first = .true.
      maxobs=0
      do i = 1, norbits

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

         nread=nread+1
!        if(kx == 120) nkx= 120
!        if(kx == 227) nkx= 181
!        if(kx == 888) nkx= 888

!        Is pblkh in convinfo file
         ikx=0
         do j=1,nconvtype
            if(kx == ictype(j)) then
               ikx=j
               exit
            end if
         end do
         if(ikx == 0) cycle

         lon_deg = lon(i)
         lat_deg = lat(i)
         if(lon_deg>= r360)lon_deg=lon_deg-r360
         if(lon_deg < zero)lon_deg=lon_deg+r360
         dlon_earth_deg=lon_deg
         dlat_earth_deg=lat_deg
         dlon_earth=lon_deg*deg2rad
         dlat_earth=lat_deg*deg2rad
         if(regional)then
            call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)    ! convert to rotated coordinate
            if(diagnostic_reg) then
               call txy2ll(dlon,dlat,rlon00,rlat00)
               ntest=ntest+1
               cdist=sin(dlat_earth)*sin(rlat00)+cos(dlat_earth)*cos(rlat00)* &
                    (sin(dlon_earth)*sin(rlon00)+cos(dlon_earth)*cos(rlon00))
               cdist=max(-one,min(cdist,one))
               disterr=acos(cdist)*rad2deg
               disterrmax=max(disterrmax,disterr)
            end if
            if(outside) cycle   ! check to see if outside regional domain
         else
            dlat = dlat_earth
            dlon = dlon_earth
            call grdcrd1(dlat,rlats,nlat,1)
            call grdcrd1(dlon,rlons,nlon,1)
         endif

!        Interpolate guess pressure profile to observation location
!        klon1= int(dlon);  klat1= int(dlat)
!        dx   = dlon-klon1; dy   = dlat-klat1
!        dx1  = one-dx;     dy1  = one-dy
!        w00=dx1*dy1; w10=dx1*dy; w01=dx*dy1; w11=dx*dy

!        klat1=min(max(1,klat1),nlat); klon1=min(max(0,klon1),nlon)
!        if (klon1==0) klon1=nlon
!        klatp1=min(nlat,klat1+1); klonp1=klon1+1
!        if (klonp1==nlon+1) klonp1=1
!        do kk=1,nsig
!           hges(kk)=w00*hgtl_full(klat1 ,klon1 ,kk) +  &
!                    w10*hgtl_full(klatp1,klon1 ,kk) + &
!                    w01*hgtl_full(klat1 ,klonp1,kk) + &
!                    w11*hgtl_full(klatp1,klonp1,kk)
!        end do
!        sin2  = sin(dlat_earth)*sin(dlat_earth)
!        termg = grav_equator * &
!           ((one+somigliana*sin2)/sqrt(one-eccentricity*eccentricity*sin2))
!        termr = semi_major_axis /(one + flattening + grav_ratio -  &
!           two*flattening*sin2)
!        termrg = (termg/grav)*termr
!        do k=1,nsig
!           zges(k) = (termr*hges(k)) / (termrg-hges(k))
!        end do

         if(offtime_data) then

!          in time correction for observations to account for analysis
!                time being different from obs file time.
           write(date,'( i10)') idate
           read (date,'(i4,3i2)') iy,im,idd,ihh
           idate5(1)=iyear
           idate5(2)=imonth
           idate5(3)=idate
           idate5(4)=ihour
           idate5(5)=0
           call w3fs21(idate5,minobs)    !  obs ref time in minutes relative to historic date
           idate5(1)=iadate(1)
           idate5(2)=iadate(2)
           idate5(3)=iadate(3)
           idate5(4)=iadate(4)
           idate5(5)=0
           call w3fs21(idate5,minan)    !  analysis ref time in minutes relative to historic date

!          Add obs reference time, then subtract analysis time to get obs time relative to analysis

           time_correction=float(minobs-minan)/60._r_kind

         else
           time_correction=zero
         end if

!        Time check
         timeobs=real(rhour(i), r_double)
         time=timeobs + time_correction
         t4dv=timeobs + toff
         zeps=1.0e-8_r_kind
         if (t4dv<zero  .and.t4dv>      -zeps) t4dv=zero
         if (t4dv>winlen.and.t4dv<winlen+zeps) t4dv=winlen
         t4dv=t4dv + time_correction
         if (l4dvar.or.l4densvar) then
           if (t4dv<zero.OR.t4dv>winlen) cycle
         else
           if((real(abs(time)) > real(ctwind(ikx)) .or. real(abs(time)) > real(twindin))) cycle 
         end if
         print*, 'read_pblkh: idate, iadate, timeobs, toff, t4dv=', idate, iadate, timeobs, toff, t4dv

         pblkhelev=sfc_elev(i)
         pblkh_calipso=pblkh(i)
         if (pblkh_calipso .lt. 0.0) then
            pblkhqm=15
            pblkhob=-999.0_r_kind
!        0=shallow ocean, 1=land, 2=coastlines, 3=shallow inland water, 4=intermittent water, 
!        5=deep inland water, 6=continental ocean, 7=deep ocean"
         else if (sfc_mask(i) == 2) then
            pblkhqm=9
            pblkhob=-999.0_r_kind
         else
            pblkhqm=0
            if (pblkh_calipso < 0.1_r_kind) pblkh_calipso=0.1_r_kind
            !pblkhob=log(pblkh_calipso)
            pblkhob=pblkh_calipso
         end if

!        if (nkx==131 .or. nkx==133 .or. nkx==135) then
!          anal_time=0
!          obs_time=0
!          tmp_time=zero
!          tmp_time(2)=timeobs
!          anal_time(1)=iadate(1)
!          anal_time(2)=iadate(2)
!          anal_time(3)=iadate(3)
!          anal_time(5)=iadate(4)
!          call w3movdat(tmp_time,anal_time,obs_time) ! observation time

!          lobs_time=0
!          tmp_time=zero
!          cenlon_tmp=hdr(2)
!          if (hdr(2) > 180.0) cenlon_tmp=hdr(2)-360.0_r_kind
!          tmp_time(2)=cenlon_tmp/15.0_r_kind
!          call w3movdat(tmp_time,obs_time,lobs_time) ! local observation time
!          ltime = lobs_time(5)+lobs_time(6)/60.0_r_kind+lobs_time(7)/3600.0_r_kind
!          if ((ltime.gt.21.0) .and. (ltime.lt.5.0)) pblkhqm=3
!        end if

!        Set usage variable
         usage = 0.
         if(icuse(ikx) <= 0) usage=150.
         if(pblkhqm == 15 .or. pblkhqm == 9) usage=150.

!        Set inflate_error logical 
         inflate_error=.false.
         if (pblkhqm == 3) inflate_error=.true.

!--Outputs

         ndata=ndata+1
         nodata=nodata+1
         iout=ndata

         if(ndata > maxobs) then
            write(6,*)'READ_PBLH:  ***WARNING*** ndata > maxobs for ',obstype
            ndata = maxobs
         end if

!        Get information from surface file necessary for conventional data here
         call deter_sfc2(dlat_earth,dlon_earth,t4dv,idomsfc,tsavg,ff10,sfcr,zz)

!        setup for sort of averaged obs 
         pblkhoe=0.1_r_kind  ! temporarily 100/1400
         if (nkx==120) pblkhoe=0.07_r_kind
         if (nkx==181) pblkhoe=0.10_r_kind
         if (inflate_error) pblkhoe=pblkhoe*1.05_r_kind

         cdata_all(1,iout)=pblkhoe                  ! pblkh error (cb)
         cdata_all(2,iout)=dlon                    ! grid relative longitude
         cdata_all(3,iout)=dlat                    ! grid relative latitude
         cdata_all(4,iout)=pblkhelev                ! pblkh obs elevation
         cdata_all(5,iout)=pblkhob                  ! pblkh obs
         cdata_all(6,iout)=r_missing
         cdata_all(7,iout)=t4dv                    ! time
         cdata_all(8,iout)=ikx                     ! type
         cdata_all(9,iout)=r_missing
         cdata_all(10,iout)=pblkhqm                 ! quality mark
         cdata_all(11,iout)=usage                  ! usage parameter
         cdata_all(12,iout)=dlon_earth_deg         ! earth relative longitude (degrees)
         cdata_all(13,iout)=dlat_earth_deg         ! earth relative latitude (degrees)
         cdata_all(14,iout)=sfc_mask(i)
         cdata_all(15:18,iout)=r_missing

      end do
!   Normal exit

!   Write observation to scratch file
     call count_obs(ndata,nreal,ilat,ilon,cdata_all,nobs)
     write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
     write(lunout) ((cdata_all(j,i),j=1,nreal),i=1,ndata)
     deallocate(cdata_all)

     deallocate(lat, lon, pblkh, sfc_elev)
     deallocate(sfc_mask, liadr_data_alt, ATB)
     deallocate(ryear, rmonth, rdate, rhour, rminute)

     end subroutine read_pblkh_calipso

     subroutine read_pblkh_wp_text(nread,ndata,nodata,infile,obstype,lunout,twindin,&
         sis,hgtl_full,nobs)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  read_pblkh     read obs from msgs in wind profiler text data files
!
! program history log:
!   2021 zhu
!
!   input argument list:
!     infile   - unit from which to read BUFR data
!     obstype  - observation type to process
!     lunout   - unit to which to write data for further processing
!     hgtl_full- 3d geopotential height on full domain grid
!
!   output argument list:
!     nread    - number of type "obstype" observations read
!     nodata   - number of individual "obstype" observations read
!     ndata    - number of type "obstype" observations retained for further processing
!     twindin  - input group time window (hours)
!     sis      - satellite/instrument/sensor indicator
!     nobs     - array of observations on each subdomain for each processor
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
      use kinds, only: r_kind,r_double,i_kind,r_single
      use constants, only: zero,one_tenth,one,deg2rad,rad2deg,three,rearth, &
           grav_equator,eccentricity,somigliana,grav_ratio,grav, &
           semi_major_axis,flattening,two
      use gridmod, only: diagnostic_reg,regional,nlon,nlat,nsig, &
           tll2xy,txy2ll,rotate_wind_ll2xy,rotate_wind_xy2ll,&
           rlats,rlons
      use convinfo, only: nconvtype,ctwind, &
           icuse,ictype,ioctype
      use gsi_4dvar, only: l4dvar,l4densvar,time_4dvar,winlen
      use obsmod, only: iadate,offtime_data,bmiss
      use deter_sfc_mod, only: deter_sfc2
      use mpimod, only: npe
      implicit none

!     Declare passed variables
      character(len=*),intent(in):: infile,obstype
      character(20),intent(in):: sis
      integer(i_kind),intent(in):: lunout
      integer(i_kind),intent(inout):: nread,ndata,nodata
      integer(i_kind),dimension(npe),intent(inout):: nobs
      real(r_kind),intent(in):: twindin
      real(r_kind),dimension(nlat,nlon,nsig),intent(in):: hgtl_full

!     Declare local parameters
      real(r_kind),parameter:: r360 = 360.0_r_kind

      integer(i_kind) lunin,idate

      real(r_kind),allocatable,dimension(:,:):: cdata_all

      character(10) date
      logical first,outside,inflate_error,lexist,more_data
      integer(i_kind) ihh,idd,iret,im,iy,istat,i,j,k
      integer(i_kind) ikx,nkx,kx,nreal,ilat,ilon,iout
      integer(i_kind) kk,klon1,klat1,klonp1,klatp1
      integer(i_kind) ntest,nchanl
      integer(i_kind) pblkhqm,maxobs,idomsfc
!     integer(i_kind),dimension(8):: obs_time,anal_time,lobs_time
!     real(r_kind) ltime,cenlon_tmp
      real(r_kind) usage
      real(r_kind) pblkhob,pblkhoe,pblkhelev
      real(r_kind) dlat,dlon,dlat_earth,dlon_earth,stnelev
      real(r_kind) dlat_earth_deg,dlon_earth_deg
      real(r_kind) cdist,disterr,disterrmax,rlon00,rlat00
      real(r_kind) :: tsavg,ff10,sfcr,zz
!     real(r_kind),dimension(5):: tmp_time
      real(r_kind) sin2,termg,termr,termrg
      real(r_kind),dimension(nsig):: zges,hges
      real(r_kind) dx,dy,dx1,dy1,w00,w10,w01,w11

      integer(i_kind) idate5(5),minobs,minan
      real(r_kind) time_correction,timeobs,time,toff,t4dv,zeps

!     real(r_single) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblkh_wp
      real(r_kind) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblkh_wp
      data lunin /50/

!     Initialize variables
!     nreal=14
      nreal=18
      ntest=0
      kx= 888

!     Check if pblkh file exists
      inquire(file=trim(infile),exist=lexist)
      if (.not.lexist) return

      open(lunin,file=trim(infile),form='formatted')
      read(lunin, *) idate

      maxobs=0
      first=.true.
      more_data=.true.
      readfile1: do while(more_data)
         read(lunin,*,iostat=istat) stationid, lat_deg, lon_deg, altitude, localtime, utctime, localday, utcday, pblkh_wp
         if (istat/=0) then
            more_data=.false.
         else
            maxobs=maxobs+1
         end if

!        Time offset
         if (first) then
            call time_4dvar(idate,toff)
            first=.false.
         end if
      end do readfile1

      if (maxobs == 0) return

      allocate(cdata_all(nreal,maxobs))
      nread=0
      nchanl=0
      ilon=2
      ilat=3
      close(lunin)
      open(lunin,file=trim(infile),form='formatted')
      read(lunin, *) idate
      print*, 'read_pblkh idate=',idate,' maxobs=',maxobs
      do i=1,maxobs

         nread=nread+1
         if(kx == 120) nkx= 120
         if(kx == 227) nkx= 181
         if(kx == 888) nkx= 888

!        Is pblkh in convinfo file
         ikx=0
         do j=1,nconvtype
            if(kx == ictype(j)) then
               ikx=j
               exit
            end if
         end do
         if(ikx == 0) cycle

         read(lunin,*) stationid, lat_deg, lon_deg, altitude, localtime, utctime, localday, utcday, pblkh_wp
         print*, 'read_pblkh:', stationid, lat_deg, lon_deg, altitude, localtime, utctime, localday, utcday, pblkh_wp

         if(lon_deg>= r360)lon_deg=lon_deg-r360
         if(lon_deg < zero)lon_deg=lon_deg+r360
         dlon_earth_deg=lon_deg
         dlat_earth_deg=lat_deg
         dlon_earth=lon_deg*deg2rad
         dlat_earth=lat_deg*deg2rad
         if(regional)then
            call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)    ! convert to rotated coordinate
            if(diagnostic_reg) then
               call txy2ll(dlon,dlat,rlon00,rlat00)
               ntest=ntest+1
               cdist=sin(dlat_earth)*sin(rlat00)+cos(dlat_earth)*cos(rlat00)* &
                    (sin(dlon_earth)*sin(rlon00)+cos(dlon_earth)*cos(rlon00))
               cdist=max(-one,min(cdist,one))
               disterr=acos(cdist)*rad2deg
               disterrmax=max(disterrmax,disterr)
            end if
            if(outside) cycle   ! check to see if outside regional domain
         else
            dlat = dlat_earth
            dlon = dlon_earth
            call grdcrd1(dlat,rlats,nlat,1)
            call grdcrd1(dlon,rlons,nlon,1)
         endif

!        Interpolate guess pressure profile to observation location
         klon1= int(dlon);  klat1= int(dlat)
         dx   = dlon-klon1; dy   = dlat-klat1
         dx1  = one-dx;     dy1  = one-dy
         w00=dx1*dy1; w10=dx1*dy; w01=dx*dy1; w11=dx*dy

         klat1=min(max(1,klat1),nlat); klon1=min(max(0,klon1),nlon)
         if (klon1==0) klon1=nlon
         klatp1=min(nlat,klat1+1); klonp1=klon1+1
         if (klonp1==nlon+1) klonp1=1
         do kk=1,nsig
            hges(kk)=w00*hgtl_full(klat1 ,klon1 ,kk) +  &
                     w10*hgtl_full(klatp1,klon1 ,kk) + &
                     w01*hgtl_full(klat1 ,klonp1,kk) + &
                     w11*hgtl_full(klatp1,klonp1,kk)
         end do
         sin2  = sin(dlat_earth)*sin(dlat_earth)
         termg = grav_equator * &
            ((one+somigliana*sin2)/sqrt(one-eccentricity*eccentricity*sin2))
         termr = semi_major_axis /(one + flattening + grav_ratio -  &
            two*flattening*sin2)
         termrg = (termg/grav)*termr
         do k=1,nsig
            zges(k) = (termr*hges(k)) / (termrg-hges(k))
         end do


         if(offtime_data) then

!          in time correction for observations to account for analysis
!                time being different from obs file time.
           write(date,'( i10)') idate
           read (date,'(i4,3i2)') iy,im,idd,ihh
           idate5(1)=iy
           idate5(2)=im
           idate5(3)=idd
           idate5(4)=ihh
           idate5(5)=0
           call w3fs21(idate5,minobs)    !  obs ref time in minutes relative to historic date
           idate5(1)=iadate(1)
           idate5(2)=iadate(2)
           idate5(3)=iadate(3)
           idate5(4)=iadate(4)
           idate5(5)=0
           call w3fs21(idate5,minan)    !  analysis ref time in minutes relative to historic date

!          Add obs reference time, then subtract analysis time to get obs time relative to analysis

           time_correction=float(minobs-minan)/60._r_kind

         else
           time_correction=zero
         end if


!        Time check
         timeobs=real(utctime, r_double)
         time=timeobs + time_correction
         t4dv=timeobs + toff
         zeps=1.0e-8_r_kind
         if (t4dv<zero  .and.t4dv>      -zeps) t4dv=zero
         if (t4dv>winlen.and.t4dv<winlen+zeps) t4dv=winlen
         t4dv=t4dv + time_correction
         if (l4dvar.or.l4densvar) then
           if (t4dv<zero.OR.t4dv>winlen) cycle
         else
           if((real(abs(time)) > real(ctwind(ikx)) .or. real(abs(time)) > real(twindin))) cycle 
         end if
         print*, 'read_pblkh: idate, iadate, timeobs, toff, t4dv=', idate, iadate, timeobs, toff, t4dv

         stnelev=altitude
         pblkhelev=stnelev
         pblkhob=pblkh_wp
         pblkhqm=0
         if (pblkhob .lt. 0.0) pblkhqm=15

!        if (nkx==131 .or. nkx==133 .or. nkx==135) then
!          anal_time=0
!          obs_time=0
!          tmp_time=zero
!          tmp_time(2)=timeobs
!          anal_time(1)=iadate(1)
!          anal_time(2)=iadate(2)
!          anal_time(3)=iadate(3)
!          anal_time(5)=iadate(4)
!          call w3movdat(tmp_time,anal_time,obs_time) ! observation time

!          lobs_time=0
!          tmp_time=zero
!          cenlon_tmp=hdr(2)
!          if (hdr(2) > 180.0) cenlon_tmp=hdr(2)-360.0_r_kind
!          tmp_time(2)=cenlon_tmp/15.0_r_kind
!          call w3movdat(tmp_time,obs_time,lobs_time) ! local observation time
!          ltime = lobs_time(5)+lobs_time(6)/60.0_r_kind+lobs_time(7)/3600.0_r_kind
!          if ((ltime.gt.21.0) .and. (ltime.lt.5.0)) pblkhqm=3
!        end if

!        Set usage variable
         usage = 0.
         if(icuse(ikx) <= 0) usage=150.
         if(pblkhqm == 15 .or. pblkhqm == 9) usage=150.

!        Set inflate_error logical 
         inflate_error=.false.
         if (pblkhqm == 3) inflate_error=.true.

!--Outputs

         ndata=ndata+1
         nodata=nodata+1
         iout=ndata

         if(ndata > maxobs) then
            write(6,*)'READ_PBLH:  ***WARNING*** ndata > maxobs for ',obstype
            ndata = maxobs
         end if


!        Get information from surface file necessary for conventional data here
         call deter_sfc2(dlat_earth,dlon_earth,t4dv,idomsfc,tsavg,ff10,sfcr,zz)

!        setup for sort of averaged obs 
         pblkhoe=100.0_r_kind  ! temporarily
         if (nkx==120) pblkhoe=50.
         if (nkx==181) pblkhoe=100.
         if (inflate_error) pblkhoe=pblkhoe*1.5_r_kind

         cdata_all(1,iout)=pblkhoe                  ! pblkh error (cb)
         cdata_all(2,iout)=dlon                    ! grid relative longitude
         cdata_all(3,iout)=dlat                    ! grid relative latitude
         cdata_all(4,iout)=pblkhelev                ! pblkh obs elevation
         cdata_all(5,iout)=pblkhob                  ! pblkh obs
         cdata_all(6,iout)=stationid               ! station id
         cdata_all(7,iout)=t4dv                    ! time
         cdata_all(8,iout)=ikx                     ! type
         cdata_all(9,iout)=localtime               ! obs local time
         cdata_all(10,iout)=pblkhqm                 ! quality mark
         cdata_all(11,iout)=usage                  ! usage parameter
         cdata_all(12,iout)=dlon_earth_deg         ! earth relative longitude (degrees)
         cdata_all(13,iout)=dlat_earth_deg         ! earth relative latitude (degrees)
         cdata_all(14,iout)=altitude               ! station elevation (m)

      end do
      close(lunin)
!   Normal exit

!   Write observation to scratch file
     call count_obs(ndata,nreal,ilat,ilon,cdata_all,nobs)
     write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
     write(lunout) ((cdata_all(j,i),j=1,nreal),i=1,ndata)
     deallocate(cdata_all)

     if (ndata == 0) then
        close(lunin)
        write(6,*)'READ_pblh:  close(',lunin,')'
     endif

     close(lunin)
     end subroutine read_pblkh_wp_text

  end module read_pblkh
