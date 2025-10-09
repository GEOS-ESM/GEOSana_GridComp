!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: gsimod  ---

!
! !INTERFACE:
!
  module read_pblsld

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
      use deter_sfc_mod, only: deter_sfc2
      use mpimod, only: npe

      implicit none

      private

! !PUBLIC ROUTINES:

      public read_pblsld_calipso_cats_icesat2
!---------------------------------------------------------------------------

   CONTAINS

     subroutine read_pblsld_calipso_cats_icesat2(nread,ndata,nodata,infile,obstype,lunout,twindin,&
         sis,nobs)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  read_pblsld_calipso_cats_icesat2     read obs from msgs in NETCDF files
!
! program history log:
!   2022-10     Y. Zhu - adapted from read_pblh
!   2023-04     E.-G. Yang - added thinning based on distance for calipso and cats
!   2023-08     E.-G. Yang - subtracted terrain height from CALIPSO PBLH (above mean sea level)
!                            in order to match a pblsld model variable (above sfc height).
!                            1) CALIPSO: terrain height from obs file
!   2023-08     E.-G. Yang - added qc for CATS to exclude night-time data over land
!   2023-10     E.-G. Yang - moving average for calipso and cats (within 12 km radius, 25 km thinning)
!                            Before moving average for CALIPSO, sfc height is subtracted.
!   2023-10     E.-G. Yang - changing the way of converting UTC to local time for obs
!                            Except for hourly obs, local time is obtained by adding float(lon/15) to UTC (YYYYMMDD HH:MI)
!                            For hourly obs (MPLNET, Wind Profiler), local time: int(lon/15) + UTC (YYYYMMDD HH:MI)
!   2023-10     E.-G. Yang - removed qc regarding LT for CALIPSO
!   2023-11     E.-G. Yang - added mixed sfc type qc for all obs 
!   2023-12     E.-G. Yang - modified distance w.r.t. analysis grid (50 km) for thinning and moving average of calipso and cats
!   2024-02     Y. Zhu - modelling observation error
!   2024-07     E.-G. Yang - added qc for removing CALIPSO and CATS over land
!   2024-08     E.-G. Yang - added different obs error model for CALIPSO, CATS
!   2024-08     E.-G. Yang - changed obs error (CALIPSO, CATS, GRWP)

      use netcdf
      use buddycheck_mod, only: gc_dist
      implicit none

!     Declare passed variables
      character(len=*),intent(in):: infile,obstype
      character(20),intent(in):: sis
      integer(i_kind),intent(in):: lunout
      integer(i_kind),intent(inout):: nread,ndata,nodata
      integer(i_kind),dimension(npe),intent(inout):: nobs
      real(r_kind),intent(in):: twindin

!     Declare local parameters
      real(r_kind),parameter:: r360 = 360.0_r_kind

      real(r_kind),allocatable,dimension(:,:):: cdata_all

      character(10) date
      logical first,outside,inflate_error,lexist,more_data
      logical thin_calipso, thin_cats, thin_mplnet
      integer(i_kind) iret,im,iy,ihh,idd,i,j,k
      integer(i_kind) ikx,nkx,kx,nreal,ilat,ilon,iout
      integer(i_kind) kk,klon1,klat1,klonp1,klatp1
      integer(i_kind) ntest,nchanl
      integer(i_kind) pblsldqm,maxobs,idomsfc
      integer(i_kind) oltime
      integer(i_kind) hr,ltime_mm,ltime_dd,ltime_hh,ltime_min
      real(r_kind) hr_min
      real(r_kind) usage
      real(r_kind) pblsldob,pblsldoe,pblsldelev,pblsldosfc
      real(r_kind) dlat,dlon,dlat_earth,dlon_earth
      real(r_kind) dlat_earth_deg,dlon_earth_deg
      real(r_kind) cdist,disterr,disterrmax,rlon00,rlat00
      real(r_kind) :: tsavg,ff10,sfcr,zz
      real(r_kind) dx,dy,dx1,dy1,w00,w10,w01,w11

      integer(i_kind) idate5(5),minobs,minan
      real(r_kind) time_correction,timeobs,time,toff,t4dv,zeps

!     real(r_single) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblsld_calipso
      real(r_kind) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblsld_obs_kind

      integer(i_kind) :: ncid,ierr,dimid1,dimid2,dimid3,norbits,nheights
      integer(i_kind) :: noobs, strlen, site_maxstrlen, instrument_maxstrlen
      integer(i_kind) :: varid0,varid1,varid2,varid3,varid4,varid5,varid6
      integer(i_kind) :: varid7,varid8,varid9,varid10,varid11,varid12,varid13
      integer(i_kind) :: varid14, varid15
      integer(i_kind) :: ioyr, iomo, iody, iohr, iomi, iodate
      integer(i_kind), dimension(1) :: ana_time
      real(r_kind), allocatable, dimension(:) :: lat, lon, pblsld
      real(r_kind), allocatable, dimension(:) :: ryear, rmonth, rdate, rhour, rminute
      ! calipso
      real(r_kind), allocatable, dimension(:) :: sfc_elev, sfc_mask_calipso, liadr_data_alt
      real(r_kind), allocatable, dimension(:,:) :: ATB
      ! cats
      real(r_kind), allocatable, dimension(:) :: loc_hour, qc_flag, day_night_flag, sfc_mask_cats
      ! Thinning
      real(r_kind) :: o1_lat, o1_lon, o2_lat, o2_lon, dist, thindist
      integer(i_kind) :: nthinobs
      ! Thinning for CALIPSO
      real(r_kind), allocatable, dimension(:) :: lat_thin_final, lon_thin_final, pblsld_mvavg_thin_final
      ! Thinning for cats
      real(r_kind) :: mvavg_rad
      real(r_kind), allocatable, dimension(:) :: lat_qc_thin_final, lon_qc_thin_final, pblsld_qc_mvavg_thin_final
      ! Observation error model
      real(r_kind) :: cpblh1, cpblh2, cerr1a, cerr2a
      real(r_kind) :: cerr1b, cerr2b
      real(r_kind) :: cerr1c, cerr2c
      real(r_kind) :: cerr1d, cerr2d

      logical :: selected

!     Initialize obs err parameters
      cpblh1 = 2000.0_r_kind
      cpblh2 = 4000.0_r_kind
      cerr1c = 350.0_r_kind ! calipso
      cerr1d = 400.0_r_kind ! cats
      cerr2c = 750.0_r_kind
      cerr2d = 850.0_r_kind

!     Check if pblsld file exists
      inquire(file=trim(infile),exist=lexist)
      if (.not.lexist) return

!     Read data
      ierr =  NF90_OPEN(trim(infile),0,ncid)
      if (ierr /= nf90_noerr) call handle_err(ierr,"open")

!     Read dimensions depending on observation types
      if (index(trim(infile),'calipso') > 0) then

         kx = 890 ! calipso
         nkx = 890 ! calipso
         ierr = NF90_INQ_DIMID(ncid,'norbits',dimid1)
         if (ierr /= nf90_noerr) call handle_err(ierr,"norbits")
         ierr = NF90_INQ_DIMID(ncid,'nheights',dimid2)
         if (ierr /= nf90_noerr) call handle_err(ierr,"nheights")

         ierr = nf90_inquire_dimension(ncid, dimid1, len = norbits)
         ierr = nf90_inquire_dimension(ncid, dimid2, len = nheights)
         print*, 'read_pblsld,',trim(infile),': norbits=', norbits, ' nheights=', nheights

      else if (index(trim(infile),'cats') > 0) then

         kx = 891 ! cats
         nkx = 891 ! cats
         ierr = NF90_INQ_DIMID(ncid,'nobs',dimid1)
         if (ierr /= nf90_noerr) call handle_err(ierr,"nobs")

         ierr = nf90_inquire_dimension(ncid, dimid1, len = noobs)
         print*, 'read_pblsld,',infile,': noobs=', noobs

      end if

      select case(kx) ! Temporary setting
         case (890) ! calipso
           maxobs=norbits 
         case (891) ! cats
            maxobs=noobs
      end select

!     Allocate
      allocate(lat(maxobs), lon(maxobs), pblsld(maxobs))
      allocate(ryear(maxobs), rmonth(maxobs), rdate(maxobs), rhour(maxobs), rminute(maxobs))
      select case(kx)
         case (890) ! calipso
            allocate(sfc_elev(norbits), sfc_mask_calipso(norbits), liadr_data_alt(nheights), ATB(norbits,nheights))
         case (891) ! cats
            allocate(loc_hour(noobs), qc_flag(noobs), day_night_flag(noobs), sfc_mask_cats(noobs))
      end select

!     Read variables

!     Analysis Time
      ierr = NF90_INQ_VARID(ncid,'Ana_Time',varid0)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid0,ana_time)
      print*, "YEGREAD kx=",kx,",ana_time=", ana_time
!     Latitude: degrees
      ierr = NF90_INQ_VARID(ncid,'lat',varid1)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid1,lat)
!     Longitude: degrees
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
!     PBL_Height: meters
      select case(kx)
         case (890) ! calipso
            ierr = NF90_INQ_VARID(ncid,'PBL_Height',varid8)
         case (891) ! cats
            ierr = NF90_INQ_VARID(ncid,'PBL_Height',varid8)
      end select

      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid8,pblsld)

!     Other variables
      select case(kx)

         case (890) ! calipso

            !SurfaceElevation: meters
            ierr = NF90_INQ_VARID(ncid,'SurfaceElevation',varid9)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid9,sfc_elev)
            !0=shallow ocean, 1=land, 2=coastlines, 3=shallow inland water, 4=intermittent water, 
            !5=deep inland water, 6=continental ocean, 7=deep ocean
            ierr = NF90_INQ_VARID(ncid,'Land_Water_Mask',varid10)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid10,sfc_mask_calipso)
            !Lidar_Data_Altitude: km
            ierr = NF90_INQ_VARID(ncid,'Lidar_Data_Altitude',varid11)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid11,liadr_data_alt)
            !Total_Attenuated_Backscatter_532: km^-1 sr^-1
            ierr = NF90_INQ_VARID(ncid,'Total_Attenuated_Backscatter_532',varid12)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid12,ATB)

         case (891) ! cats

            !Local_Hour: hour
            ierr = NF90_INQ_VARID(ncid,'Local_Hour',varid9)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid9,loc_hour)
            ierr = NF90_INQ_VARID(ncid,'QC_Flag',varid10)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid10,qc_flag)
            ierr = NF90_INQ_VARID(ncid,'Day_Night_Flag',varid11)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid11,day_night_flag)
            ierr = NF90_INQ_VARID(ncid,'Land_Water_Mask',varid12)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid12,sfc_mask_cats)

      end select

      ierr = NF90_CLOSE(ncid)
      if (ierr /= nf90_noerr) call handle_err(ierr,"close")

      nreal=17 ! 16 ! 15 ! Temporary setting

      !-------------------------------
      ! Spatial thinning for CALIPSO and CATS
      !-------------------------------
      thin_calipso=.true.
      thin_cats=.true.
      !-------------------------------
      ! COUNT for thinning: CALIPSO (subtracting sfc_elev and moving average)
      !-------------------------------
      if (thin_calipso .and. kx==890) then

         !count the number of obs after thinning 
         thindist=50.0_r_kind ! [km] GSI Grid, analysis grid (lower resolution than model grid)
         mvavg_rad=25.0_r_kind !12! [km] moving average radius 25 km
         call count_elev_mvavg_thin_obs(maxobs,lat,lon,sfc_elev,pblsld,thindist,mvavg_rad,nthinobs,lat_thin_final,lon_thin_final,pblsld_mvavg_thin_final)

         allocate(cdata_all(nreal,nthinobs))
         print*, "yeg_thinning_calipso: maxobs=",maxobs,",nthinobs=",nthinobs

      !-------------------------------
      ! COUNT for thinning: CATS (considering qc and moving average)
      !-------------------------------
      elseif (thin_cats .and. kx==891) then

         !count the number of obs after thinning 
         thindist=50.0_r_kind ! [km] GSI Grid
         mvavg_rad=25.0_r_kind ! [km] moving average radius 25 km
         !Thinning based on distance as well as qc
         call count_qc_mvavg_thin_obs(maxobs,lat,lon,qc_flag,pblsld,thindist,mvavg_rad,nthinobs,lat_qc_thin_final,lon_qc_thin_final,pblsld_qc_mvavg_thin_final)

         allocate(cdata_all(nreal,nthinobs))
         print*, "yeg_thinning_cats: maxobs=",maxobs,",nthinobs=",nthinobs

!      else
!         allocate(cdata_all(nreal,maxobs))

      end if
      nchanl=0
      nread=0

      first = .true.
      !maxobs=0
      !do i = 1, norbits
      do i = 1, maxobs

         print*, " -------------------------"
!        Time offset         
!        Obs Time
         ioyr = ryear(i)
         iomo = rmonth(i)
         iody = rdate(i)
         iohr = rhour(i)
         iomi = rminute(i)

         if (ioyr==r_missing) then
            write (date,'(i4,3i2.2)') 0,0,0,0
         else
            write (date,'(i4,3i2.2)') ioyr,iomo,iody,iohr
         end if
         read (date,'( i10)') iodate

         print*, 'i=',i,', kx=',kx,',iodate=', iodate

         if (first) then
            !call time_4dvar(idate,toff) 
            call time_4dvar(ana_time(1),toff) 
            ! toff: time since start of 4D-Var window (hours)
            print*, 'YEG12: iodate,toff=',iodate,toff
            first=.false.
         end if

         nread=nread+1
         !kx=890
         if(kx == 120) nkx= 120
         if(kx == 890) nkx= 890
         if(kx == 891) nkx= 891
         if(kx == 892) nkx= 892
         if(kx == 893) nkx= 893

         !----------------------------------------------
         ! thinning calipso depending on distance
         if (thin_calipso .and. kx==890) then

            selected=.False.
            findloopp: do j = 1, nthinobs

               if (lat(i)==lat_thin_final(j) .and. lon(i)==lon_thin_final(j)) then
                  selected=.True.
                  pblsld(i)=pblsld_mvavg_thin_final(j)
                  exit findloopp
               end if

            end do findloopp

            if (.not.selected) cycle

         end if

         !----------------------------------------------
         ! thinning cats depending on distance as well as qc
         if (thin_cats .and. kx==891) then

            selected=.False.
            findloop: do j = 1, nthinobs

               if (lat(i)==lat_qc_thin_final(j) .and. lon(i)==lon_qc_thin_final(j)) then
                  selected=.True.
                  pblsld(i)=pblsld_qc_mvavg_thin_final(j)
                  exit findloop
               end if

            end do findloop

            if (.not.selected) cycle

         end if


!        Is pblsld in convinfo file
         ikx=0
         loop_convinfo_pblsld: do j=1,nconvtype
           if (trim(ioctype(j)) /= trim(obstype))cycle loop_convinfo_pblsld
           if (nkx == ictype(j)) then
              ikx=j
              exit
           end if
         end do loop_convinfo_pblsld

         print*, "YEG_read_pblsld L819:nkx=",nkx,", ikx=",ikx
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

!           in time correction for observations to account for analysis
!                time being different from obs file time.
            write(date,'( i10)') iodate
            read (date,'(i4,3i2)') iy,im,idd,ihh
            !idate5(1)=iyear
            !idate5(2)=imonth
            !idate5(3)=idate
            !idate5(4)=ihour
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

!           Add obs reference time, then subtract analysis time to get obs time relative to analysis

            time_correction=float(minobs-minan)/60._r_kind

         else
            time_correction=zero
         end if

!-----------------------------------------------------------------------
!        Time check
!-----------------------------------------------------------------------
!        timeobs (obs - analysis time [hrs]) should be calculated here. 

!        in time correction for observations to account for analysis
!              time being different from obs file time.
         ! obs time
         idate5(1)=ioyr
         idate5(2)=iomo
         idate5(3)=iody
         idate5(4)=iohr
         idate5(5)=iomi
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
         print*, 'read_pblsld: iodate, iadate, timeobs, toff, t4dv=', iodate, iadate, timeobs, toff, t4dv
         print*, 'YEG GRID: t4dv =', t4dv

         !------------------------------------------
         ! Local Time (only available for 2015/Aug)
         !------------------------------------------
         ! Converting UTC to Local Time (UTC + hr)
         ! longitude (degree East) -> East: lon>0, West: lon<0
         ! 1) hourly obs (MPLNET and Radar wind profiler)
         ! 2) non-hourly obs (CALIPSO, CATS, raob)

         select case(kx)

            case(890) ! CALIPSO
               ! Converting UTC to Local Time (UTC + hr_min)
               ! longitude (degree East) -> East: lon>0, West: lon<0
               if (lon(i)>180) then
                  hr_min = (lon(i)-360.0)/15.
               else
                  hr_min = lon(i)/15.
               end if

               print*, "yeg_read_pblsld 890 L964 iodate=",iodate,", lon(i)=",lon(i),", hr_min=",hr_min
               call convert_localtime_min(ltime_mm, ltime_dd, ltime_hh, ltime_min, hr_min, iomo, iody, iohr,iomi)
               oltime = ioyr*1000000 + ltime_mm*10000 + ltime_dd*100 + ltime_hh
               print*, "yeg_read_pblsld 890 L966 Local time, oltime=",oltime

            case(891) ! CATS
               !oltime = r_missing ! using localtime from cats or newly calculate ??????????? 
               ! Converting UTC to Local Time (UTC + hr_min)
               ! longitude (degree East) -> East: lon>0, West: lon<0
               if (lon(i)>180) then
                  hr_min = (lon(i)-360.0)/15.
               else
                  hr_min = lon(i)/15.
               end if

               print*, "yeg_read_pblsld 891 iodate=",iodate,", lon(i)=",lon(i),", hr_min=",hr_min
               call convert_localtime_min(ltime_mm, ltime_dd, ltime_hh, ltime_min, hr_min, iomo, iody, iohr,iomi)
               oltime = ioyr*1000000 + ltime_mm*10000 + ltime_dd*100 + ltime_hh
               print*, "yeg_read_pblsld 891 Local time, oltime=",oltime

         end select

!        Get information from surface file necessary for conventional data here
         call deter_sfc2(dlat_earth,dlon_earth,t4dv,idomsfc,tsavg,ff10,sfcr,zz)
         !-----------------------
         ! Sfc Elevation from obs (CALIPSO, MPLNET, RAOB) or model at obs location (CATS)
         !-----------------------
         !pblsldelev=sfc_elev(i)
         select case(kx)
            case(890)
               pblsldelev=sfc_elev(i)
            case(891) ! no info from obs file (CATS)
               !pblsldelev=r_missing
               pblsldelev=zz ! surface height (zz) from model at obs location
         end select

         !--------------------------------------------------
         ! Subtract terrain height from CALIPSO
         ! in order to match model pblsld (above sfc height)
         ! CALIPSO = height above mean sea level
         !--------------------------------------------------
         ! Subtract sfc height from CALIPSO (will be done in moving average and thinning)
         !if (kx==890) then 

         !   pblsld_obs_kind = pblsld(i) - pblsldelev
         !   print*, "YEG: kx=",kx,",pblsld(i)=",pblsld(i)
         !   print*, "YEG: pblsldelev=",pblsldelev
         !   print*, "YEG: pblsld-elev=",pblsld_obs_kind

         !else

         !   pblsld_obs_kind=pblsld(i)

         !end if
         pblsld_obs_kind=pblsld(i)
         !------------------------------
         ! Obs sfc type
         ! Make the consistent land_water_mask for all obs types
         ! pblsldosfc: 1=Land, 2=water, 3=ice and snow, 4=mixed (coastlines)
         !------------------------------

         select case(kx)

            case(890) ! CALIPSO
            !0=shallow ocean, 1=land, 2=coastlines, 3=shallow inland water, 4=intermittent water, 
            !5=deep inland water, 6=continental ocean, 7=deep ocean"

               if (sfc_mask_calipso(i) == 1) then ! land
                  pblsldosfc=1
               else if (sfc_mask_calipso(i) == 2) then ! coastlines
                  pblsldosfc=4
               else ! water 
                  pblsldosfc=2
               end if

            case(891) ! CATS
            !; 17=water; 1=evergreen needleleaf forest, 2=evergreen broadleaf forest,
            !; 3=deciduous needleleaf forest, 4=deciduous broadleaf forest, 5=mixed forest,
            !; 6=closed shrublands, 7=open shrublands, 8=woody savannas, 9=savannas,
            !; 10=grasslands, 11=permanent wetlands, 12=croplands, 13=urban,
            !; 14=Cropland/Natural Vegetation Mosaic, 15=Permanent Snow and Ice,
            !; 16=Barren or Sparsely Vegetated
               if (sfc_mask_cats(i) == 17) then ! water
                  pblsldosfc=2
               else if (sfc_mask_cats(i) == 15) then ! permanent snow and ice
                  pblsldosfc=3
               else ! land (for CATS, there is no coastline category for sfc type)
                  pblsldosfc=1
               end if

         end select

         !----------------------
         ! QC
         !----------------------
         ! Initialize pblsldqm
         pblsldqm=0

         !------------
         ! CATS
         !------------
         !For CATS, recommended to use Qflag > 1 (do not use flag < 2).
         if (kx == 891) then
            if (qc_flag(i) < 2) then
               pblsldqm=9
            ! QC for night-time over land
            else if (mod(oltime,100) < 10 .or. mod(oltime,100) > 17) then ! night-time
               ! land
               if (pblsldosfc /= 2) then ! pblsldosfc = 2 (sfc_mask_cats(i)=17) : water
                                        ! pblsldosfc /= 2 : land
                  pblsldqm=10
                  print*, "yeg_cats_qc: sfc=",pblsldosfc
                  print*, "yeg_cats_qc: mod(oltime,100)=",mod(oltime,100)
               end if

            ! QC over land
            ! pblsldosfc: 1=Land, 2=water, 3=ice and snow, 4=mixed (coastlines)
            !else if (pblsldosfc /= 2) then ! pblsldosfc = 2 (sfc_mask_cats(i)=17) : water
            !                            ! pblsldosfc /= 2 : land, ice and snow
            !   pblsldqm=10
            !   print*, "yeg_cats_qc land: 1: land, 3:ice,snow, pblsldosfc=",pblsldosfc
            end if

         end if

         !----------
         ! CALIPSO
         !----------
         if (kx == 890) then ! calipso
            if (pblsldosfc == 4) then ! coastline from obs data
               pblsldqm=9

            ! QC over land
            ! pblsldosfc: 1=Land, 2=water, 3=ice and snow, 4=mixed (coastlines)
            !else if (pblsldosfc /= 2) then ! pblsldosfc = 2 : water
            !                            ! pblsldosfc /= 2 : land (no 3:ice_and_snow for calipso)
            !   pblsldqm=10
            !   print*, "yeg_calipso_qc land: 1: land, pblsldosfc=",pblsldosfc

            !QC for CALIPSO: Latitude 
            !only consider -60 < lat < 60 for calipso
            else if (lat(i)>60.0_r_kind .or. lat(i) < (-1.0)*60.0_r_kind) then
               pblsldqm=11
               print*, "yeg_calipso_qc: lat(i)=",lat(i)
            !We are not considering localtime for CALIPSO any more.
            end if
         end if

         !----------
         ! ALL OBS
         !----------
         !?????? Discuss 6000 is proper or not.
         if (pblsld_obs_kind .le. 0.0 .or. pblsld_obs_kind .gt. 6000_r_kind) then
            pblsldqm=15
            print*, 'pblsldqm=15, pblh<0 or pblh>6000, pblsld_obs_kind=',pblsld_obs_kind
         else if (idomsfc>=3) then ! dominate model sfc type >=3 <-- mixed if any surrounding grid has different sfc type (water(0),land(1),ice(2))
            pblsldqm=20
            print*, 'pblsldqm=20, idomsfc>=3 (mixed sfc type), idomsfc=',idomsfc
         end if

         if (pblsld_obs_kind <= zero) then
            pblsldob=r_missing
         else
            pblsldob=pblsld_obs_kind
         end if

!        Set usage variable
         usage = 0.
         print*, 'YEG18: icuse(ikx)=',icuse(ikx), ',ikx=',ikx
         if(icuse(ikx) <= 0) usage=150.
         ! qm=15 for all obs (when pblh < 0 or pblh > 6000)
         ! qm=20 for all obs (when model sfc type = mixed)
         ! qm=9 (coastlines),11 (high lat) for CALIPSO
         ! qm=9 (qc flag),10 (night-time land) for CATS
         ! qm=9 (night-time) for MPLNET
         ! qm=9 (qc flag),10 (night-tiime) for wind profiler
         if(pblsldqm == 9 .or. pblsldqm == 10 .or. pblsldqm == 11 .or. pblsldqm == 15 .or. pblsldqm == 20) usage=150.
         print*, 'YEG19: usage=',usage

!        Set inflate_error logical 
         inflate_error=.false.
         if (pblsldqm == 3) inflate_error=.true.

!--Outputs

         ndata=ndata+1
         nodata=nodata+1
         iout=ndata

         write(6,*) 'YEG1: maxobs, ndata, nodata, iout =',maxobs,ndata,nodata,iout

         if(ndata > maxobs) then
            write(6,*)'READ_PBLRI:  ***WARNING*** ndata > maxobs for ',obstype
            ndata = maxobs
         end if

         write(6,*) 'YEG2: maxobs, ndata, nodata, iout =',maxobs, ndata,nodata,iout

!        setup for sort of averaged obs 
         !---------------
         ! Obs Error
         !---------------
         !pblsldoe=0.1_r_kind  ! temporarily 100/1400
         !pblsldoe=1.31_r_kind  ! temporarily 100/1400
         !pblsldoe=1.5_r_kind ! 2.0_r_kind  ! Mar/10/2023
         !pblsldoe=0.075_r_kind ! 2.0_r_kind  ! Mar/10/2023

         !pblsldoe=0.1_r_kind*pblsldob ! 1.5_r_kind was used
         !pblsldoe=0.05_r_kind*pblsldob ! 1.5_r_kind was used
         if (pblsldob == r_missing) then
            pblsldoe = r_missing
         else
            if (nkx==890) then ! CALIPSO
               if(pblsldob <= cpblh1) then
                  pblsldoe = cerr1c
               else if(pblsldob > cpblh1 .and. pblsldob < cpblh2) then
                  pblsldoe = cerr1c + (pblsldob-cpblh1)* &
                           (cerr2c-cerr1c)/(cpblh2-cpblh1)
               else 
                  pblsldoe = cerr2c
               endif         
            else if (nkx==891) then ! CATS
               if(pblsldob <= cpblh1) then
                  pblsldoe = cerr1d
               else if(pblsldob > cpblh1 .and. pblsldob < cpblh2) then
                  pblsldoe = cerr1d + (pblsldob-cpblh1)* &
                           (cerr2d-cerr1d)/(cpblh2-cpblh1)
               else 
                  pblsldoe = cerr2d
               endif         
            endif         
         endif         
    
         if (inflate_error) pblsldoe=pblsldoe*1.05_r_kind

         cdata_all(1,iout)=pblsldoe                 ! pblsld error (cb)
         cdata_all(2,iout)=dlon                    ! grid relative longitude
         cdata_all(3,iout)=dlat                    ! grid relative latitude
         cdata_all(4,iout)=pblsldelev               ! pblsld obs elevation
         cdata_all(5,iout)=pblsldob                 ! pblsld obs
         cdata_all(6,iout)=r_missing               ! index of station id
         cdata_all(7,iout)=t4dv                    ! time
         cdata_all(8,iout)=ikx                     ! type
         print*, 'YEG_read_pblsld L1087: ikx=',ikx
         cdata_all(9,iout)=oltime                  ! local time (YYYYMMDDHH) for all obs types!!!!!!!!!!!!!!!!
         cdata_all(10,iout)=pblsldqm                ! quality mark
         print*, 'YEG_read_pblsld L1090: pblsldqm=',pblsldqm,',usage=',usage
         cdata_all(11,iout)=usage                  ! usage parameter
         cdata_all(12,iout)=dlon_earth_deg         ! earth relative longitude (degrees)
         cdata_all(13,iout)=dlat_earth_deg         ! earth relative latitude (degrees)

         !sfc type
         cdata_all(14,iout)=pblsldosfc
         cdata_all(15,iout)=idomsfc
         cdata_all(16,iout)=r_missing
         cdata_all(17,iout)=kx                     ! 120 or 890-893

      end do
!     Normal exit

      ilat=3
      ilon=2
      write(6,*) 'YEG3: ndata, nodata, iout =',ndata,nodata,iout

!     Write observation to scratch file
      call count_obs(ndata,nreal,ilat,ilon,cdata_all,nobs)
      write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
      write(lunout) ((cdata_all(j,i),j=1,nreal),i=1,ndata)
      deallocate(cdata_all)
 
      deallocate(lat, lon, pblsld)
      deallocate(ryear, rmonth, rdate, rhour, rminute)
      select case(kx)
         case(890) ! calipso
            deallocate(sfc_elev, sfc_mask_calipso, liadr_data_alt, ATB)
         case(891) ! cats
            deallocate(loc_hour, qc_flag, day_night_flag, sfc_mask_cats)
      end select

      end subroutine read_pblsld_calipso_cats_icesat2

     subroutine convert_localtime_hr(ltime_mm, ltime_dd, ltime_hh, interval, ana_mm, ana_dd, ana_hh)
     use kinds, only: i_kind
     integer(i_kind), intent(out) :: ltime_mm, ltime_dd, ltime_hh
     integer(i_kind), intent(in ) :: interval, ana_mm, ana_dd, ana_hh
     ltime_mm = ana_mm
     ltime_dd = ana_dd
     ltime_hh = ana_hh + interval

     ! localtime = UTC + hr
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

     subroutine count_elev_mvavg_thin_obs(nobs,lat,lon,sfc_elev,pblsld,mesh,mvavg_rad,nthinobs,lat_thin,lon_thin,pblsld_mvavg_thin)
     use kinds, only: i_kind,r_kind
     use buddycheck_mod, only: gc_dist
     implicit none
     integer(i_kind), intent(in)  :: nobs ! the number of obs before thinning
     real(r_kind), dimension(nobs), intent(in)  :: lat, lon
     real(r_kind), dimension(nobs), intent(in)  :: sfc_elev,pblsld
     real(r_kind), intent(in)  :: mesh ! thinning size (km)
     real(r_kind), intent(in)  :: mvavg_rad! moving average radius
     integer(i_kind), intent(out) :: nthinobs ! the number of obs after thinning
     real(r_kind), dimension(:), allocatable, intent(out) :: lat_thin, lon_thin
     real(r_kind), dimension(:), allocatable, intent(out)  :: pblsld_mvavg_thin

     integer(i_kind) :: i,j,k
     integer(i_kind) :: begin,ending
     integer(i_kind) :: no_obs_mvavg
     real(r_kind) :: o1_lat, o1_lon, o2_lat, o2_lon, dist
     real(r_kind) :: tobs_mvavg, sum_pblsld

     ! FOR CALIPSO 
     ! Count the number of obs after
     ! 1) Thinning depending on distance (mesh)
     ! 2) Substracting sfc_elev from pblsld CALIPSO
     ! 3) calculating moving average within radius (mvavg_rad) for thinned obs
 
     ! 1) Thinning depending on distance (mesh)
     ! Count the number of obs after thinning depending on distance (mesh)
     do i = 1, nobs
        if (i == 1) then
           o1_lat = lat(i)
           o1_lon = lon(i)
           nthinobs=1
        else
           o2_lat = lat(i)
           o2_lon = lon(i)
           !function gc_dist(inlat1,inlon1,inlat2,inlon2) [meters]
           !lat1,lon1,lat2,lon2 in degrees
           ! *0.001_r_kind
           dist = gc_dist(o1_lat,o1_lon,o2_lat,o2_lon)*0.001_r_kind
           !print*,"yeg_thinning, i=",i,",dist[km]=",dist

           if (dist >= mesh) then
              nthinobs = nthinobs+1
              !print*,"yeg_thinning, dist[km]=",dist,",nthinobs=",nthinobs
              o1_lat = o2_lat
              o1_lon = o2_lon
           end if
        end if
     end do

     allocate(lat_thin(nthinobs),lon_thin(nthinobs),pblsld_mvavg_thin(nthinobs))

     ! 2) Save thinned obs and substract sfc_elev from pblsld
     ! 3) calculating moving average within radius (mvavg_rad) for thinned obs
     ! Save thinned QCed obs (lat,lon,pblsld)
     no_obs_mvavg=80 !40 !80 !40 ! moving average --> +- 80 obs along the track (25 km, 330m between obs)
     tobs_mvavg=0 ! total number of obs for moving avg
     sum_pblsld=0

     do i = 1, nobs
        if (i == 1) then
           nthinobs=1
           o1_lat = lat(i)
           o1_lon = lon(i)

           lat_thin(nthinobs)=lat(i)
           lon_thin(nthinobs)=lon(i)
           sum_pblsld = sum_pblsld + pblsld(i) - sfc_elev(i) ! subtracting sfc height
           tobs_mvavg = tobs_mvavg + 1
        else
           !--------------------------
           ! Decide thinning or not
           o2_lat = lat(i)
           o2_lon = lon(i)
           !function gc_dist(inlat1,inlon1,inlat2,inlon2) [meters]
           !lat1,lon1,lat2,lon2 in degrees
           ! *0.001_r_kind
           dist = gc_dist(o1_lat,o1_lon,o2_lat,o2_lon)*0.001_r_kind
           !print*,"yeg_thinning, i=",i,",dist[km]=",dist

           if (dist >= mesh) then

              tobs_mvavg=0 ! total number of obs for moving avg
              sum_pblsld=0
              nthinobs = nthinobs+1

              lat_thin(nthinobs)=lat(i)
              lon_thin(nthinobs)=lon(i)
              sum_pblsld = sum_pblsld + pblsld(i) - sfc_elev(i)
              tobs_mvavg = tobs_mvavg + 1

              o1_lat = o2_lat
              o1_lon = o2_lon
           else
              cycle
           end if

        end if

        ! Consider near obs for moving average
        if (i <= no_obs_mvavg) then ! j=1,2,3,4,5
           begin=i*(-1)+1
           ending=no_obs_mvavg
        else if (i > nobs-no_obs_mvavg) then  !j=16,17,18,19,20
           begin=no_obs_mvavg*(-1)
           ending=nobs-i
        else
           begin=no_obs_mvavg*(-1)
           ending=no_obs_mvavg
        end if

        ! for +- obs near thinned obs
        do k = begin,ending
           if (k.ne.0) then
              o2_lat = lat(i+k)
              o2_lon = lon(i+k)
              dist = gc_dist(o1_lat,o1_lon,o2_lat,o2_lon)*0.001_r_kind
              if (dist < mvavg_rad) then ! select obs for moving average within radius 25 km
                 sum_pblsld = sum_pblsld + pblsld(i+k) - sfc_elev(i)
                 tobs_mvavg = tobs_mvavg + 1
              end if
           end if
        end do

        ! Save moving average pblh
        pblsld_mvavg_thin(nthinobs) = sum_pblsld / tobs_mvavg
        print*, 'YEG_CALIPSO_MVAVG: nthinobs=',nthinobs
        print*, 'YEG_CALIPSO_MVAVG: sum_pblsld,tobs_mvavg,pblh_avg=',sum_pblsld,tobs_mvavg,pblsld_mvavg_thin(nthinobs)

     end do

     end subroutine count_elev_mvavg_thin_obs

     subroutine count_qc_mvavg_thin_obs(nobs,lat,lon,qc,pblsld,mesh,mvavg_rad,n_qc_thin,lat_qc_thin,lon_qc_thin,pblsld_qc_mvavg_thin)
     use kinds, only: i_kind,r_kind
     use buddycheck_mod, only: gc_dist
     implicit none
     integer(i_kind), intent(in)  :: nobs ! the number of obs before thinning
     real(r_kind), dimension(nobs), intent(in)  :: lat, lon, qc, pblsld
     real(r_kind), dimension(:), allocatable :: lat_qc, lon_qc, pblsld_qc
     real(r_kind), dimension(:), allocatable, intent(out) :: lat_qc_thin, lon_qc_thin
     real(r_kind), dimension(:), allocatable, intent(out)  :: pblsld_qc_mvavg_thin
     real(r_kind), intent(in)  :: mesh ! thinning size (km)
     real(r_kind), intent(in)  :: mvavg_rad! moving average radius
     integer(i_kind), intent(out) :: n_qc_thin ! the number of obs after thinning

     integer(i_kind) :: i,j,k
     integer(i_kind) :: begin,ending
     integer(i_kind) :: n_qc, no_obs_mvavg
     real(r_kind) :: o1_lat, o1_lon, o2_lat, o2_lon, dist, dist_prev
     real(r_kind) :: tobs_mvavg, sum_pblsld
 
     ! FOR CATS 
     ! Count the number of obs after
     ! 1) Considering qc_flag>1 (good obs)
     ! 2) Thinning depending on distance (mesh)
     ! 3) calculating moving average within radius (mvavg_rad) for thinned obs

     n_qc=0

     ! 1) Save QC>1 (good obs)

     ! Count the number for qc>1 (n_qc), qc<=1 (n_bqc)
     do i = 1, nobs
        if (qc(i) > 1) then
           n_qc=n_qc+1
        end if 
     end do

     allocate(lat_qc(n_qc),lon_qc(n_qc),pblsld_qc(n_qc))

     n_qc=0
     ! Create lat, lon, qc for good qc
     do i = 1, nobs
        if (qc(i) > 1) then
           n_qc=n_qc+1
           lat_qc(n_qc) = lat(i)
           lon_qc(n_qc) = lon(i)
           pblsld_qc(n_qc) = pblsld(i)
        end if 
     end do

     print*, "YEG_CATS_QC_THIN, 1) total_obs,n_qc(qc>1)=",nobs,n_qc

     ! 2) Thinning for qc > 1 obs (good obs)

     ! Count the number of thinned obs for qc>1
     do i = 1, n_qc
        if (i == 1) then
           o1_lat = lat_qc(i)
           o1_lon = lon_qc(i)
           n_qc_thin=1
        else
           o2_lat = lat_qc(i)
           o2_lon = lon_qc(i)
           !function gc_dist(inlat1,inlon1,inlat2,inlon2) [meters]
           !lat1,lon1,lat2,lon2 in degrees
           ! *0.001_r_kind
           dist = gc_dist(o1_lat,o1_lon,o2_lat,o2_lon)*0.001_r_kind

           if (dist >= mesh) then
              n_qc_thin = n_qc_thin+1
              o1_lat = o2_lat
              o1_lon = o2_lon
           end if

        end if
     end do

     allocate(lat_qc_thin(n_qc_thin),lon_qc_thin(n_qc_thin),pblsld_qc_mvavg_thin(n_qc_thin))

     ! 3) calculating moving average within radius (mvavg_rad) for thinned obs
     ! Save thinned QCed obs (lat,lon,pblsld)
     no_obs_mvavg=10 ! 5 !10 !5 ! moving average --> +- 10 obs (25 km) along the track
     tobs_mvavg=0 ! total number of obs for moving avg
     sum_pblsld=0

     do i = 1, n_qc

        if (i == 1) then
           n_qc_thin=1
           j=i
           o1_lat = lat_qc(i)
           o1_lon = lon_qc(i)

           lat_qc_thin(n_qc_thin) = lat_qc(i)
           lon_qc_thin(n_qc_thin) = lon_qc(i)
           sum_pblsld = sum_pblsld + pblsld_qc(i)
           tobs_mvavg = tobs_mvavg + 1
        else
           !--------------------------
           ! Decide thinning or not
           o2_lat = lat_qc(i)
           o2_lon = lon_qc(i)
           !function gc_dist(inlat1,inlon1,inlat2,inlon2) [meters]
           !lat1,lon1,lat2,lon2 in degrees
           ! *0.001_r_kind
           dist = gc_dist(o1_lat,o1_lon,o2_lat,o2_lon)*0.001_r_kind

           if (dist >= mesh) then

              tobs_mvavg=0 ! total number of obs for moving avg
              sum_pblsld=0
              n_qc_thin = n_qc_thin+1
              j=i

              lat_qc_thin(n_qc_thin) = lat_qc(i)
              lon_qc_thin(n_qc_thin) = lon_qc(i)
              sum_pblsld = sum_pblsld + pblsld_qc(i)
              tobs_mvavg = tobs_mvavg + 1

              o1_lat = o2_lat
              o1_lon = o2_lon

           else
              cycle
           end if
        end if

        ! Consider near obs for moving average
        if (j <= no_obs_mvavg) then ! j=1,2,3,4,5
           begin=j*(-1)+1
           ending=no_obs_mvavg
        else if (j > n_qc-no_obs_mvavg) then  !j=16,17,18,19,20
           begin=no_obs_mvavg*(-1)
           ending=n_qc-j
        else
           begin=no_obs_mvavg*(-1)
           ending=no_obs_mvavg
        end if

        ! for +- obs near thinned obs
        do k = begin,ending
           if (k.ne.0) then
              o2_lat = lat_qc(j+k)
              o2_lon = lon_qc(j+k)
              dist = gc_dist(o1_lat,o1_lon,o2_lat,o2_lon)*0.001_r_kind
              if (dist < mvavg_rad) then ! select obs for moving average within radius 25 km
                 sum_pblsld = sum_pblsld + pblsld_qc(j+k)
                 tobs_mvavg = tobs_mvavg + 1
              end if
           end if
        end do

        ! Save moving average pblh
        pblsld_qc_mvavg_thin(n_qc_thin) = sum_pblsld / tobs_mvavg
        print*, 'YEG_CATS_MVAVG: n_qc_thin=',n_qc_thin
        print*, 'YEG_CATS_MVAVG: sum_pblsld,tobs_mvavg,pblh_avg=',sum_pblsld,tobs_mvavg,pblsld_qc_mvavg_thin(n_qc_thin)

     end do

     print*, "YEG_CATS_QC_THIN, 2) n_qc_thin=",n_qc_thin

     end subroutine count_qc_mvavg_thin_obs

  end module read_pblsld
