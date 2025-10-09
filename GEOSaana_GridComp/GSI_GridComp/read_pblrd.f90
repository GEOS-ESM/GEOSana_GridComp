!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.3, GMAO  !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: gsimod  ---

!
! !INTERFACE:
!
  module read_pblrd

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

      public read_pblrd_rwp
!---------------------------------------------------------------------------

   CONTAINS

     subroutine read_pblrd_rwp(nread,ndata,nodata,infile,obstype,lunout,twindin,&
         sis,nobs)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  read_pblrd_rwp                      read obs from msgs in NETCDF files
!
! program history log:
!   2022-10     Y. Zhu - adapted from read_pblh
!   2023-10     E.-G. Yang - added Radar Wind Profiler
!   2023-10     E.-G. Yang - changing the way of converting UTC to local time for obs
!                            Except for hourly obs, local time is obtained by adding float(lon/15) to UTC (YYYYMMDD HH:MI)
!                            For hourly obs (MPLNET, Wind Profiler), local time: int(lon/15) + UTC (YYYYMMDD HH:MI)
!   2023-10     E.-G. Yang - added qc for WP
!   2023-11     E.-G. Yang - added mixed sfc type qc for all obs 
!   2023-12     E.-G. Yang - added night-time qc for WP (6 PM - 10 AM)
!   2024-02     Y. Zhu - modelling observation error
!   2024-08     E.-G. Yang - added different obs error model for GRWP
!   2024-08     E.-G. Yang - changed local time QC for GRWP
!   2024-08     E.-G. Yang - added inflation of obs error for GRWP depending on averaged obs number (not varying one)
!   2024-08     E.-G. Yang - changed obs error 

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
      integer(i_kind) pblrdqm,maxobs,idomsfc
      integer(i_kind) oltime
      integer(i_kind) hr,ltime_mm,ltime_dd,ltime_hh,ltime_min
      real(r_kind) hr_min
      real(r_kind) usage
      real(r_kind) pblrdob,pblrdoe,pblrdelev,pblrdosfc
      real(r_kind) dlat,dlon,dlat_earth,dlon_earth
      real(r_kind) dlat_earth_deg,dlon_earth_deg
      real(r_kind) cdist,disterr,disterrmax,rlon00,rlat00
      real(r_kind) :: tsavg,ff10,sfcr,zz
      real(r_kind) dx,dy,dx1,dy1,w00,w10,w01,w11

      integer(i_kind) idate5(5),minobs,minan
      real(r_kind) time_correction,timeobs,time,toff,t4dv,zeps

!     real(r_single) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblrd_calipso
      real(r_kind) stationid,lat_deg,lon_deg,altitude,localtime,utctime,localday,utcday,pblrd_obs_kind

      integer(i_kind) :: ncid,ierr,dimid1,dimid2,dimid3,norbits,nheights
      integer(i_kind) :: noobs, strlen, site_maxstrlen, instrument_maxstrlen
      integer(i_kind) :: varid0,varid1,varid2,varid3,varid4,varid5,varid6
      integer(i_kind) :: varid7,varid8,varid9,varid10,varid11,varid12,varid13
      integer(i_kind) :: varid14, varid15
      integer(i_kind) :: ioyr, iomo, iody, iohr, iomi, iodate
      integer(i_kind), dimension(1) :: ana_time
      real(r_kind), allocatable, dimension(:) :: lat, lon, pblrd
      real(r_kind), allocatable, dimension(:) :: ryear, rmonth, rdate, rhour, rminute
      ! wind profiler

      real(r_double), allocatable, dimension(:) :: rsecond
      real(r_kind),allocatable,dimension(:) :: surface_alti, cloud_fraction, error, negative_snr_indicator
      character(5),allocatable,dimension(:) :: sid_wp
      character(5),allocatable,dimension(:,:) :: sid_wp_org
      character(5) :: c_station_id_wp
      real(r_double) :: r_station_id_wp
      ! Observation error model
      real(r_kind) :: cpblh1, cpblh2, cerr1a, cerr2a
      real(r_kind) :: cerr1b, cerr2b
      real(r_kind) :: cerr1c, cerr2c
      real(r_kind) :: cerr1d, cerr2d

      logical :: selected
!     equivalence to handle character names
      equivalence(r_station_id_wp,c_station_id_wp)

!     Initialize obs err parameters
      cpblh1 = 2000.0_r_kind
      cpblh2 = 4000.0_r_kind
      cerr1b = 250.0_r_kind ! grwp
      cerr2b = 700.0_r_kind

!     Check if pblrd file exists
      inquire(file=trim(infile),exist=lexist)
      if (.not.lexist) return

!     Read data
      ierr =  NF90_OPEN(trim(infile),0,ncid)
      if (ierr /= nf90_noerr) call handle_err(ierr,"open")

!     Read dimensions depending on observation types
      if (index(trim(infile),'wp') > 0) then ! Wind profiler

         kx = 893 ! wind profiler 
         nkx = 893 ! wind profiler

         ierr = NF90_INQ_DIMID(ncid,'nobs',dimid1)
         if (ierr /= nf90_noerr) call handle_err(ierr,"nobs")
         ierr = nf90_inquire_dimension(ncid, dimid1, len = noobs)

         ierr = NF90_INQ_DIMID(ncid,'Sid_maxstrlen',dimid2)
         if (ierr /= nf90_noerr) call handle_err(ierr,"Sid_maxstrlen")
         ierr = nf90_inquire_dimension(ncid, dimid2, len = strlen)

         print*, 'read_pblrd,',trim(infile),': noobs=', noobs

      end if

      select case(kx) ! Temporary setting
         case (893) ! wp
            maxobs=noobs
      end select

!     Allocate
      allocate(lat(maxobs), lon(maxobs), pblrd(maxobs))
      allocate(ryear(maxobs), rmonth(maxobs), rdate(maxobs), rhour(maxobs), rminute(maxobs))
      select case(kx)
         case (893) ! wind profiler
            allocate(rsecond(noobs), surface_alti(noobs), cloud_fraction(noobs), error(noobs), negative_snr_indicator(noobs))
            allocate(sid_wp_org(noobs,strlen), sid_wp(noobs))
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
         case (893) ! wind profiler
            ierr = NF90_INQ_VARID(ncid,'PBL_Height',varid8)
      end select

      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid8,pblrd)

!     Other variables
      select case(kx)

         case (893) ! wind profiler

            ierr = NF90_INQ_VARID(ncid,'Second',varid9)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid9,rsecond)
            ierr = NF90_INQ_VARID(ncid,'Altitude',varid10)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid10,surface_alti)
            ierr = NF90_INQ_VARID(ncid,'Cloud_Fraction',varid11)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid11,cloud_fraction)
            ierr = NF90_INQ_VARID(ncid,'Vert_Resolution_2x',varid12)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid12,error)
            ierr = NF90_INQ_VARID(ncid,'Negative_SNR_Indicator',varid13)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid13,negative_snr_indicator)

            ierr = NF90_INQ_VARID(ncid,'Station_ID',varid14)
            if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid14,sid_wp_org)

            do i = 1, noobs
               sid_wp(i)=sid_wp_org(i,1)
            end do

      end select

      ierr = NF90_CLOSE(ncid)
      if (ierr /= nf90_noerr) call handle_err(ierr,"close")

      nreal=17 ! 16 ! 15 ! Temporary setting

      allocate(cdata_all(nreal,maxobs))

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

         if (first) then
            !call time_4dvar(idate,toff) 
            call time_4dvar(ana_time(1),toff) 
            ! toff: time since start of 4D-Var window (hours)
            first=.false.
         end if

         nread=nread+1
         if(kx == 893) nkx= 893

!        Is pblrd in convinfo file
         ikx=0
         loop_convinfo_pblrd: do j=1,nconvtype
           if (trim(ioctype(j)) /= trim(obstype))cycle loop_convinfo_pblrd
           if (nkx == ictype(j)) then
              ikx=j
              exit
           end if
         end do loop_convinfo_pblrd

         print*, "YEG_read_pblrd L819:nkx=",nkx,", ikx=",ikx
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
         print*, 'read_pblrd: iodate, iadate, timeobs, toff, t4dv=', iodate, iadate, timeobs, toff, t4dv
         print*, 'YEG GRID: t4dv =', t4dv

         !------------------------------------------
         ! Local Time (only available for 2015/Aug)
         !------------------------------------------
         ! Converting UTC to Local Time (UTC + hr)
         ! longitude (degree East) -> East: lon>0, West: lon<0
         ! 1) hourly obs (MPLNET and Radar wind profiler)
         ! 2) non-hourly obs (CALIPSO, CATS, raob)

         select case(kx)

           case(893) ! RWP <- int(lon/15)
               ! Converting UTC to Local Time (UTC + hr)
               ! longitude (degree East) -> East: lon>0, West: lon<0
               if (lon(i)>180) then
                  hr = nint((lon(i)-360)/15.)
               else
                  hr = nint(lon(i)/15.)
               end if

               call convert_localtime_hr(ltime_mm, ltime_dd, ltime_hh, hr, iomo, iody, iohr)
               oltime = ioyr*1000000 + ltime_mm*10000 + ltime_dd*100 + ltime_hh

         end select

!        Get information from surface file necessary for conventional data here
         call deter_sfc2(dlat_earth,dlon_earth,t4dv,idomsfc,tsavg,ff10,sfcr,zz)
         !-----------------------
         ! Sfc Elevation from obs (CALIPSO, MPLNET, RAOB) or model at obs location (CATS)
         !-----------------------
         !pblrdelev=sfc_elev(i)
         select case(kx)
            case(893)
               pblrdelev=surface_alti(i)
         end select

         pblrd_obs_kind=pblrd(i)
         !------------------------------
         ! Obs sfc type
         ! Make the consistent land_water_mask for all obs types
         ! pblrdosfc: 1=Land, 2=water, 3=ice and snow, 4=mixed (coastlines)
         !------------------------------

         select case(kx)
            case(893) ! wp
               pblrdosfc=r_missing ! Check if this is okay !!!!!!!
         end select

         !----------------------
         ! QC
         !----------------------
         ! Initialize pblrdqm
         pblrdqm=0

         !------------
         ! RWP
         !------------
         ! Radar wind profiler: SNR qc flag = 0 good obs
         !if (kx == 893 .and. negative_snr_indicator(i).ne.0) then ! WP QC 0: good flag
         if (kx == 893) then
            if (negative_snr_indicator(i).ne.0) then
               pblrdqm=9
            ! QC for local time (no DA 5PM-9AM, time > 4 pm or < 10 am) 
            !else if (mod(oltime,100) < 10 .or. mod(oltime,100) > 16) then
            !   pblrdqm=10
            end if
         end if

         !----------
         ! ALL OBS
         !----------
         if (pblrd_obs_kind .le. 0.0 .or. pblrd_obs_kind .gt. 6000_r_kind) then
            pblrdqm=15
            print*, 'pblrdqm=15, pblh<0 or pblh>6000, pblrd_obs_kind=',pblrd_obs_kind
         else if (idomsfc>=3) then ! dominate model sfc type >=3 <-- mixed if any surrounding grid has different sfc type (water(0),land(1),ice(2))
            pblrdqm=20
            print*, 'pblrdqm=20, idomsfc>=3 (mixed sfc type), idomsfc=',idomsfc
         end if

         if (pblrd_obs_kind <= zero) then
            pblrdob=r_missing
         else
            pblrdob=pblrd_obs_kind
         end if

!        Set usage variable
         usage = 0.
         print*, 'YEG18: icuse(ikx)=',icuse(ikx), ',ikx=',ikx
         if(icuse(ikx) <= 0) usage=150.
         ! qm=15 for all obs (when pblh < 0 or pblh > 6000)
         ! qm=20 for all obs (when model sfc type = mixed)
         ! qm=9 (qc flag),10 (night-tiime) for wind profiler
         if(pblrdqm == 9 .or. pblrdqm == 10 .or. pblrdqm == 15 .or. pblrdqm == 20) usage=150.
         print*, 'YEG19: usage=',usage

!        Set inflate_error logical 
         inflate_error=.false.
         if (pblrdqm == 3) inflate_error=.true.

!--Outputs

         ndata=ndata+1
         nodata=nodata+1
         iout=ndata

         if(ndata > maxobs) then
            write(6,*)'READ_PBLRI:  ***WARNING*** ndata > maxobs for ',obstype
            ndata = maxobs
         end if

!        setup for sort of averaged obs 
         !---------------
         ! Obs Error
         !---------------
         !pblrdoe=0.1_r_kind  ! temporarily 100/1400
         !pblrdoe=1.31_r_kind  ! temporarily 100/1400
         !pblrdoe=1.5_r_kind ! 2.0_r_kind  ! Mar/10/2023
         !pblrdoe=0.075_r_kind ! 2.0_r_kind  ! Mar/10/2023

         !pblrdoe=0.1_r_kind*pblrdob ! 1.5_r_kind was used
         !pblrdoe=0.05_r_kind*pblrdob ! 1.5_r_kind was used
         if (pblrdob == r_missing) then
            pblrdoe = r_missing
         else
            if(pblrdob <= cpblh1) then
               pblrdoe = cerr1b
            else if(pblrdob > cpblh1 .and. pblrdob < cpblh2) then
               pblrdoe = cerr1b + (pblrdob-cpblh1)* &
                        (cerr2b-cerr1b)/(cpblh2-cpblh1)
            else
               pblrdoe = cerr2b
            endif
         endif         
    
         if (nkx==893) then ! MPLNET or GRWP

            !pblrdoe=1.0_r_kind ! smaller than BEC PBLRI (1.2)
!            pblrdoe=0.6_r_kind ! smaller than BEC PBLRI (1.2): latest
            ! inflate obs error of MPLNET according to no. of obs in DA window at each station
!            if (pblrd_no_hr_in_da_window(i) > 1.0) then
!               !sqrt 2 to 6 -> 1.4, 1.7, 2.0, 2.2, 2.4
!               pblrdoe=pblrdoe*sqrt(pblrd_no_hr_in_da_window(i))
!               print*, 'mplnet_oerror=',pblrdoe
!               print*, 'mplnet_no_obs=',pblrd_no_hr_in_da_window(i)
!            end if
         ! inflate obs error of MPLNET or GRWP according to averaged no. of obs in DA window
            pblrdoe=pblrdoe*sqrt(3.0_r_kind)
            if (nkx==893) print*, 'grwp_oerror=',pblrdoe
         end if

         if (inflate_error) pblrdoe=pblrdoe*1.05_r_kind

         ! Station id for wp
         if (nkx==893) then
            c_station_id_wp=sid_wp(i)
         end if

         cdata_all(1,iout)=pblrdoe                 ! pblrd error (cb)
         cdata_all(2,iout)=dlon                    ! grid relative longitude
         cdata_all(3,iout)=dlat                    ! grid relative latitude
         cdata_all(4,iout)=pblrdelev               ! pblrd obs elevation
         cdata_all(5,iout)=pblrdob                 ! pblrd obs
         cdata_all(6,iout)=r_station_id_wp         ! index of station id
         cdata_all(7,iout)=t4dv                    ! time
         cdata_all(8,iout)=ikx                     ! type
         cdata_all(9,iout)=oltime                  ! local time (YYYYMMDDHH) for all obs types!!!!!!!!!!!!!!!!
         cdata_all(10,iout)=pblrdqm                ! quality mark
         print*, 'YEG_read_pblrd L1090: pblrdqm=',pblrdqm,',usage=',usage
         cdata_all(11,iout)=usage                  ! usage parameter
         cdata_all(12,iout)=dlon_earth_deg         ! earth relative longitude (degrees)
         cdata_all(13,iout)=dlat_earth_deg         ! earth relative latitude (degrees)

         !sfc type
         cdata_all(14,iout)=pblrdosfc
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
 
      deallocate(lat, lon, pblrd)
      deallocate(ryear, rmonth, rdate, rhour, rminute)
      select case(kx)
         case(893) ! wp
            deallocate(surface_alti, rsecond)
            deallocate(cloud_fraction, error, negative_snr_indicator, sid_wp, sid_wp_org)
      end select

      end subroutine read_pblrd_rwp

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


     subroutine count_elev_mvavg_thin_obs(nobs,lat,lon,sfc_elev,pblrd,mesh,mvavg_rad,nthinobs,lat_thin,lon_thin,pblrd_mvavg_thin)
     use kinds, only: i_kind,r_kind
     use buddycheck_mod, only: gc_dist
     implicit none
     integer(i_kind), intent(in)  :: nobs ! the number of obs before thinning
     real(r_kind), dimension(nobs), intent(in)  :: lat, lon
     real(r_kind), dimension(nobs), intent(in)  :: sfc_elev,pblrd
     real(r_kind), intent(in)  :: mesh ! thinning size (km)
     real(r_kind), intent(in)  :: mvavg_rad! moving average radius
     integer(i_kind), intent(out) :: nthinobs ! the number of obs after thinning
     real(r_kind), dimension(:), allocatable, intent(out) :: lat_thin, lon_thin
     real(r_kind), dimension(:), allocatable, intent(out)  :: pblrd_mvavg_thin

     integer(i_kind) :: i,j,k
     integer(i_kind) :: begin,ending
     integer(i_kind) :: no_obs_mvavg
     real(r_kind) :: o1_lat, o1_lon, o2_lat, o2_lon, dist
     real(r_kind) :: tobs_mvavg, sum_pblrd

     ! FOR CALIPSO 
     ! Count the number of obs after
     ! 1) Thinning depending on distance (mesh)
     ! 2) Substracting sfc_elev from pblrd CALIPSO
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

     allocate(lat_thin(nthinobs),lon_thin(nthinobs),pblrd_mvavg_thin(nthinobs))

     ! 2) Save thinned obs and substract sfc_elev from pblrd
     ! 3) calculating moving average within radius (mvavg_rad) for thinned obs
     ! Save thinned QCed obs (lat,lon,pblrd)
     no_obs_mvavg=80 !40 !80 !40 ! moving average --> +- 80 obs along the track (25 km, 330m between obs)
     tobs_mvavg=0 ! total number of obs for moving avg
     sum_pblrd=0

     do i = 1, nobs
        if (i == 1) then
           nthinobs=1
           o1_lat = lat(i)
           o1_lon = lon(i)

           lat_thin(nthinobs)=lat(i)
           lon_thin(nthinobs)=lon(i)
           sum_pblrd = sum_pblrd + pblrd(i) - sfc_elev(i) ! subtracting sfc height
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
              sum_pblrd=0
              nthinobs = nthinobs+1

              lat_thin(nthinobs)=lat(i)
              lon_thin(nthinobs)=lon(i)
              sum_pblrd = sum_pblrd + pblrd(i) - sfc_elev(i)
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
                 sum_pblrd = sum_pblrd + pblrd(i+k) - sfc_elev(i)
                 tobs_mvavg = tobs_mvavg + 1
              end if
           end if
        end do

        ! Save moving average pblh
        pblrd_mvavg_thin(nthinobs) = sum_pblrd / tobs_mvavg
        print*, 'YEG_CALIPSO_MVAVG: nthinobs=',nthinobs
        print*, 'YEG_CALIPSO_MVAVG: sum_pblrd,tobs_mvavg,pblh_avg=',sum_pblrd,tobs_mvavg,pblrd_mvavg_thin(nthinobs)

     end do

     end subroutine count_elev_mvavg_thin_obs

     subroutine count_qc_mvavg_thin_obs(nobs,lat,lon,qc,pblrd,mesh,mvavg_rad,n_qc_thin,lat_qc_thin,lon_qc_thin,pblrd_qc_mvavg_thin)
     use kinds, only: i_kind,r_kind
     use buddycheck_mod, only: gc_dist
     implicit none
     integer(i_kind), intent(in)  :: nobs ! the number of obs before thinning
     real(r_kind), dimension(nobs), intent(in)  :: lat, lon, qc, pblrd
     real(r_kind), dimension(:), allocatable :: lat_qc, lon_qc, pblrd_qc
     real(r_kind), dimension(:), allocatable, intent(out) :: lat_qc_thin, lon_qc_thin
     real(r_kind), dimension(:), allocatable, intent(out)  :: pblrd_qc_mvavg_thin
     real(r_kind), intent(in)  :: mesh ! thinning size (km)
     real(r_kind), intent(in)  :: mvavg_rad! moving average radius
     integer(i_kind), intent(out) :: n_qc_thin ! the number of obs after thinning

     integer(i_kind) :: i,j,k
     integer(i_kind) :: begin,ending
     integer(i_kind) :: n_qc, no_obs_mvavg
     real(r_kind) :: o1_lat, o1_lon, o2_lat, o2_lon, dist, dist_prev
     real(r_kind) :: tobs_mvavg, sum_pblrd
 
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

     allocate(lat_qc(n_qc),lon_qc(n_qc),pblrd_qc(n_qc))

     n_qc=0
     ! Create lat, lon, qc for good qc
     do i = 1, nobs
        if (qc(i) > 1) then
           n_qc=n_qc+1
           lat_qc(n_qc) = lat(i)
           lon_qc(n_qc) = lon(i)
           pblrd_qc(n_qc) = pblrd(i)
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

     allocate(lat_qc_thin(n_qc_thin),lon_qc_thin(n_qc_thin),pblrd_qc_mvavg_thin(n_qc_thin))

     ! 3) calculating moving average within radius (mvavg_rad) for thinned obs
     ! Save thinned QCed obs (lat,lon,pblrd)
     no_obs_mvavg=10 ! 5 !10 !5 ! moving average --> +- 10 obs (25 km) along the track
     tobs_mvavg=0 ! total number of obs for moving avg
     sum_pblrd=0

     do i = 1, n_qc

        if (i == 1) then
           n_qc_thin=1
           j=i
           o1_lat = lat_qc(i)
           o1_lon = lon_qc(i)

           lat_qc_thin(n_qc_thin) = lat_qc(i)
           lon_qc_thin(n_qc_thin) = lon_qc(i)
           sum_pblrd = sum_pblrd + pblrd_qc(i)
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
              sum_pblrd=0
              n_qc_thin = n_qc_thin+1
              j=i

              lat_qc_thin(n_qc_thin) = lat_qc(i)
              lon_qc_thin(n_qc_thin) = lon_qc(i)
              sum_pblrd = sum_pblrd + pblrd_qc(i)
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
                 sum_pblrd = sum_pblrd + pblrd_qc(j+k)
                 tobs_mvavg = tobs_mvavg + 1
              end if
           end if
        end do

        ! Save moving average pblh
        pblrd_qc_mvavg_thin(n_qc_thin) = sum_pblrd / tobs_mvavg
        print*, 'YEG_CATS_MVAVG: n_qc_thin=',n_qc_thin
        print*, 'YEG_CATS_MVAVG: sum_pblrd,tobs_mvavg,pblh_avg=',sum_pblrd,tobs_mvavg,pblrd_qc_mvavg_thin(n_qc_thin)

     end do

     print*, "YEG_CATS_QC_THIN, 2) n_qc_thin=",n_qc_thin

     end subroutine count_qc_mvavg_thin_obs

     subroutine hourly_mean_obs(nobs,lat,lon,ryear,rmonth,rdate,rhour,rminute,pblrd,& 
                              surface_alti,flag_cloud_screen,site,instrument,use_cloud,&
                              nthinobs,lat_hr,lon_hr,ryear_hr,rmonth_hr,rdate_hr,rhour_hr,rminute_hr,pblrd_hr,&
                              surface_alti_hr,flag_use_cloud_hr,site_hr,instrument_hr,pblrd_no_hr_in_da_window)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:  hourly_mean_obs    calculate hourly mean obs for MPLNET PBLH
!
! program history log:
!   2023-07-04  eyang   -  
!
!   input argument list:
!     nobs          - the number of obs
!     lat
!     lon
!     ryear
!     rmonth
!     rdate
!     rhour
!     rminute
!     pblrd           - pbl height
!     surface_alti    - surface altitude [meters]
!     flag_cloud_screen - 1: cloud free 2: cloud_fraction>0% 4: cloud_detection_fail 8: no_could_product
!     site            - site information (e.g., GSFC)
!     instrument      - instrument information (e.g., MPL44104)
!     use_cloud       - use all data (including cloud-affected data) or use only cloud free
!                       0: consider only cloud free 1: consider all data (including cloud affected data)
!   output argument list:
!     nthinobs        - the number of thinned obs (temporally)
!     lat_hr
!     lon_hr
!     ryear_hr
!     rmonth_hr
!     rdate_hr
!     rhour_hr
!     rminute_hr
!     pblrd_hr        - pbl height
!     surface_alti_hr - surface altitude [meters]
!     flag_use_cloud_hr - 0: use only cloud free data 1: use all data available (including cloud affected data) 
!     site_hr         - site information (e.g., GSFC)
!     instrument_hr   - instrument information (e.g., MPL44104)
!     pblrd_no_hr_in_da_window   - no. of obs in DA window at each station (0 - 6)
!
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

     use kinds, only: i_kind,r_kind,r_double
     use, INTRINSIC :: IEEE_ARITHMETIC
     implicit none
     integer(i_kind), intent(in)  :: nobs ! the number of obs before thinning
     real(r_kind), dimension(nobs), intent(in)  :: lat, lon, ryear,rmonth,rdate,rhour,rminute, pblrd
     real(r_kind), dimension(nobs), intent(in)  :: surface_alti,flag_cloud_screen
     integer(i_kind), intent(in)  :: use_cloud
     character(19),dimension(nobs), intent(in)  :: site
     character(8), dimension(nobs), intent(in)  :: instrument
     real(r_kind), dimension(:), allocatable, intent(out)  :: lat_hr, lon_hr, ryear_hr,rmonth_hr,rdate_hr,rhour_hr
     real(r_kind), dimension(:), allocatable, intent(out)  :: rminute_hr, pblrd_hr
     real(r_kind), dimension(:), allocatable, intent(out)  :: surface_alti_hr,flag_use_cloud_hr
     character(19),dimension(:), allocatable, intent(out)  :: site_hr
     character(8), dimension(:), allocatable, intent(out)  :: instrument_hr
     real(r_kind), dimension(:), allocatable, intent(out)  :: pblrd_no_hr_in_da_window !no. of obs in DA window at each station (0 - 6)
     real(r_kind), dimension(:), allocatable :: lat_qc, lon_qc
     integer(i_kind), intent(out) :: nthinobs ! the number of obs after thinning

     character(19),dimension(nobs) :: site_unique_temp
     character(19),dimension(:), allocatable :: site_unique
     character(8),dimension(:), allocatable :: instrument_unique
     real(r_kind), dimension(6) :: hour_unique
     ! temp array
     real(r_kind), dimension(:,:), allocatable :: pblrd_temp, pblrd_no
     real(r_kind), dimension(:), allocatable :: pblrd_no_hr_in_da_window_temp
     real(r_kind), dimension(:,:), allocatable :: lat_temp, lon_temp, ryear_temp, rmonth_temp, rdate_temp
     real(r_kind), dimension(:,:), allocatable :: surface_alti_temp
     integer(i_kind) :: i,j,k,m
     integer(i_kind) :: site_count, hour_count
     real(r_kind) :: o1_lat, o1_lon

     ! Calculate hourly mean pblh for each station
     ! 1) Find unique sites and hours
     ! 2) collect data for hourly mean depending on sites
     ! 3) Calculate average
     ! 4) Assign hourly mean value to 1-d array

     ! 1) Find unique sites and hours

     site_unique_temp(1)=site(1)
     site_count=1

     do i = 2, nobs

        do j = 1, site_count

           if (site(i)==site_unique_temp(j)) then
              exit
           end if

        end do

        if (j == site_count+1) then
           site_count = site_count + 1
           site_unique_temp(site_count) = site(i)
        end if

     end do

     ! Store unique sites and hours
     ! a) sites and instruments
     allocate(site_unique(site_count))
     allocate(instrument_unique(site_count))
     site_unique(1)=site(1)
     instrument_unique(1)=instrument(1)
     site_count=1

     do i = 2, nobs

        do j = 1, site_count

           if (site(i)==site_unique(j)) then
              exit
           end if

        end do

        if (j == site_count+1) then
           site_count = site_count + 1
           site_unique(site_count) = site(i)
           instrument_unique(site_count) = instrument(i)
        end if

     end do

     ! b) hours
     hour_unique(1)=rhour(1)
     hour_count=1

     do i = 2, nobs

        do j = 1, hour_count

           if (rhour(i)==hour_unique(j)) then
              exit
           end if

        end do

        if (j == hour_count+1) then
           hour_count = hour_count + 1
           hour_unique(hour_count) = rhour(i)
        end if

     end do


     ! 2) collect data for hourly mean depending on sites
     allocate(pblrd_temp(site_count, 6))
     allocate(pblrd_no(site_count, 6))
     allocate(pblrd_no_hr_in_da_window_temp(site_count))
     allocate(lat_temp(site_count,6), lon_temp(site_count,6), ryear_temp(site_count,6),rmonth_temp(site_count,6))
     allocate(rdate_temp(site_count,6),surface_alti_temp(site_count,6))
     pblrd_temp = 0.
     pblrd_no = 0.

     ! read data 
     do i = 1, nobs ! all obs
        do j = 1, site_count ! sites
           do k = 1, 6 ! hour bin

              if (site(i) == site_unique(j) .and. rhour(i) == hour_unique(k)) then

                 lat_temp(j,k)=lat(i)
                 lon_temp(j,k)=lon(i)
                 ryear_temp(j,k)=ryear(i)
                 rmonth_temp(j,k)=rmonth(i)
                 rdate_temp(j,k)=rdate(i)
                 surface_alti_temp(j,k)=surface_alti(i)

                 if (IEEE_IS_NAN(pblrd(i))) then
                    cycle
                 else
                    if (use_cloud==1) then ! use all data (including cloud-affected data)
                       pblrd_temp(j, k) = pblrd_temp(j,k) + pblrd(i)
                       pblrd_no(j, k) = pblrd_no(j,k) + 1
                    else ! use_cloud==0 (use only cloud-free data)
                       if (flag_cloud_screen(i) == 1) then ! cloud-free
                          pblrd_temp(j, k) = pblrd_temp(j,k) + pblrd(i)
                          pblrd_no(j, k) = pblrd_no(j,k) + 1
                       end if
                    end if
                 end if

              end if

           end do
        end do
     end do

     ! 3) Calculate average and 
     !    the number of obs in DA window at each station
     !    pblrd_no_hr_in_da_window   - no. of obs in DA window at each station (0 - 6)
     pblrd_no_hr_in_da_window_temp=0.
     do j = 1, site_count
        do k = 1, 6 ! hour bin

           if (pblrd_no(j,k)==0) then
              pblrd_temp(j,k) = -9999
           else
              pblrd_temp(j,k) = pblrd_temp(j,k)/pblrd_no(j,k)
              pblrd_no_hr_in_da_window_temp(j) = pblrd_no_hr_in_da_window_temp(j) + 1
              print*, 'j=',j,',k=',k,',pblrdnoinda(j)=',pblrd_no_hr_in_da_window_temp(j)
           end if
           print*, 'j=',j,'k=',k,',pblrd_no=',pblrd_no(j,k)
           print*, 'j=',j,'k=',k,',pblrd_temp=',pblrd_temp(j,k)

        end do
        print*, 'j=',j,',pblrdnoinda(j)=',pblrd_no_hr_in_da_window_temp(j)
     end do
              
     ! 4) Assign hourly mean value to 1-d array 

     nthinobs=site_count*6
     allocate(lat_hr(nthinobs),lon_hr(nthinobs),ryear_hr(nthinobs),rmonth_hr(nthinobs),rdate_hr(nthinobs))
     allocate(rhour_hr(nthinobs),rminute_hr(nthinobs),pblrd_hr(nthinobs),surface_alti_hr(nthinobs),flag_use_cloud_hr(nthinobs))
     allocate(site_hr(nthinobs),instrument_hr(nthinobs), pblrd_no_hr_in_da_window(nthinobs))
     m=1
     do j = 1, site_count
        do k = 1, 6 ! hour bin

           lat_hr(m)=lat_temp(j,k)
           lon_hr(m)=lon_temp(j,k)
           ryear_hr(m)=ryear_temp(j,k)
           rmonth_hr(m)=rmonth_temp(j,k)
           rdate_hr(m)=rdate_temp(j,k)
           rhour_hr(m)=hour_unique(k)
           rminute_hr(m)=30.0
           pblrd_hr(m)=pblrd_temp(j,k)
           surface_alti_hr(m)=surface_alti_temp(j,k)
           flag_use_cloud_hr(m)=use_cloud ! 1: use all data 0: use only cloud-free data
           site_hr(m)=site_unique(j)
           instrument_hr(m)=instrument_unique(j)
           pblrd_no_hr_in_da_window(m)=pblrd_no_hr_in_da_window_temp(j)
           m=m+1

        end do
     end do

     end subroutine hourly_mean_obs

  end module read_pblrd
