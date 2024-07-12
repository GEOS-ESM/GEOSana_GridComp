!-------------------------------------------------------------------------
!  NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GMAO  !
!-------------------------------------------------------------------------
!BOP

! ! IPROGRAM: calculate pblh from radiosonde using diagnostic files

! ! INTERFACE:

   program calculate_pblh_from_raob_diag

   use constants, only: zero,one_tenth,one,deg2rad,rad2deg,three, half, r_missing
   use kinds, only: r_double,i_kind,r_kind
   use netcdf
   implicit none

!$$$  main program documentation block
!                .      .    .                                       .
! main program: GSI_ANL
!   PRGMMR: DERBER           ORG: NP23                DATE: 1999-08-20
!
! abstract: This program calculates the planetary boundary layer (PBL) height 
!   from radiosonde observations using diagnostic files. In this code, PBLH is
!   defined based on Bulk Richardson number (Ri). In addition, two additional
!   PBLH definitions are calculated: thetav-gradient and q-gradient PBLH
!   0. Read analysis time (as argument)
!   1. Read obs diagnostic files
!   2. Extract radiosonde observations
!      a) Extract radiosonde observations and convert obstime [0,6] to [-3,3] 
!         (merra-2: 201508-09 [0,6]->[-3,3] vs r21c: 202112 -> already [-3,3])
!      b) Extract unique raob station id (remove dups)
!         & sort obs in a descending order according to press
!      c) convert utc to local time using longitude
!         local time = nint(lon/15) + UTC
!   3. Extract and assign raob_q (specific humidity) and
!      raob_w (mixing ratio) for corresponding the same P-levels
!      of raob temperature (sid)
!      (vertically linear interpolation - log(pressure) )
!   4. Convert specific humidity to mixing ratio
!   5. Convert temperature to virtual temperature
!   6. Convert virtual temperature to virtual potential temperature
!   7. Convert temperature to potential temperature
!   8. Convert Pressure -> geopotential height
!   9. Vertically linearly interpolate observed P levels of wind
!      to observed P levels of T
!  10. Assign sfc level values to height, u, v, theta_v
!  11. Calculate Bulk Richardson number
!                   g * (theta_v_z - theta_v_sfc) * (z - z_sfc)
!      RiB = -------------------------------------------------------
!             (theta_v_sfc) * ( (u_z - u_sfc)**2 + (v_z - v_sfc)**2 )
!      here, surface wind is considered as zero likewise ‘pblri’ model variable.
!  12. Calculate thetav-gradient and q-gradient PBLH
!      a) Consider the levels where p>=500 hPa
!      b) averaged value between two p levels:
!         pblh_thetav, lat, lon, time
!      c) value at the lower level among two p levels:
!         ltime_pblh_thetav
!  13. Calculate RiB-based PBLH
!      a) QC flag: if less than 10 p-levels below 500 hPa
!                  then qc = 5, otherwise qc = 0
!      b) if ri > 0.25 for 2nd obs level,
!         pblh_ri should be assigned as height at 2nd level (AGL) and QC_Flag=50
!      c) vertically linearly interpolated value between two p levels:
!         pblh_ri, lat, lon, time
!      d) value at the lower level among two p levels:
!         ltime_pblh_ri
!  14. Save obstime for YYYY, MM, DD, HH, and Min format
!  15. Write nc file
!
! program history log:
!   2023-05-17 eyang    - change geometric altitude to geopotential height
!   2023-08-29 eyang    - read the part of diag nc file name as argument
!   2023-10-31 eyang    - remove pblh > 6000 m and assign it as qc=99
!   2023-12-05 eyang    - add the argument for conversion time (merra2: -3, r21c: +0) 
!   2023-12-20 eyang    - add one more condition to obtain q gradient-PBLH above 25 m
!   2024-01-07 eyang    - change minimum q-grad pblh from 25 m to 50 m
!   2024-01-09 eyang    - remove minimum q grad pblh
!
! usage:
!
!   input files or argument:
!
!     infile             - diagnostic file (diag_conv_VAR_ges.nc4, VAR: t, q, uv)
!     timechar_argument  - analysis time (this argument is required for the executable
!                                         e.g., raob_pblh.x 2015082300 r21c_pbl 3)
!     fnamechar_argument - the part of nc diagnostic file name 
!                                        (2nd argument is required for the executable
!                                         e.g., raob_pblh.x 2015082300 r21c_pbl 3)
!
!     offsetchar_argument - the part of time conversion (offset) (merra2: [0,6]->[-3,3], r21c: [-3,3]->[-3,3])
!                                        (3rd argument is required for the executable
!                                         e.g., raob_pblh.x 2015082300 r21c_pbl 3)
!
!   output files:  (including scratch files)
!
!     raob_pblh.TIME.nc4 - output radiosonde pblh file 
!
!   subprograms called:
!
!     modules:
!         From GSI:
!           constants, kinds
!
!     libraries:
!         NETCDF    - the netcdf library
!
! ! attributes:
!   language: f90
!
!$$$

!    NOTE: 
!
!
!==================================================================================================


!  Declare passed variables
   character(len=225):: infile,obstype
   character(len=2):: obsvar

!  Declare local parameters
   integer(i_kind)  :: raob_type_tq=120
   integer(i_kind)  :: raob_type_uv=220
   real(r_double)   :: tcri_crit=0.25
   character(len=225):: filename

!  Declare local variables
   integer(i_kind),allocatable,dimension(:) :: obs_type_t
   integer(i_kind),allocatable,dimension(:) :: obs_type_q
   integer(i_kind),allocatable,dimension(:) :: obs_type_uv
   integer(i_kind),allocatable,dimension(:) :: sid_plevels_t
   integer(i_kind),allocatable,dimension(:) :: sid_plevels_t_500
   integer(i_kind),allocatable,dimension(:) :: sid_plevels_q
   integer(i_kind),allocatable,dimension(:) :: sid_plevels_q_500
   integer(i_kind),allocatable,dimension(:) :: sid_plevels_uv
   integer(i_kind),allocatable,dimension(:) :: sid_plevels_uv_500
   real(r_double),allocatable,dimension(:)  :: lat_t, lon_t
   real(r_double),allocatable,dimension(:)  :: lat_q, lon_q
   real(r_double),allocatable,dimension(:)  :: lat_uv,lon_uv
   real(r_double),allocatable,dimension(:)  :: lat,lon ! raob
   real(r_double),allocatable,dimension(:)  :: stn_elev_t, pressure_t
   real(r_double),allocatable,dimension(:)  :: stn_elev_q, pressure_q
   real(r_double),allocatable,dimension(:)  :: stn_elev_uv,pressure_uv
   real(r_double),allocatable,dimension(:)  :: stn_elev,pressure ! raob
   real(r_double),allocatable,dimension(:)  :: height_t, prep_qc_mark_t
   real(r_double),allocatable,dimension(:)  :: height_q, prep_qc_mark_q
   real(r_double),allocatable,dimension(:)  :: height_uv,prep_qc_mark_uv
   real(r_double),allocatable,dimension(:)  :: height,prep_qc_mark! raob
   real(r_double),allocatable,dimension(:)  :: height_converted ! converted from pressure
   real(r_double),allocatable,dimension(:)  :: height_converted_sfc ! converted from pressure
   !real(r_double),allocatable,dimension(:)  :: geo_altitude_converted ! geometric altitude: converted from geo. height
   !real(r_double),allocatable,dimension(:)  :: geo_altitude_converted_sfc ! geometric altitude: converted from geo. height
   real(r_double),allocatable,dimension(:)  :: obs_t, obs_q
   real(r_double),allocatable,dimension(:)  :: obs_u, obs_v
   real(r_double),allocatable,dimension(:)  :: raob_t, raob_q ! raob
   real(r_double),allocatable,dimension(:)  :: raob_tv, raob_w ! raob
   real(r_double),allocatable,dimension(:)  :: raob_theta, raob_thetav ! potential temperature
   real(r_double),allocatable,dimension(:)  :: raob_thetav_sfc
   real(r_double),allocatable,dimension(:)  :: raob_u, raob_v ! raob
   real(r_double),allocatable,dimension(:)  :: raob_u_sfc, raob_v_sfc
   real(r_double),allocatable,dimension(:)  :: time_t, time_q, time_uv
   real(r_double),allocatable,dimension(:)  :: errinv_input_t
   real(r_double),allocatable,dimension(:)  :: errinv_input_q
   real(r_double),allocatable,dimension(:)  :: errinv_input_uv
   real(r_double),allocatable,dimension(:)  :: time ! raob
   integer(i_kind),allocatable,dimension(:) :: ltime ! local time = nint(lon/15)
   real(r_double),allocatable,dimension(:)  :: raob_ri
   real(r_double),allocatable,dimension(:)  :: raob_pblh_ri
   real(r_double),allocatable,dimension(:)  :: raob_pblh_q
   real(r_double),allocatable,dimension(:)  :: raob_pblh_thetav
   real(r_double),allocatable,dimension(:)  :: time_raob_pblh_ri
   real(r_double),allocatable,dimension(:)  :: time_raob_pblh_q
   real(r_double),allocatable,dimension(:)  :: time_raob_pblh_thetav
   integer(i_kind)                          :: ltime_mm,ltime_dd,ltime_hh,ltime_yy,ltime_min
   real(r_kind)                             :: hr_min
   integer(i_kind)                          :: ana_mm0,ana_dd0,ana_hh0,ana_yy0,ana_min0
   integer(i_kind)                          :: ana_mm1,ana_dd1,ana_hh1,ana_yy1,ana_min1
   integer(i_kind)                          :: time_argument
   integer(i_kind)                          :: offset_argument
   integer(i_kind)                          :: ana_time
   ! For each sid
!   integer(i_kind),allocatable,dimension(:) :: ana_ltime ! local time = nint(lon/15) + ana_time (UTC)
   real(r_double),allocatable,dimension(:)  :: pblh_ri
   real(r_double),allocatable,dimension(:)  :: pblh_q
   real(r_double),allocatable,dimension(:)  :: pblh_thetav
   real(r_double),allocatable,dimension(:)  :: pblh_q_lowlevel
   real(r_double),allocatable,dimension(:)  :: pblh_thetav_lowlevel
   real(r_double),allocatable,dimension(:)  :: time_pblh_ri
   real(r_double),allocatable,dimension(:)  :: time_pblh_q
   real(r_double),allocatable,dimension(:)  :: time_pblh_thetav
   real(r_double),allocatable,dimension(:)  :: yyyy_pblh_ri
   real(r_double),allocatable,dimension(:)  :: mm_pblh_ri
   real(r_double),allocatable,dimension(:)  :: dd_pblh_ri
   real(r_double),allocatable,dimension(:)  :: hh_pblh_ri
   real(r_double),allocatable,dimension(:)  :: min_pblh_ri
   integer(i_kind),allocatable,dimension(:) :: ltime_pblh_ri
   integer(i_kind),allocatable,dimension(:) :: ltime_pblh_q
   integer(i_kind),allocatable,dimension(:) :: ltime_pblh_thetav
   real(r_double),allocatable,dimension(:)  :: lat_pblh_ri
   real(r_double),allocatable,dimension(:)  :: lat_pblh_q
   real(r_double),allocatable,dimension(:)  :: lat_pblh_thetav
   real(r_double),allocatable,dimension(:)  :: lon_pblh_ri
   real(r_double),allocatable,dimension(:)  :: lon_pblh_q
   real(r_double),allocatable,dimension(:)  :: lon_pblh_thetav
   integer(i_kind),allocatable,dimension(:)  :: qc_mark_pblh
   real(r_double),allocatable,dimension(:)  :: stn_elev_pblh

   real(r_double),allocatable,dimension(:)  :: lat_so, lon_so, stn_elev_so
   real(r_double),allocatable,dimension(:)  :: pressure_so, height_so, prep_qc_mark_so
   real(r_double),allocatable,dimension(:)  :: raob_t_so, time_so 

   real(r_double),allocatable,dimension(:)  :: latq, lonq, stn_elevq, pressureq, heightq
   real(r_double),allocatable,dimension(:)  :: prep_qc_markq,raob_q_org,timeq
   real(r_double),allocatable,dimension(:)  :: latuv, lonuv, stn_elevuv, pressureuv, heightuv
   real(r_double),allocatable,dimension(:)  :: prep_qc_markuv,raob_u_org,raob_v_org,timeuv

   character(8),allocatable,dimension(:)  :: sid_t
   character(8),allocatable,dimension(:)  :: sid_q
   character(8),allocatable,dimension(:)  :: sid_uv
   character(8),allocatable,dimension(:)  :: sid ! raob
   character(8),allocatable,dimension(:)  :: sidq 
   character(8),allocatable,dimension(:)  :: siduv 
   character(8),allocatable,dimension(:)  :: sid_so ! sorted 
   character(8),allocatable,dimension(:)  :: raob_sid_list ! unique raob sid list
   character(8),allocatable,dimension(:)  :: raob_sid_list_q ! unique raob sid list
   character(8),allocatable,dimension(:)  :: raob_sid_list_uv ! unique raob sid list

   character(10)   :: timechar_argument
   character(20)   :: fnamechar_argument
   character(1)    :: offsetchar_argument
   character(12)   :: fn_yyyymmdd_hhz
   integer(i_kind) :: yyyymmdd
   integer(i_kind) :: date_time_t
   integer(i_kind) :: date_time_q
   integer(i_kind) :: date_time_uv
   integer(i_kind) :: n_raob_t ! the number of raob T
   integer(i_kind) :: n_raob_q ! the number of raob q
   integer(i_kind) :: n_raob_uv ! the number of raob UV
   integer(i_kind) :: i,j,k,m,q
   integer(i_kind) :: np ! the number of pressure level depending on sid
   integer(i_kind) :: location, indexx
   real(r_kind)    :: h, dz
!   real(r_kind)   :: geom_alti
   real(r_kind)    :: slope_q, slope_u, slope_v ! linear vertical interpolation (q and wind)
   real(r_kind)    :: height_sfc, u_sfc, v_sfc, thetav_sfc
   real(r_kind)    :: denom, numer
   real(r_kind)    :: pblh_ri_temporary
   real(r_kind)    :: maxdthvdz, dthvdz
   real(r_kind)    :: mindqdz, dqdz

   integer(i_kind) :: ncid,ierr,dimid1,dimid2,dimid3,dimid4
   integer(i_kind) :: varid1,varid2,varid3,varid4,varid5,varid6                         
   integer(i_kind) :: varid7,varid8,varid9
   integer(i_kind) :: varid10,varid11,varid12,varid13,varid14,varid15,varid16
   integer(i_kind) :: varid17,varid18,varid19
   integer(i_kind) :: varid20,varid21,varid22,varid23,varid24,varid25,varid26
   integer(i_kind) :: varid27,varid28,varid29
   integer(i_kind) :: varid30,varid31,varid32,varid33,varid34,varid35,varid36
   integer(i_kind) :: varid37,varid38,varid39
   integer(i_kind) :: varid40,varid41,varid42,varid43,varid44,varid45,varid46
   integer(i_kind) :: varid47,varid48,varid49
   integer(i_kind) :: varid50,varid51,varid52,varid53,varid54,varid55,varid56
   integer(i_kind) :: varid57,varid58,varid59

!  Declare derived constants
   real(r_kind):: rv, rd, cp, rd_over_cp, fv, grav,rdog

!  Numeric Constants
   real(r_kind),parameter:: r1000 = 1000.0_r_kind

!  Derived atmospheric constants
   rd     = 2.8705e+2_r_kind
   rv     = 4.6150e+2_r_kind
   cp     = 1.0046e+3_r_kind            !  specific heat of air @pressure  (J/kg/K) 
   grav   = 9.80665e+0_r_kind
   fv  = rv/rd-one
   rd_over_cp = rd/cp
   rdog = rd/grav
!  https://github.com/NCAR/rda-prepbufr-decode/blob/main/README.md      

!----------------------------------------------------------
!  0. Read analysis time and nc diag file name and conversion time (three arguments)
!----------------------------------------------------------

   ! Get argument for analysis time.
   !if (command_argument_count().ne.1) then
   if (command_argument_count().ne.3) then
      print*, "Error: Analysis time and nc diag file name are required as three arguments."
   end if

   call get_command_argument(1,timechar_argument)
   read(timechar_argument,*) time_argument
   call get_command_argument(2,fnamechar_argument)
   call get_command_argument(3,offsetchar_argument)
   read(offsetchar_argument,*) offset_argument

   print*, "==================================================="
   print*, "The argument for analysis time is", time_argument
   print*, "==================================================="

   ! Define
   ana_time = time_argument
   print*, "ana_time=", ana_time
!  set analysis time
   ana_hh0 = mod(ana_time, 100)
   ana_dd0 = ( mod(ana_time, 10000) - ana_hh0 ) / 100
   ana_mm0 = ( mod(ana_time, 1000000) - ana_dd0*100 - ana_hh0 ) / 10000
   ana_yy0 = ( ana_time - ana_mm0*10000 - ana_dd0*100 - ana_hh0 ) / 1000000
   ana_yy1 = ana_yy0

   yyyymmdd = ana_yy0*10000+ana_mm0*100+ana_dd0 

   if (ana_hh0<10) then
      write(fn_yyyymmdd_hhz, '(I8,A2,I1,A1)') yyyymmdd,'_0',ana_hh0,'z'
   else
      write(fn_yyyymmdd_hhz, '(I8,A1,I2,A1)') yyyymmdd,'_',ana_hh0,'z'
   end if

!--------------------------------
!  1. Read obs diagnostic files
!--------------------------------
!  Temperature
   obsvar='t' ! 'uv' !'t' ! 'q', 'uv'
   infile = trim(fnamechar_argument)//'.diag_conv_'//trim(obsvar)//'_ges.'//trim(fn_yyyymmdd_hhz)//'.nc4'
   !infile = 'd5124_m2_jan10.diag_conv_'//trim(obsvar)//'_ges.'//trim(fn_yyyymmdd_hhz)//'.nc4'
   !infile = 'd5124_m2_jan10.diag_conv_'//trim(obsvar)//'_ges.20150822_00z.nc4'
   print*, "Read file:", infile
   call read_diag(infile,obsvar,obs_type_t,lat_t,lon_t,stn_elev_t,pressure_t,&
        height_t,prep_qc_mark_t,obs_t,time_t,errinv_input_t,sid_t,date_time_t)

!  Specific humidity
   obsvar='q' ! 'uv' !'t' ! 'q', 'uv'
   infile = trim(fnamechar_argument)//'.diag_conv_'//trim(obsvar)//'_ges.'//trim(fn_yyyymmdd_hhz)//'.nc4'
   !infile = 'd5124_m2_jan10.diag_conv_'//trim(obsvar)//'_ges.'//trim(fn_yyyymmdd_hhz)//'.nc4'
   print*, "Read file:", infile
   call read_diag(infile,obsvar,obs_type_q,lat_q,lon_q,stn_elev_q,pressure_q,&
        height_q,prep_qc_mark_q,obs_q,time_q,errinv_input_q,sid_q,date_time_q)

!  Wind
   obsvar='uv' ! 'uv' !'t' ! 'q', 'uv'
   infile = trim(fnamechar_argument)//'.diag_conv_'//trim(obsvar)//'_ges.'//trim(fn_yyyymmdd_hhz)//'.nc4'
   !infile = 'd5124_m2_jan10.diag_conv_'//trim(obsvar)//'_ges.'//trim(fn_yyyymmdd_hhz)//'.nc4'
   print*, "Read file:", infile
   call read_diag(infile,obsvar,obs_type_uv,lat_uv,lon_uv,stn_elev_uv,pressure_uv,&
        height_uv,prep_qc_mark_uv,obs_u,time_uv,errinv_input_uv,sid_uv,date_time_uv,obs_v)

   print*, "---Complete Reading Diagnostic Files"

!-----------------------------
!  2_1. Extract radiosonde observations
!       & convert obstime [0,6] to [-3,3] (for merra2: 201508-09)
!        (c.f., for r21c (202112) no need to convert obstime, it is already [-3,3]
!-----------------------------

   ! Temperature
   n_raob_t=0
   do i = 1, size(obs_t)
      if (obs_type_t(i) == raob_type_tq) then ! if obs is raob
         n_raob_t=n_raob_t+1
      end if
   end do
   print*, "n_raob_t=", n_raob_t
 
   ! Specific humidity
   n_raob_q=0
   do i = 1, size(obs_q)
      if (obs_type_q(i) == raob_type_tq) then ! if obs is raob
         n_raob_q=n_raob_q+1
      end if
   end do
   print*, "n_raob_q=", n_raob_q
 
   ! Wind
   n_raob_uv=0
   do i = 1, size(obs_u)
      if (obs_type_uv(i) == raob_type_uv) then ! if obs is raob
         n_raob_uv=n_raob_uv+1
      end if
   end do
   print*, "n_raob_uv=", n_raob_uv
 
   ! allocate for T
   allocate(lat(n_raob_t), lon(n_raob_t), stn_elev(n_raob_t), pressure(n_raob_t), &
      height(n_raob_t), prep_qc_mark(n_raob_t), &
      raob_t(n_raob_t), raob_q(n_raob_t), & 
      time(n_raob_t), sid(n_raob_t))
   lat=r_missing
   lon=r_missing
   stn_elev=r_missing
   pressure=r_missing
   height=r_missing
   prep_qc_mark=r_missing
   raob_t=r_missing
   raob_q=r_missing
   time=r_missing
   sid=''
   ! allocate for q
   allocate(latq(n_raob_q), lonq(n_raob_q), stn_elevq(n_raob_q), pressureq(n_raob_q), &
      prep_qc_markq(n_raob_q), raob_q_org(n_raob_q), heightq(n_raob_q), &
      timeq(n_raob_q), sidq(n_raob_q))
   latq=r_missing
   lonq=r_missing
   stn_elevq=r_missing
   pressureq=r_missing
   heightq=r_missing
   prep_qc_markq=r_missing
   raob_q_org=r_missing
   timeq=r_missing
   sidq=''

   ! allocate for uv
   allocate(latuv(n_raob_uv), lonuv(n_raob_uv), stn_elevuv(n_raob_uv), pressureuv(n_raob_uv), &
      prep_qc_markuv(n_raob_uv), raob_u_org(n_raob_uv), raob_v_org(n_raob_uv),heightuv(n_raob_uv), &
      timeuv(n_raob_uv), siduv(n_raob_uv))
   latuv=r_missing
   lonuv=r_missing
   stn_elevuv=r_missing
   pressureuv=r_missing
   heightuv=r_missing
   prep_qc_markuv=r_missing
   raob_u_org=r_missing
   raob_v_org=r_missing
   raob_q=r_missing
   timeuv=r_missing
   siduv=''

   !allocate(raob_sid_list(n_raob_t),raob_sid_list_q(n_raob_q), raob_sid_list_uv(n_raob_uv))
 
   !---------------------------------------
   ! Extract only radiosonde observations
   ! & convert obstime [0,6] to [-3,3]
   ! for merra-2 (201508-09), '-3', but for r21c (202112), '0'
   !---------------------------------------

   ! For Temepeature
   j = 1
   do i = 1, size(obs_t)

      if (obs_type_t(i) == raob_type_tq) then ! if obs is raob

         lat(j)=lat_t(i)
         lon(j)=lon_t(i)
         stn_elev(j)=stn_elev_t(i)
         pressure(j)=pressure_t(i)

         if (height_t(i) > 1000000) then ! missing value
            height(j)=r_missing
         else
            height(j)=height_t(i)
         end if

         prep_qc_mark(j)=prep_qc_mark_t(i)
         !setup_qc_mark(j)=setup_qc_mark_t(i)
         raob_t(j)=obs_t(i)
         ! Convert obstime (hr) relative to the beginning of DA window
         ! to obstime relative to analysis time
         ! for merra-2 (201508-09), '-3', but for r21c (202112), '0'
         !time(j)=time_t(i)-3
         time(j)=time_t(i)-offset_argument
         sid(j)=sid_t(i)
         print*, "j=",j,",lat(j)=",lat(j),",lon(j)=",lon(j)
         print*, "        sid(j)=",sid(j),",pressure(j)=",pressure(j)
         j = j + 1
 
      end if
 
   end do
 
   print*, "size(lat_t)=",size(lat_t)  
   print*, "n_raob_t=",n_raob_t
   print*, "size(lat)=",size(lat)  
 
   deallocate(obs_type_t,lat_t,lon_t,stn_elev_t,pressure_t,height_t,prep_qc_mark_t,obs_t,time_t,sid_t, &
         errinv_input_t)
 
   !---------------------------------------
   ! For specific Humidity
   j = 1
   do i = 1, size(obs_q)

      if (obs_type_q(i) == raob_type_tq) then ! if obs is raob

         latq(j)=lat_q(i)
         lonq(j)=lon_q(i)
         stn_elevq(j)=stn_elev_q(i)
         pressureq(j)=pressure_q(i)
         heightq(j)=height_q(i)
         prep_qc_markq(j)=prep_qc_mark_q(i)
         !setup_qc_mark(j)=setup_qc_mark_t(i)
         raob_q_org(j)=obs_q(i)
         ! Convert obstime (hr) relative to the beginning of DA window
         ! to obstime relative to analysis time
!         timeq(j)=time_q(i)-3
         timeq(j)=time_q(i)-offset_argument
         sidq(j)=sid_q(i)
         j = j + 1
 
      end if
 
   end do
 
   deallocate(obs_type_q,lat_q,lon_q,stn_elev_q,pressure_q,height_q,prep_qc_mark_q,obs_q,time_q,sid_q, &
         errinv_input_q)
 
   !---------------------------------------
   ! For wind
   j = 1
   do i = 1, size(obs_u)
 
      if (obs_type_uv(i) == raob_type_uv) then ! if obs is raob
 
         latuv(j)=lat_uv(i)
         lonuv(j)=lon_uv(i)
         stn_elevuv(j)=stn_elev_uv(i)
         pressureuv(j)=pressure_uv(i)
         heightuv(j)=height_uv(i)
         prep_qc_markuv(j)=prep_qc_mark_uv(i)
         !setup_qc_mark(j)=setup_qc_mark_t(i)
         raob_u_org(j)=obs_u(i)
         raob_v_org(j)=obs_v(i)
         ! Convert obstime (hr) relative to the beginning of DA window
         ! to obstime relative to analysis time
!         timeuv(j)=time_uv(i)-3
         timeuv(j)=time_uv(i)-offset_argument
         siduv(j)=sid_uv(i)
         j = j + 1

      end if
 
   end do

   deallocate(obs_type_uv,lat_uv,lon_uv,stn_elev_uv,pressure_uv,height_uv,prep_qc_mark_uv,obs_u,obs_v,time_uv,sid_uv, &
         errinv_input_uv)
 
 !------------------------------------------------------------------
 !  2_2. Extract unique raob station id (remove dups) 
 !       & Sort obs in a descending order according to press
 !------------------------------------------------------------------

   ! Temperature
   obsvar='t'
   !print*, "Before sorting T sid=",sid
   call sorting(obsvar,lat,lon,stn_elev,pressure,height,prep_qc_mark,time,sid,raob_sid_list,sid_plevels_t,sid_plevels_t_500,raob_t)
   !print*, "After sorting T sid=",sid
   !print*, "After sorting, T sid_plevels_t=",sid_plevels_t, "sid_plevels_t_500=",sid_plevels_t_500
 
   ! Specific humidity
   obsvar='q'
   !print*, "Before sorting q sid=",sidq
   call sorting(obsvar,latq,lonq,stn_elevq,pressureq,heightq,prep_qc_markq,timeq,sidq,raob_sid_list_q,sid_plevels_q,sid_plevels_q_500,raob_q_org)
   !print*, "After sorting q sid=",sidq
 
   ! Wind
   obsvar='uv'
   !print*, "Before sorting uv sid=",siduv
   call sorting(obsvar,latuv,lonuv,stn_elevuv,pressureuv,heightuv,prep_qc_markuv,timeuv,siduv,raob_sid_list_uv,sid_plevels_uv,sid_plevels_uv_500,raob_u_org,raob_v_org)
   !print*, "After sorting uv sid=",siduv
   !print*, "After sorting raob_u_org=",raob_u_org
 
 !-------------------------------------------------
 !  2_3. convert utc to local time using longitude
 !       local time = nint(lon/15) + UTC
 !-------------------------------------------------

   allocate(ltime(n_raob_t))
   ltime=r_missing
   do i = 1, size(raob_t)
 
      ! Analysis time
      ana_hh0 = mod(date_time_t, 100)
      ana_dd0 = ( mod(date_time_t, 10000) - ana_hh0 ) / 100 
      ana_mm0 = ( mod(date_time_t, 1000000) - ana_dd0*100 - ana_hh0 ) / 10000 
      ana_yy0 = ( date_time_t - ana_mm0*10000 - ana_dd0*100 - ana_hh0 ) / 1000000

      ! Converting UTC to Local Time (UTC + hr_min)
      ! longitude (degree East) -> East: lon>0, West: lon<0
      if (lon(i)>180) then
         hr_min = (lon(i)-360)/15.
      else
         hr_min = lon(i)/15. 
      end if

      !print*, "LOCAL: L337, lon(i)=",lon(i),", hr=",hr
      ! local time = UTC(ana_time) + nint(lon/15) 
      call convert_localtime_min(ltime_yy,ltime_mm, ltime_dd, ltime_hh, ltime_min, hr_min, ana_yy0, ana_mm0, ana_dd0, ana_hh0,ana_min0)
      !print*, "LOCAL: L340, ltime_hh=",ltime_hh
  
      ltime(i) = ltime_yy*1000000 + ltime_mm*10000 + ltime_dd*100 + ltime_hh 
      print*,"i=",i,", UTC(analysis time)=",date_time_t, ", hr_min=",hr_min," -> local time=",ltime(i)

   end do

!------------------------------------------------------------------
!  3. Extract and assign raob_q (specific humidity) and 
!     raob_w (mixing ratio) for corresponding the same P-levels 
!     of raob temperature (sid)
!     - If there is no same observed pressure level of q compared to T,
!     vertically linear interpolation - log(pressure)
!------------------------------------------------------------------
   allocate(raob_w(n_raob_t))
   raob_w=r_missing
   iloop: do i = 1, size(raob_t)

      jloop: do j = 1, size(raob_q_org)
      !jloop: do j = 1, size(obs_q)

 !       if raob and same sid
         if (sidq(j) == sid(i)) then 
  
            if (pressureq(j) == pressure(i)) then
 
               raob_q(i) = raob_q_org(j)
 
               exit jloop
 
            else ! if there is no same observed pressure level of q compared to T
                 ! then, do linearly vertical log pressure level interpolation
 
               if (j<size(raob_q_org)) then
 
                  if (sidq(j)==sidq(j+1) .and. pressure(i) < pressureq(j) .and. pressure(i) > pressureq(j+1)) then
 
                     slope_q = (raob_q_org(j)-raob_q_org(j+1)) / (log(pressureq(j))-log(pressureq(j+1)))
                     raob_q(i) = raob_q_org(j+1) + (log(pressure(i)) - log(pressureq(j+1))) * slope_q
  
                     !print*, '------------q interpolation to observed p level of T ----------'
                     !print*, 'there is no same observed pressure level of q compared to Temperature******'
                     !print*, 'sidq(j)=',sidq(j), 'pressure(i)=',pressure(i)
                     !print*, 'pressureq(j)=',pressureq(j),"pressureq(j+1)=", pressureq(j+1)
                     !print*, 'raob_q(i)=',raob_q(i), 'obs_q(j), obs_q(j+1)', raob_q_org(j),raob_q_org(j+1)
 
                     exit jloop
 
                  end if
 
               end if
   
            end if
 
         end if

      end do jloop
 
      !------------------------------------------------------------------
      !  4. Convert specific humidity to mixing ratio
      !------------------------------------------------------------------
      raob_w(i) = raob_q(i)/(one-raob_q(i)) 
 
   end do iloop
 
 !------------------------------------------------------------------
 ! 5. Convert temperature to virtual temperature
 ! 6. Convert virtual temperature to virtual potential temperature
 ! 7. Convert temperature to potential temperature
 ! 8. Convert Pressure -> geopotential height 
 !------------------------------------------------------------------
   allocate(height_converted(n_raob_t))
   !allocate(height_converted(n_raob_t), geo_altitude_converted(n_raob_t))
   allocate(raob_tv(n_raob_t), raob_theta(n_raob_t), raob_thetav(n_raob_t))
   height_converted=r_missing
   !geo_altitude_converted=r_missing
   raob_tv=r_missing
   raob_theta=r_missing
   raob_thetav=r_missing
   iloop_convert: do i = 1, size(raob_t)
 
      !------------------------------------------------------------------
      !  5. Convert temperature to virtual temperature
      !  tv(j,i,k) = t(j,i,k) * (one + fv*w(j,i,k))
      !  fv = rv/rd-one    ! used in virtual temperature equation
      !  w = mixing ratio
      !------------------------------------------------------------------
      raob_tv(i) = raob_t(i) * (one + fv * raob_w(i))
      !------------------------------------------------------------------
      !  6. Convert virtual temperature to virtual potential temperature
      !  theta_v = tv * (r1000 / pressure) ** rd_over_cp
      !------------------------------------------------------------------
      raob_thetav(i) = raob_tv(i) * (r1000 / pressure(i)) ** rd_over_cp
      !------------------------------------------------------------------
      !  7. Convert temperature to potential temperature
      !  theta = t * (r1000 / pressure) ** rd_over_cp <-- rd can be used?
      !------------------------------------------------------------------
      raob_theta(i) = raob_t(i) * (r1000 / pressure(i)) ** rd_over_cp

      !------------------------------------------------------------------
      ! 8. Convert pressure -> geopotential height 
      !------------------------------------------------------------------
 
      if (height(i) > 0 .and. height(i) < stn_elev(i)) then
 
         height_converted(i) = height(i) - stn_elev(i) ! Above Ground Level (AGL)
         !call convert_geop_height_to_geom_altitude(lat(i), height_converted(i), geom_alti)
         !geo_altitude_converted(i) = geom_alti ! Above Ground Level (AGL)
         print*, "-----------------------------------"
         print*, "*** height is lower than stn_elev : height(i) < stn_elev(i) ***"
         print*, "i=",i, "height(i)=",height(i),"stn_elev(i)=",stn_elev(i)
         print*, "        pressure(i)=",pressure(i), "height_converted(i)=",height_converted(i)
         print*, "-----------------------------------"

      else if (height(i) == stn_elev(i)) then
 
         height_converted(i) = zero ! Above Ground Level (AGL) 
         !call convert_geop_height_to_geom_altitude(lat(i), height_converted(i), geom_alti)
         !geo_altitude_converted(i) = geom_alti ! Above Ground Level (AGL)

      else if (sid(i).ne.sid(i-1)) then ! first level for sid

         if (height(i) == r_missing) then ! if first level height is missing

            print*, "-----------------------------------"
            print*, "---i=",i,",missing height(first level), stn_elev(i)=",stn_elev(i)
            print*, "---height(i) is missing--- sid:",sid(i), "p=",pressure(i)
            height_converted(i) = r_missing
            !geo_altitude_converted(i) = r_missing ! Above Ground Level (AGL)
            print*, "-----------------------------------"

         ! height > stn_elev but the lowest level is above sfc level
         else if (height(i) > stn_elev(i)) then
 
            print*, "-----------------------------------"
            print*, "---i=",i,"--else (height > stn_elev but the lowest level is above sfc level)"
            print*, "---lowest level height(i) >stn elev--- sid:",sid(i), "p=",pressure(i)
            height_converted(i) = height(i) - stn_elev(i) ! Above Ground Level (AGL)
            !call convert_geop_height_to_geom_altitude(lat(i), height_converted(i), geom_alti)
            !geo_altitude_converted(i) = geom_alti ! Above Ground Level (AGL)
            print*, "    height(i)=",height(i),"stn_elev(i)=",stn_elev(i)
            print*, "    height_converted(i)=",height_converted(i)
            !print*, "    geo_altitude_converted(i) =",geo_altitude_converted(i)
            print*, "-----------------------------------"

         end if

      else if (sid(i) == sid(i-1)) then  ! For most cases (from 2nd level)

         if (height(i) == stn_elev(i)) then ! if stn_elev is found above the lowest level

            height_converted(i) = zero ! Above Ground Level (AGL)
            !call convert_geop_height_to_geom_altitude(lat(i), height_converted(i), geom_alti)
            !geo_altitude_converted(i) = geom_alti ! Above Ground Level (AGL)

         else if (height_converted(i-1) == r_missing) then ! if the previous level height is missing

            print*, "--------------------------------"
            print*, "---- previous level height is missing"
            print*, "---- stn_elev(i)=",stn_elev(i)
            print*, "---- sid:",sid(i), "p(i)=",pressure(i),"p(i-1)=",pressure(i-1)
            print*, "--------------------------------"
            height_converted(i) = r_missing ! ! Above Ground Level (AGL)
            !geo_altitude_converted(i) = r_missing ! Above Ground Level (AGL)
            print*, "-----------------------------------"

         else
 
            print*, "-----------------------------------"
            print*, "-------sid(i) == sid(i-1), calculate h from two p levels"
            print*, "i=",i, "sid(i)=",sid(i),"sid(i-1)=",sid(i-1)
 
            h = rdog * half * (raob_tv(i-1) + raob_tv(i))
            dz = h * log(pressure(i-1)/pressure(i))
            height_converted(i) = height_converted(i-1) + dz
            print*, "p(i)=",pressure(i),",height_converted(i)=",height_converted(i)
            print*, "p(i-1)=",pressure(i-1),",height_converted(i-1)=",height_converted(i-1)
            print*, "dz=",dz
 
            !call convert_geop_height_to_geom_altitude(lat(i), height_converted(i), geom_alti)
            !geo_altitude_converted(i) = geom_alti ! Above Ground Level (AGL)
            !print*, "geo_altitude_converted(i) =",geo_altitude_converted(i)

         end if

      end if
 
   end do iloop_convert

   !------------------------------------------------------------
   ! 9. Vertically linearly interpolate observed P levels of wind 
   ! to observed P levels of T
   !------------------------------------------------------------
   allocate(raob_u(n_raob_t), raob_v(n_raob_t), raob_ri(n_raob_t), raob_pblh_ri(n_raob_t), &
        raob_u_sfc(n_raob_t), raob_v_sfc(n_raob_t), raob_thetav_sfc(n_raob_t),&
        height_converted_sfc(n_raob_t))
        !geo_altitude_converted_sfc(n_raob_t))
   allocate(raob_pblh_thetav(n_raob_t))
   allocate(raob_pblh_q(n_raob_t))
   raob_u=r_missing
   raob_v=r_missing
   raob_ri=r_missing
   raob_pblh_ri=r_missing
   raob_u_sfc=r_missing
   raob_v_sfc=r_missing
   raob_thetav_sfc=r_missing
   height_converted_sfc=r_missing
   !geo_altitude_converted_sfc=r_missing
   raob_pblh_thetav=r_missing
   raob_pblh_q=r_missing
   iloop_t_interp: do i = 1, size(raob_t)
 
      jloop_uv_interp: do j = 1, size(raob_u_org)

  !      if raob and same sid
         if (siduv(j) == sid(i)) then
 
            if (pressureuv(j) == pressure(i)) then
 
               raob_u(i) = raob_u_org(j)
               raob_v(i) = raob_v_org(j)
               !print*, '--------------------------------------------------------'
               !print*, 'UV: j=,',j,',siduv(j)=',siduv(j), '** same pressure(i)=',pressure(i)
               !print*, 'raob_u_org(j)', raob_u_org(j)
               !print*, 'raob_v_org(j)', raob_v_org(j)
               !print*, '--------------------------------------------------------'
               exit jloop_uv_interp
 
            else

               if (j<size(raob_u_org)) then
 
                  if (siduv(j)==siduv(j+1) .and. pressure(i) < pressureuv(j) .and. pressure(i) > pressureuv(j+1)) then
                  
                     slope_u = (raob_u_org(j)-raob_u_org(j+1)) / (log(pressureuv(j))-log(pressureuv(j+1)))
                     raob_u(i) = raob_u_org(j+1) + (log(pressure(i)) - log(pressureuv(j+1))) * slope_u
                       
                     slope_v = (raob_v_org(j)-raob_v_org(j+1)) / (log(pressureuv(j))-log(pressureuv(j+1)))
                     raob_v(i) = raob_v_org(j+1) + (log(pressure(i)) - log(pressureuv(j+1))) * slope_v
                     !print*, 'siduv(j)=',siduv(j), 'pressure(i)=',pressure(i)
                     !print*, 'pressureuv(j)=',pressureuv(j),"pressureuv(j+1)=", pressureuv(j+1)
                     !print*, 'raob_u(i)=',raob_u(i), 'raob_u_org(j), raob_u_org(j+1)', raob_u_org(j),raob_u_org(j+1)
                       
                     exit jloop_uv_interp
 
                  end if

               end if

            end if

         end if

      end do jloop_uv_interp

   end do iloop_t_interp

   !------------------------------------------------------------------
   !  10. Assign sfc level values to height, u, v, theta_v
   !------------------------------------------------------------------

   iloop_sid_list: do i = 1, size(raob_sid_list)

      height_sfc=0
      u_sfc=0
      v_sfc=0
      thetav_sfc=0
      jloop_t: do j = 1, size(raob_t)

         if (sid(j) == raob_sid_list(i)) then

            if (height(j) == stn_elev(j)) then 

               height_sfc = height_converted(j) 
               u_sfc = raob_u(j) 
               v_sfc = raob_v(j) 
               thetav_sfc = raob_thetav(j) 

            else if (height(j) > 0 .and. height(j) < stn_elev(j)) then ! ?

               height_sfc = height_converted(j) 
               u_sfc = raob_u(j) 
               v_sfc = raob_v(j) 
               thetav_sfc = raob_thetav(j) 

            else if (j > 1 .and. sid(j).ne.sid(j-1)) then ! first level for sid

               if (height(j) > stn_elev(j)) then

                  height_sfc = height_converted(j) 
                  u_sfc = raob_u(j) 
                  v_sfc = raob_v(j) 
                  thetav_sfc = raob_thetav(j) 

               end if

            end if

            height_converted_sfc(j) = height_sfc
            raob_u_sfc(j) = u_sfc
            raob_v_sfc(j) = v_sfc
            raob_thetav_sfc(j) = thetav_sfc

            if (thetav_sfc == 0) then !?
               print*, "thetav_sfc == 0, height(j)=",height(j), ",stn_elev(j)=",stn_elev(j)
            end if

         end if

      end do jloop_t
 
   end do iloop_sid_list
   !------------------------------------------------------------------
   !  11. Calculate Bulk Richardson number
   !------------------------------------------------------------------

   !  Ri_b =     g * (theta_v_z - theta_v_sfc) * (z - z_sfc)
   !        -------------------------------------------------------
   !        (theta_v_sfc) * ( (u_z - u_sfc)**2 + (v_z - v_sfc)**2 )

   do i = 1, size(raob_t)

      numer = (raob_thetav(i) - raob_thetav_sfc(i)) * (height_converted(i))
      !numer = (raob_thetav(i) - raob_thetav_sfc(i)) * (geo_altitude_converted(i))
      !numer = (raob_thetav(i) - raob_thetav_sfc(i)) * (geo_altitude_converted(i)-geo_altitude_converted_sfc(i))
      ! i) consider sfc wind 
      !denom = (raob_u(i) - raob_u_sfc(i))**2 + (raob_v(i) - raob_v_sfc(i))**2
      ! ii) consider sfc wind as zero 
      denom = raob_u(i)**2 + raob_v(i)**2
      !raob_ri(i) = grav / raob_thetav_sfc(i) * numer / denom
      if (numer == 0) then ! at sfc level -> RiB = zero/zero ??????????
         raob_ri(i) = zero/zero
      else
         raob_ri(i) = grav / raob_thetav_sfc(i) * numer / denom
      end if

   end do

   !------------------------------------------------------------------
   !  12. Calculate thetav-gradient and q-gradient PBLH 
   !------------------------------------------------------------------
   !raob_pblh_thetav

   allocate(pblh_q(size(raob_sid_list)))
   allocate(pblh_q_lowlevel(size(raob_sid_list)))
   allocate(lat_pblh_q(size(raob_sid_list)))
   allocate(lon_pblh_q(size(raob_sid_list)))
   allocate(time_pblh_q(size(raob_sid_list)))
   allocate(ltime_pblh_q(size(raob_sid_list)))

   allocate(pblh_thetav(size(raob_sid_list)))
   allocate(pblh_thetav_lowlevel(size(raob_sid_list)))
   allocate(lat_pblh_thetav(size(raob_sid_list)))
   allocate(lon_pblh_thetav(size(raob_sid_list)))
   allocate(time_pblh_thetav(size(raob_sid_list)))
   allocate(ltime_pblh_thetav(size(raob_sid_list)))

   pblh_q=r_missing
   pblh_q_lowlevel=r_missing
   lat_pblh_q=r_missing
   lon_pblh_q=r_missing
   time_pblh_q=r_missing
   ltime_pblh_q=r_missing

   pblh_thetav=r_missing
   pblh_thetav_lowlevel=r_missing
   lat_pblh_thetav=r_missing
   lon_pblh_thetav=r_missing
   time_pblh_thetav=r_missing
   ltime_pblh_thetav=r_missing

   ! thetav-gradient-based PBLH
   !
   ! 1) averaged value between two p levels:
   !    pblh_thetav, lat_pblh_thetav, lon_pblh_thetav, time_pblh_thetav
   ! 2) value at the lower level among two p levels:
   !    ltime_pblh_thetav
   !    pblh_thetav_lowlevel 

   do i = 1, size(raob_sid_list)

      maxdthvdz = zero
      mindqdz = one
      do j = 1, size(raob_t)

         if (sid(j) == raob_sid_list(i)) then

            !if ( geo_altitude_converted(j) <= 5000 ) then
            if ( pressure(j) >= 500 ) then

               if (j > 1 .and. sid(j) == sid(J-1)) then

                  dthvdz = (raob_thetav(j) - raob_thetav(j-1)) / ( height_converted(j) - height_converted(j-1)) 
                  dqdz = (raob_q(j) - raob_q(j-1)) / ( height_converted(j) - height_converted(j-1)) 
                  
                  ! thetav-gradient 
                  if (dthvdz > maxdthvdz) then
   !              if (height1(j-1) > stn_elev1(j-1)) then ! if higher than 2nd level,
                     maxdthvdz = dthvdz
                     pblh_thetav(i) = 0.5 * ( height_converted(j-1) + height_converted(j))
                     pblh_thetav_lowlevel(i) = height_converted(j-1)

                     ! lat, lon, and time
                     lat_pblh_thetav(i) = 0.5 * (lat(j-1) + lat(j))
                     lon_pblh_thetav(i) = 0.5 * (lon(j-1) + lon(j))
                     time_pblh_thetav(i) = 0.5 * (time(j-1) + time(j))
                     ltime_pblh_thetav(i) = ltime(j-1) ! local time is integer
                  
                  end if
                  
                  ! q-gradient  
                  ! only consider q gradient based pblh > 50 m
                  !if (dqdz < mindqdz .and.  0.5 * ( height_converted(j-1) + height_converted(j)) > 50.0) then
                  if (dqdz < mindqdz) then
   !              if (height1(j-1) > stn_elev1(j-1)) then ! if higher than 2nd level,
                     mindqdz = dqdz
                     pblh_q(i) = 0.5 * ( height_converted(j-1) + height_converted(j))
                     pblh_q_lowlevel(i) = height_converted(j-1)

                     ! lat, lon, and time
                     lat_pblh_q(i) = 0.5 * (lat(j-1) + lat(j))
                     lon_pblh_q(i) = 0.5 * (lon(j-1) + lon(j))
                     time_pblh_q(i) = 0.5 * (time(j-1) + time(j))
                     ltime_pblh_q(i) = ltime(j-1) ! local time is integer

                  end if

               end if

            end if

         end if

      end do 

   end do 

   do i = 1, size(raob_sid_list)

      do j = 1, size(raob_t)

         if (sid(j) == raob_sid_list(i)) then

            raob_pblh_thetav(j) = pblh_thetav(i)
            raob_pblh_q(j) = pblh_q(i)

         end if

      end do 

   end do

   !------------------------------------------------------------------
   !  13. Calculate RiB-based PBLH 
   !     ! QC Flag, local time, 
   !------------------------------------------------------------------
   !raob_pblh_ri
   !raob_pblh_q
   !raob_pblh_thetav
   !pblh_qc_mark

   allocate(pblh_ri(size(raob_sid_list)))
   allocate(lat_pblh_ri(size(raob_sid_list)))
   allocate(lon_pblh_ri(size(raob_sid_list)))
   allocate(time_pblh_ri(size(raob_sid_list)))
   allocate(ltime_pblh_ri(size(raob_sid_list)))

   allocate(qc_mark_pblh(size(raob_sid_list)))
   allocate(stn_elev_pblh(size(raob_sid_list)))

   pblh_ri=r_missing
   lat_pblh_ri=r_missing
   lon_pblh_ri=r_missing
   time_pblh_ri=r_missing
   ltime_pblh_ri=r_missing
   qc_mark_pblh=r_missing
   stn_elev_pblh=r_missing

   ! Ri-based PBLH
   ! 1) vertically linearly interpolated value between two p levels:
   !    pblh_ri, lat, lon, time
   ! 2) value at the lower level among two p levels:
   !    ltime_pblh_ri

   iloop_sid: do i = 1, size(raob_sid_list)

      ! QC Flag

      ! if less than 10 p-levels below 500 hPa
      ! -> qc = 5, otherwise qc = 0
      if (sid_plevels_t_500(i) < 10) then
         qc_mark_pblh(i) = 5
      else
         qc_mark_pblh(i) = 0
      end if

      jloop_tlist: do j = 1, size(raob_t)

         if (sid(j) == raob_sid_list(i)) then

            stn_elev_pblh(i) = stn_elev(j)

            if (j > 1 .and. sid(j) == sid(j-1)) then

               if (raob_ri(j) >= tcri_crit) then

                  ! PBLH
                  ! for 2nd obs level, pblh_ri should be assigned as height at 2nd level (AGL)
                  ! QC_flag = 50
                  if (raob_thetav(j-1) == raob_thetav_sfc(j-1)) then ! if 2nd level, pblh_ri should be 2nd level
                     pblh_ri(i) = height_converted(j)
                     lat_pblh_ri(i) = lat(j)
                     lon_pblh_ri(i) = lon(j)
                     time_pblh_ri(i) = time(j)
                     ltime_pblh_ri(i) = ltime(j)

                     if (pblh_ri(i) > 6000) then
                        pblh_ri(i) = r_missing
                        qc_mark_pblh(i) = 99 ! pblh>6000m : qc=99
                     else
                        qc_mark_pblh(i) = 50
                     end if

                     exit jloop_tlist

                  else

                     pblh_ri(i) = height_converted(j-1) + (tcri_crit - raob_ri(j-1)) / &
                           (raob_ri(j) - raob_ri(j-1)) * (height_converted(j) - height_converted(j-1))
                     ! lat, lon, and time
                     lat_pblh_ri(i) = lat(j-1) + (tcri_crit - raob_ri(j-1)) / &
                           (raob_ri(j) - raob_ri(j-1)) * (lat(j) - lat(j-1))
                     lon_pblh_ri(i) = lon(j-1) + (tcri_crit - raob_ri(j-1)) / &
                           (raob_ri(j) - raob_ri(j-1)) * (lon(j) - lon(j-1))
                     time_pblh_ri(i) = time(j-1) + (tcri_crit - raob_ri(j-1)) / &
                           (raob_ri(j) - raob_ri(j-1)) * (time(j) - time(j-1))

                     ! local time: lower level value
                     ltime_pblh_ri(i) = ltime(j-1) ! local time is integer

                     if (pblh_ri(i) > 6000) then
                        pblh_ri(i) = r_missing
                        qc_mark_pblh(i) = 99 ! pblh>6000m : qc=99
                     end if

                     exit jloop_tlist

                  end if

               end if

            end if

         end if

      end do jloop_tlist

   end do iloop_sid


   do i = 1, size(raob_sid_list)

      do j = 1, size(raob_t)

         if (sid(j) == raob_sid_list(i)) then

            raob_pblh_ri(j) = pblh_ri(i)

         end if

      end do 

   end do

   !------------------------------------------------------------------
   !  14. Save obstime for YYYY, MM, DD, HH, and Min format
   !------------------------------------------------------------------

   allocate(yyyy_pblh_ri(size(raob_sid_list)))
   allocate(mm_pblh_ri(size(raob_sid_list)))
   allocate(dd_pblh_ri(size(raob_sid_list)))
   allocate(hh_pblh_ri(size(raob_sid_list)))
   allocate(min_pblh_ri(size(raob_sid_list)))

   yyyy_pblh_ri = r_missing
   mm_pblh_ri = r_missing
   dd_pblh_ri = r_missing
   hh_pblh_ri = r_missing
   min_pblh_ri = r_missing

   do i = 1, size(raob_sid_list)

      ! time_pblh_ri = obstime relative to analysis time

      ! Save obstime depending on yyyy, mm, dd, hh, and min
      !ana_time(i) + time_pblh_ri(i) -> obstime relative to analysis time 
      !convert_obstime_min(obstime_mm, obstime_dd, obstime_hh,obstime_min, interval, ana_mm0, ana_dd0, ana_hh0)

      if (time_pblh_ri(i) .ne. r_missing .and. .not.isnan(time_pblh_ri(i))) then
         call convert_obstime_min(yyyy_pblh_ri(i),mm_pblh_ri(i),dd_pblh_ri(i),hh_pblh_ri(i),min_pblh_ri(i), time_pblh_ri(i), ana_yy0, ana_mm0, ana_dd0, ana_hh0)
      end if

      !print*, "======= i = ",i, ",raob_sid_list(q)=",raob_sid_list(q), ", ana_time(q)=",ana_time(q)
      !print*, "time_pblh_ri(i)=",time_pblh_ri(i)
      !print*, "mm_pblh_ri(i)=",mm_pblh_ri(i)
      !print*, "dd_pblh_ri(i)=",dd_pblh_ri(i)
      !print*, "hh_pblh_ri(i)=",hh_pblh_ri(i)
      !print*, "min_pblh_ri(i)=",min_pblh_ri(i)

   end do
   !------------------------------------------------------------------
   !print*, "RESULTS"
   !------------------------------------------------------------------
   !do i = 1, size(raob_sid_list)
      !print*, "====== LOCAL: i=",i, ", SID=",raob_sid_list(i)
      !print*, "ltime_pblh_ri=",ltime_pblh_ri(i)
      !print*, "pblh_ri=",pblh_ri(i)
      !print*, "pblh_thetav=",pblh_thetav(i)
      !print*, "pblh_q=",pblh_q(i)
      !print*, "qc_mark_pblh=",qc_mark_pblh(i)
      !print*, "sid_plevels_t=",sid_plevels_t(i)
      !print*, "sid_plevels_t_500=",sid_plevels_t_500(i)
      !print*, "lat_pblh_ri=",lat_pblh_ri(i)
      !print*, "lon_pblh_ri=",lon_pblh_ri(i)
   !end do

   do i = 1, size(raob_t)

      print*, "FINAL: i=",i, ",** sid=",sid(i)
      print*, "            ,** p=",pressure(i)," hPa"
      print*, "              height=",height(i), "stn_elev=",stn_elev(i) ! [m]
      print*, "              height_converted=", height_converted(i)
      !print*, "              altitude_converted=", geo_altitude_converted(i)
      !print*, "              geo_altitude_converted_sfc=", geo_altitude_converted_sfc(i)
      print*, "              height_converted_sfc=", height_converted_sfc(i)
      print*, "                raob_u_sfc=",raob_u_sfc(i)
      print*, "                raob_v_sfc=",raob_v_sfc(i)
      print*, "                raob_u=",raob_u(i)
      print*, "                raob_v=",raob_v(i)
      print*, "                raob_thetav_sfc=",raob_thetav_sfc(i)
      print*, "--------------------------------------------------"
      print*, "                raob_thetav=",raob_thetav(i)
      print*, "                raob_theta=",raob_theta(i)
      print*, "                raob_ri=",raob_ri(i)
      print*, "                raob_q=",raob_q(i)
      print*, "--------------------------------------------------"
      print*, "                ltime=",ltime(i)
      print*, "--------------------------------------------------"
      print*, "                raob_pblh_ri=",raob_pblh_ri(i)
      print*, "                raob_pblh_thetav=",raob_pblh_thetav(i)
      print*, "                raob_pblh_q=",raob_pblh_q(i)
      print*, "--------------------------------------------------"
!     print*, "                raob_w=",raob_w(i)
!     print*, "                raob_t=",raob_t(i)
!     print*, "                raob_tv=",raob_tv(i)
!     print*, "                raob_theta=",raob_theta(i)
!     print*, "                raob_thetav=",raob_thetav(i)

   end do

   !------------------------------------------------------------------
   !  15. Write nc files
   !------------------------------------------------------------------

   !  Write into netcdf files for each analysis window

   write(filename, '(A10,I10,A4)') 'raob_pblh.', date_time_t, '.nc4'
   !  ierr2 =  NF90_CREATE(trim(filename2(i)), NF90_CLOBBER, ncid2(i)) ! overwrite this file if it already exists
   ierr =  NF90_CREATE(trim(filename), NF90_NETCDF4, ncid)
   if (ierr /= nf90_noerr) call handle_err(ierr,"create")
   !  ierr2 = nf90_open(trim(filename2(i)),nf90_write, ncid2(i))
   !  if (ierr2 /= nf90_noerr) call handle_err(ierr2,"open")
   print*, "filename = ", filename

   !-----------------------
   ! Define the dimensions.
   ierr = nf90_def_dim(ncid, 'Time', 1, dimid3)
   ierr = nf90_def_dim(ncid, 'nobs', n_raob_t, dimid1)
   ierr = nf90_def_dim(ncid, 'nsid', size(raob_sid_list), dimid2)
   ierr = nf90_def_dim(ncid, 'Station_ID_maxstrlen', 8, dimid4)

   !-----------------------
   ! Define variables
   ierr = nf90_def_var(ncid, 'Ana_Time', NF90_INT64, dimid3, varid1)
   !ierr = nf90_def_var(ncid, 'Loc_Time_ri', NF90_INT64, dimid2, varid2)
   ierr = nf90_def_var(ncid, 'lat', NF90_DOUBLE, dimid2, varid3)
   ierr = nf90_def_var(ncid, 'lon', NF90_DOUBLE, dimid2, varid4)
   ierr = nf90_def_var(ncid, 'Year', NF90_DOUBLE, dimid2, varid5)
   ierr = nf90_def_var(ncid, 'Month', NF90_DOUBLE, dimid2, varid6)
   ierr = nf90_def_var(ncid, 'Date', NF90_DOUBLE, dimid2, varid7)
   ierr = nf90_def_var(ncid, 'Hour', NF90_DOUBLE, dimid2, varid8)
   ierr = nf90_def_var(ncid, 'Minute', NF90_DOUBLE, dimid2, varid9)
   ierr = nf90_def_var(ncid, 'PBL_Height_ri', NF90_DOUBLE, dimid2, varid10)
   ierr = nf90_def_var(ncid, 'PBL_Height_thetav', NF90_DOUBLE, dimid2, varid11)
   ierr = nf90_def_var(ncid, 'PBL_Height_thetav_lowlevel', NF90_DOUBLE, dimid2, varid12)
   ierr = nf90_def_var(ncid, 'PBL_Height_q', NF90_DOUBLE, dimid2, varid13)
   ierr = nf90_def_var(ncid, 'PBL_Height_q_lowlevel', NF90_DOUBLE, dimid2, varid14)
   ierr = nf90_def_var(ncid, 'Sid', NF90_CHAR, (/dimid4, dimid2/), varid15)
   ierr = nf90_def_var(ncid, 'Loc_Time_ri', NF90_INT64, dimid2, varid16)
   ierr = nf90_def_var(ncid, 'Loc_Time_thetav', NF90_INT64, dimid2, varid17)
   ierr = nf90_def_var(ncid, 'Loc_Time_q', NF90_INT64, dimid2, varid18)
   ierr = nf90_def_var(ncid, 'QC_Mark_PBLH', NF90_INT64, dimid2, varid19)
   ierr = nf90_def_var(ncid, 'Number_of_PLevels', NF90_DOUBLE, dimid2, varid20)
   ierr = nf90_def_var(ncid, 'Number_of_PLevels_500', NF90_DOUBLE, dimid2, varid21)
   ierr = nf90_def_var(ncid, 'Stn_Elev', NF90_DOUBLE, dimid2, varid22)

   ! Vertical Profiles (dimid1 = nobs)
   ierr = nf90_def_var(ncid, 'All_Pressure', NF90_DOUBLE, dimid1, varid30)
   ierr = nf90_def_var(ncid, 'All_Sid', NF90_CHAR, (/dimid4, dimid1/), varid31)
   ierr = nf90_def_var(ncid, 'All_Altitude_AGL', NF90_DOUBLE, dimid1, varid32)
   ierr = nf90_def_var(ncid, 'All_Stn_Elev', NF90_DOUBLE, dimid1, varid33)
   ierr = nf90_def_var(ncid, 'All_Obs_t', NF90_DOUBLE, dimid1, varid34)
   ierr = nf90_def_var(ncid, 'All_Obs_q', NF90_DOUBLE, dimid1, varid35)
   ierr = nf90_def_var(ncid, 'All_Obs_w', NF90_DOUBLE, dimid1, varid36)
   ierr = nf90_def_var(ncid, 'All_Obs_theta', NF90_DOUBLE, dimid1, varid37)
   ierr = nf90_def_var(ncid, 'All_Obs_thetav', NF90_DOUBLE, dimid1, varid38)
   ierr = nf90_def_var(ncid, 'All_Obs_u', NF90_DOUBLE, dimid1, varid39)
   ierr = nf90_def_var(ncid, 'All_Obs_v', NF90_DOUBLE, dimid1, varid40)
   ierr = nf90_def_var(ncid, 'All_Ri', NF90_DOUBLE, dimid1, varid41)
   ierr = nf90_def_var(ncid, 'All_Loc_Time', NF90_INT64, dimid1, varid42)
   ierr = nf90_def_var(ncid, 'All_Obs_Time', NF90_DOUBLE, dimid1, varid43)
   ierr = nf90_def_var(ncid, 'All_QC_Mark_t', NF90_DOUBLE, dimid1, varid44)
   ierr = nf90_def_var(ncid, 'All_PBLH_ri', NF90_DOUBLE, dimid1, varid45)
   ierr = nf90_def_var(ncid, 'All_PBLH_thetav', NF90_DOUBLE, dimid1, varid46)
   ierr = nf90_def_var(ncid, 'All_PBLH_q', NF90_DOUBLE, dimid1, varid47)
   ierr = nf90_def_var(ncid, 'All_Lat', NF90_DOUBLE, dimid1, varid48)
   ierr = nf90_def_var(ncid, 'All_Lon', NF90_DOUBLE, dimid1, varid49)
   !ierr = nf90_def_var(ncid, 'SurfaceElevation', NF90_DOUBLE, dimid1, varid13)
   !ierr = nf90_def_var(ncid, 'Land_Water_Mask', NF90_DOUBLE, dimid1, varid14)
   !ierr = nf90_def_var(ncid, 'Lidar_Data_Altitude', NF90_DOUBLE, dimid2, varid11)
   !ierr = nf90_def_var(ncid, 'Total_Attenuated_Backscatter_532', NF90_DOUBLE, (/ dimid1,dimid2 /), varid12)

   ierr = nf90_put_att(ncid, varid1, 'units', 'YYYYMMDDHH')
   ierr = nf90_put_att(ncid, varid3, 'units', 'degrees_north')
   ierr = nf90_put_att(ncid, varid4, 'units', 'degrees_east')
   ierr = nf90_put_att(ncid, varid10, 'units', 'meters')
   ierr = nf90_put_att(ncid, varid11, 'units', 'meters')
   ierr = nf90_put_att(ncid, varid12, 'units', 'meters')
   ierr = nf90_put_att(ncid, varid13, 'units', 'meters')
   ierr = nf90_put_att(ncid, varid14, 'units', 'meters')
   ierr = nf90_put_att(ncid, varid16, 'units', 'YYYYMMDDHH')
   ierr = nf90_put_att(ncid, varid17, 'units', 'YYYYMMDDHH')
   ierr = nf90_put_att(ncid, varid18, 'units', 'YYYYMMDDHH')
   ierr = nf90_put_att(ncid, varid22, 'units', 'meters')

   ! Vertical Profiles
   ierr = nf90_put_att(ncid, varid30, 'units', 'hPa')
   ierr = nf90_put_att(ncid, varid32, 'units', 'meters')
   ierr = nf90_put_att(ncid, varid33, 'units', 'meters')
   ierr = nf90_put_att(ncid, varid34, 'units', 'K')
   ierr = nf90_put_att(ncid, varid35, 'units', 'kg/kg')
   ierr = nf90_put_att(ncid, varid36, 'units', 'kg/kg')
   ierr = nf90_put_att(ncid, varid37, 'units', 'K')
   ierr = nf90_put_att(ncid, varid38, 'units', 'K')
   ierr = nf90_put_att(ncid, varid39, 'units', 'm/s')
   ierr = nf90_put_att(ncid, varid40, 'units', 'm/s')
   ierr = nf90_put_att(ncid, varid42, 'units', 'YYYYMMDDHH')
   ierr = nf90_put_att(ncid, varid43, 'units', 'hours relative to analysis time')
   ierr = nf90_put_att(ncid, varid45, 'units', 'meters')
   ierr = nf90_put_att(ncid, varid46, 'units', 'meters')
   ierr = nf90_put_att(ncid, varid47, 'units', 'meters')
   ierr = nf90_put_att(ncid, varid48, 'units', 'degrees_north') 
   ierr = nf90_put_att(ncid, varid49, 'units', 'degrees_east')  

   ierr = nf90_enddef(ncid)

   !-----------------
   ! Put variables
   ierr = nf90_put_var(ncid, varid1, ana_time)
   ierr = nf90_put_var(ncid, varid3, lat_pblh_ri)
   ierr = nf90_put_var(ncid, varid4, lon_pblh_ri)
   ierr = nf90_put_var(ncid, varid5, yyyy_pblh_ri)
   ierr = nf90_put_var(ncid, varid6, mm_pblh_ri)
   ierr = nf90_put_var(ncid, varid7, dd_pblh_ri)
   ierr = nf90_put_var(ncid, varid8, hh_pblh_ri)
   ierr = nf90_put_var(ncid, varid9, min_pblh_ri)
   ierr = nf90_put_var(ncid, varid10, pblh_ri)
   ierr = nf90_put_var(ncid, varid11, pblh_thetav)
   ierr = nf90_put_var(ncid, varid12, pblh_thetav_lowlevel)
   ierr = nf90_put_var(ncid, varid13, pblh_q)
   ierr = nf90_put_var(ncid, varid14, pblh_q_lowlevel)
   ierr = nf90_put_var(ncid, varid15, raob_sid_list)
   ierr = nf90_put_var(ncid, varid16, ltime_pblh_ri)
   ierr = nf90_put_var(ncid, varid17, ltime_pblh_thetav)
   ierr = nf90_put_var(ncid, varid18, ltime_pblh_q)
   ierr = nf90_put_var(ncid, varid19, qc_mark_pblh)
   ierr = nf90_put_var(ncid, varid20, sid_plevels_t)
   ierr = nf90_put_var(ncid, varid21, sid_plevels_t_500)
   ierr = nf90_put_var(ncid, varid22, stn_elev_pblh)
   ! Vertical Profile
   ierr = nf90_put_var(ncid, varid30, pressure)
   ierr = nf90_put_var(ncid, varid31, sid)
   ierr = nf90_put_var(ncid, varid32, height_converted)
   ierr = nf90_put_var(ncid, varid33, stn_elev)
   ierr = nf90_put_var(ncid, varid34, raob_t)
   ierr = nf90_put_var(ncid, varid35, raob_q)
   ierr = nf90_put_var(ncid, varid36, raob_w)
   ierr = nf90_put_var(ncid, varid37, raob_theta)
   ierr = nf90_put_var(ncid, varid38, raob_thetav)
   ierr = nf90_put_var(ncid, varid39, raob_u)
   ierr = nf90_put_var(ncid, varid40, raob_v)
   ierr = nf90_put_var(ncid, varid41, raob_ri)
   ierr = nf90_put_var(ncid, varid42, ltime)
   ierr = nf90_put_var(ncid, varid43, time)
   ierr = nf90_put_var(ncid, varid44, prep_qc_mark)
   ierr = nf90_put_var(ncid, varid45, raob_pblh_ri)
   ierr = nf90_put_var(ncid, varid46, raob_pblh_thetav)
   ierr = nf90_put_var(ncid, varid47, raob_pblh_q)
   ierr = nf90_put_var(ncid, varid48, lat)
   ierr = nf90_put_var(ncid, varid49, lon)
   ierr = NF90_CLOSE(ncid)
   !--------------
   ! Deallocate
   !--------------
   deallocate(lat,lon,stn_elev,pressure,height,prep_qc_mark,raob_t,raob_q,time,sid)
   deallocate(latq,lonq,stn_elevq,pressureq,heightq,prep_qc_markq,raob_q_org,timeq,sidq)
   deallocate(latuv,lonuv,stn_elevuv,pressureuv,heightuv,prep_qc_markuv,raob_u_org,raob_v_org,timeuv,siduv)
   deallocate(ltime,raob_w)
   deallocate(height_converted)
   deallocate(raob_tv,raob_theta,raob_thetav)
   deallocate(raob_u,raob_v,raob_ri,raob_pblh_ri,raob_u_sfc,raob_v_sfc,raob_thetav_sfc)
   deallocate(height_converted_sfc)
   deallocate(raob_pblh_thetav, raob_pblh_q)
   deallocate(pblh_q,pblh_q_lowlevel,lat_pblh_q,lon_pblh_q,time_pblh_q,ltime_pblh_q)
   deallocate(pblh_thetav,pblh_thetav_lowlevel,lat_pblh_thetav,lon_pblh_thetav,time_pblh_thetav,ltime_pblh_thetav)
   deallocate(pblh_ri,lat_pblh_ri,lon_pblh_ri,time_pblh_ri,ltime_pblh_ri,qc_mark_pblh)
   deallocate(yyyy_pblh_ri,mm_pblh_ri,dd_pblh_ri,hh_pblh_ri,min_pblh_ri)
   deallocate(stn_elev_pblh)

   contains

!*****************************************************************
   subroutine sorting(obsvar,lat,lon,stn_elev,pressure,height,prep_qc_mark,time,sid,raob_sid_list,sid_plevels,sid_plevels_500,raob,raob_v)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    sorting    sort radiosonde observations for descending pressure (bottom to top) 
!   prgmmr: eyang        org: GMAO/NASA               date: 2023-06-26
!
! program history log:
!   2023-06-26  eyang    
!
!   input argument list:
!     obsvar    - observation variable (t, q, uv)
!     lat       - observation latitude
!     lon       - observation longitude
!     stn_elev  - station elevation [m]
!     pressure  - observation pressure
!     height    - observation height
!     prep_qc_mark - qc flag
!     raob      - observation 
!     time      - observation time
!     sid       - station id
!     raob_v    - v observation (optional)
!
!   output argument list:
!     lat       - observation latitude
!     lon       - observation longitude
!     stn_elev  - station elevation [m]
!     pressure  - observation pressure
!     height    - observation height
!     prep_qc_mark - qc flag
!     raob      - observation 
!     time      - observation time
!     sid       - station id
!     raob_v    - v observation (optional)
!     raob_sid_list   -  station id list
!     sid_plevels     -  the number of p-levels
!     sid_plevels_500 -  the number of p-levels below 500 hPa
!
! attributes:
!   language: f90
!
!$$$

   implicit none

!  Declare passed variables
   character(len=2)                                ,intent(in   )  :: obsvar
   real(r_double),allocatable,dimension(:)         ,intent(inout)  :: lat,lon,stn_elev,pressure,height
   real(r_double),allocatable,dimension(:)         ,intent(inout)  :: prep_qc_mark,raob,time
   character(8)  ,allocatable,dimension(:)         ,intent(inout)  :: sid
   character(8)  ,allocatable,dimension(:)         ,intent(out  )  :: raob_sid_list
   integer(i_kind),allocatable,dimension(:)        ,intent(out  )  :: sid_plevels ! the number of p-levels
   integer(i_kind),allocatable,dimension(:)        ,intent(out  )  :: sid_plevels_500 ! the number of p-levels below 500 hPa
   real(r_double),allocatable,dimension(:),optional,intent(inout)  :: raob_v

   real(r_double),allocatable,dimension(:) :: lat_so,lon_so,stn_elev_so,pressure_so,height_so
   real(r_double),allocatable,dimension(:) :: prep_qc_mark_so, raob_so, time_so
   real(r_double),allocatable,dimension(:) :: raob_v_so
   character(8),allocatable,dimension(:)   :: sid_so

!  newly sorted
   real(r_double),allocatable,dimension(:)  :: lat1,lon1,stn_elev1,pressure1,height1
   real(r_double),allocatable,dimension(:)  :: prep_qc_mark1,raob1,time1
   real(r_double),allocatable,dimension(:)  :: raob_v1
   character(8),allocatable,dimension(:)    :: sid1

   integer(i_kind)    :: k, indexx
   integer(i_kind)    :: number_of_sid

!------------------------------------------------------------------
!  Raob station ID list (remove dups)
!------------------------------------------------------------------

   if (allocated(raob_sid_list)) deallocate(raob_sid_list)
   ! Get the number of unique sid
   k = 1
   allocate(raob_sid_list(size(sid)))
   raob_sid_list = ''
   raob_sid_list(1) = sid(1)

   do i = 2, size(sid)

!     if sid already exist in raob_sid_list, check next
      if (any( raob_sid_list == sid(i) )) cycle
      ! No match found so add it to the output
      k = k + 1
      raob_sid_list(k) = sid(i)

   end do

   print*, raob_sid_list(1:k)
   number_of_sid = k
   print*, "Sorting: The number of unique sid is", number_of_sid

   deallocate(raob_sid_list)

   ! Assign unique sid & 
   ! get the number of press levels for each sid
   allocate(raob_sid_list(number_of_sid))
   raob_sid_list=""

   ! the number of p levels
   if (allocated(sid_plevels)) deallocate(sid_plevels)
   allocate(sid_plevels(number_of_sid))
   sid_plevels=0 ! Initialization is necessary

   ! the number of p levels below 500 hPa
   if (allocated(sid_plevels_500)) deallocate(sid_plevels_500)
   allocate(sid_plevels_500(number_of_sid))
   sid_plevels_500=0 ! Initialization is necessary

!  when i = 1
   raob_sid_list(1) = sid(1)
   sid_plevels(1)=1
   k = 1

   if (pressure(1)>=500) then
      sid_plevels_500(1)=1
   else
      sid_plevels_500(1)=0
   end if

!  when i > 1
   do i = 2, size(sid)
!     if sid already exist in raob_sid_list
      if (any( raob_sid_list == sid(i) )) then

         m=findloc(raob_sid_list, sid(i), dim=1)
         sid_plevels(m) = sid_plevels(m) + 1

         if (pressure(i) >= 500) then
            sid_plevels_500(m) = sid_plevels_500(m) + 1
         end if

      ! if new sid,
      else

         k = k + 1
         raob_sid_list(k) = sid(i)
         sid_plevels(k) = sid_plevels(k) + 1

         if (pressure(i) >= 500) then
            sid_plevels_500(k) = sid_plevels_500(k) + 1
         end if
 
      end if
   end do

   print*, "Sorting: the number of p-levels for each sid (sid_plevels)=", sid_plevels
   print*, "         the number of p-levels below 500 hPa for each sid (sid_plevels_500)=", sid_plevels_500

   do i = 1, number_of_sid
      print*, "Plevels: i=",i,",raob_sid_list(i)=",raob_sid_list(i)
      print*, "         sid_plevels(all levels)=",sid_plevels(i),",sid_plevels_500(below500hPa)=",sid_plevels_500(i)
   end do
!------------------------------------------------------------------
!  Sort depending on pressure levels for each sid and reassign
!------------------------------------------------------------------

   allocate(lat1(size(sid)), lon1(size(sid)), stn_elev1(size(sid)), pressure1(size(sid)))
   allocate(prep_qc_mark1(size(sid)), raob1(size(sid)), height1(size(sid)))
   allocate(time1(size(sid)), sid1(size(sid)))
   if (obsvar=='uv') allocate(raob_v1(size(sid)))

   do i = 1, number_of_sid ! the number of raob_sid_list
    
      np=sid_plevels(i) ! the number of pressure level depending on sid
      print*, "i=,",i, "np (the number of pressure level for each sid) =", np
      print*, "raob_sid_list(i)=", raob_sid_list(i)
      
      ! Temporary array for sorting according to pressure level
      allocate(lat_so(np), lon_so(np), stn_elev_so(np), pressure_so(np), &
         height_so(np), prep_qc_mark_so(np), raob_so(np), time_so(np), &
         sid_so(np))
      if (obsvar=='uv') allocate(raob_v_so(np))

      k = 1
      do j = 1, size(sid)

         ! Put together the same id information
         if (sid(j) == raob_sid_list(i)) then ! for the same sid

            lat_so(k)=lat(j)
            lon_so(k)=lon(j)
            stn_elev_so(k)=stn_elev(j)
            pressure_so(k)=pressure(j)
            height_so(k)=height(j)
            prep_qc_mark_so(k)=prep_qc_mark(j)
            !setup_qc_mark(j)=setup_qc_mark_t(i)
            raob_so(k)=raob(j)
            if (obsvar=='uv') raob_v_so(k)=raob_v(j)
            time_so(k)=time(j)
            sid_so(k)=sid(j)
            k = k + 1

         end if

      end do
      print*, "pressure_so (before sorting) =", pressure_so
      ! Sorting descending order for pressure level
      do k = 1, np-1 ! except for the last

         location = FindMaximum(pressure_so, k, np) ! find max from this to last
         call Swap(pressure_so(k), pressure_so(location)) ! swap this and the maximum
         call Swap(lat_so(k), lat_so(location)) ! swap this and the maximum
         call Swap(lon_so(k), lon_so(location)) ! swap this and the maximum
         call Swap(stn_elev_so(k), stn_elev_so(location)) ! swap this and the maximum
         call Swap(height_so(k), height_so(location)) ! swap this and the maximum
         call Swap(prep_qc_mark_so(k), prep_qc_mark_so(location)) ! swap this and the maximum
         call Swap(raob_so(k), raob_so(location)) ! swap this and the maximum
         if (obsvar=='uv') call Swap(raob_v_so(k), raob_v_so(location))
         call Swap(time_so(k), time_so(location)) ! swap this and the maximum
         !call Swap(sid_so(k), sid_so(location)) ! swap this and the maximum

      end do
      print*, "Sorting: pressure_so=",pressure_so  
      ! Reassign sorted array

      if (i == 1) then ! for first sid
         indexx=0
      end if
   
      do k = 1, np 

         lat1(indexx+k)=lat_so(k)
         lon1(indexx+k)=lon_so(k)
         stn_elev1(indexx+k)=stn_elev_so(k)
         pressure1(indexx+k)=pressure_so(k)
         height1(indexx+k)=height_so(k)
         prep_qc_mark1(indexx+k)=prep_qc_mark_so(k)
         !setup_qc_mark(j)=setup_qc_mark_t(i)
         raob1(indexx+k)=raob_so(k)
         if (obsvar=='uv') raob_v1(indexx+k)=raob_v_so(k)
         time1(indexx+k)=time_so(k)
         sid1(indexx+k)=sid_so(k)
         if (i>370) print*, "Sorting: k=",k,", sid1(indexx+k)=",sid1(indexx+k),",sid_so(k)=",sid_so(k)

      end do

      indexx=indexx+np
 
      deallocate(lat_so, lon_so, stn_elev_so, pressure_so, &
         height_so, prep_qc_mark_so, raob_so, time_so, &
         sid_so)
      if (obsvar=='uv') deallocate(raob_v_so)

   end do

   ! Reassign the original array with sorted array
   do i = 1, size(sid)
      
      lat(i)=lat1(i)
      lon(i)=lon1(i)
      stn_elev(i)=stn_elev1(i)
      pressure(i)=pressure1(i)
      height(i)=height1(i)
      prep_qc_mark(i)=prep_qc_mark1(i)
      raob(i)=raob1(i)
      if (obsvar=='uv') raob_v(i)=raob_v1(i)
      time(i)=time1(i)
      sid(i)=sid1(i)

   end do

   deallocate(lat1,lon1,stn_elev1,pressure1,height1,prep_qc_mark1,raob1,time1,sid1)
   if (obsvar=='uv') deallocate(raob_v1)

   end subroutine sorting

! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMaximum():
!    This function returns the location of the maximum in the section
! between Start and End.
! --------------------------------------------------------------------

   INTEGER FUNCTION FindMaximum(x, Start, End)
   use kinds, only: r_double,i_kind,r_kind

   IMPLICIT  NONE
   REAL(r_double), ALLOCATABLE,DIMENSION(:), INTENT(IN)  :: x
   INTEGER(i_kind), INTENT(IN)                :: Start, End
   REAL(r_double)                             :: Maximum
   INTEGER(i_kind)                            :: Location
   INTEGER(i_kind)                            :: i

   Maximum  = x(Start)		! assume the first is the max
   Location = Start			! record its position
   DO i = Start+1, End		! start with next elements
      IF (x(i) > Maximum) THEN	!   if x(i) greater than the max?
         Maximum  = x(i)		!      Yes, a new maximum found
         Location = i                !      record its position
      END IF
   END DO
   FindMaximum = Location        	! return the position
   END FUNCTION  FindMaximum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
   use kinds, only: r_double,i_kind,r_kind
   IMPLICIT  NONE
   REAL(r_double), INTENT(INOUT) :: a, b
   REAL(r_double)                :: Temp

   Temp = a
   a    = b
   b    = Temp
   END SUBROUTINE  Swap

!----------------------------------------------------------------------
! Subroutines for converting geopotential height to geometric altitude
!----------------------------------------------------------------------
   subroutine convert_geop_height_to_geom_altitude(latitude, geop_height, geom_altitude)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    convert_geop_height_to_geom_altitude   convert geopotential height to geometric height
!   prgmmr: eyang        org: GMAO/NASA               date: 2023-06-26
!
! program history log:
!   2023-06-26  eyang
!
!   input argument list:
!     latitude    
!     geop_height  - geopotehtial height
!
!   output argument list:
!     geom_altitude - geometric altitude
!
! attributes:
!   language: f90
!
!$$$

   use constants, only: deg2rad, grav_equator, somigliana, eccentricity, grav_ratio, grav
   use constants, only: two
   use constants, only: semi_major_axis, flattening
   use constants, only: init_constants_derived
   use constants, only: init_constants
   implicit none

!  Declare passed variables
   real(r_kind)         ,intent(in)  :: latitude
   real(r_kind)         ,intent(in)  :: geop_height
   real(r_kind)         ,intent(out) :: geom_altitude

!  Declare local variables
   real(r_kind)  :: sin2
   !grav_equator, somigliana, eccentricity
   real(r_kind)  :: termg, termr, termrg
   logical       :: regional
!  Convert geopotential height at layer midpoints to geometric height using
!  equations (17, 20, 23) in MJ Mahoney's note "A discussion of various
!  measures of altitude" (2001).  Available on the web at
!  http://mtp.jpl.nasa.gov/notes/altitude/altitude.html
!
!  termg  = equation 17
!  termr  = equation 21
!  termrg = first term in the denominator of equation 23
!  zges   = equation 23

   regional=.false.
   call init_constants_derived
   call init_constants(regional)

   sin2 = sin(latitude*deg2rad)**2 
   termg = grav_equator * &
           ((one+somigliana*sin2)/sqrt(one-eccentricity*eccentricity*sin2))
   
   termr = semi_major_axis / (one + flattening + grav_ratio - two*flattening*sin2)
   termrg = (termg/grav)*termr

   geom_altitude = (termr*geop_height) / (termrg-geop_height)  ! eq (23) at interface (topo corrected)
   !print*, "geom_altitude=",geom_altitude
   !print*, "geop_height=",geop_height

   end subroutine convert_geop_height_to_geom_altitude

!**********************************************************************
!----------------------------------------------------------------------
! Subroutines for reading diagnostic files
!----------------------------------------------------------------------
   subroutine read_diag(infile,obsvar,obs_type,lat,lon,stn_elev,pressure,&
     height,prep_qc_mark,obs,time,errinv_input,sid,date_time,obs_v)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_diag    read diagnostic file
!   prgmmr: eyang        org: GMAO/NASA               date: 2023-06-26
!
! program history log:
!   2023-06-26  eyang
!
!   input argument list:
!     infile    - diagnostic file name
!     obsvar    - observation variable (t, q, uv)
!
!   output argument list:
!     obs_type  - observation type
!     lat       - observation latitude (degrees)
!     lon       - observation longitude (degrees)
!     stn_elev  - station elevation (meters)
!     pressure  - observation pressure (hPa)
!     height    - observation height (meters)
!     prep_qc_mark - qc flag
!     obs       - observation
!     obs_v     - v observation (optional)
!     time      - observation time (relative to beginning of DA window) [0,6]
!     errinv_input - inverse observation error 
!     sid       - station id
!     date_time - YYYYMMDDHH
!
! attributes:
!   language: f90
!
!$$$

   implicit none

!  Declare passed variables
   character(len=225)                              ,intent(in ) :: infile
   character(len=2)                                ,intent(in ) :: obsvar
   integer(i_kind),allocatable,dimension(:)        ,intent(out) :: obs_type
   real(r_double),allocatable,dimension(:)         ,intent(out) :: lat, lon
   real(r_double),allocatable,dimension(:)         ,intent(out) :: stn_elev, pressure
   real(r_double),allocatable,dimension(:)         ,intent(out) :: height, prep_qc_mark
   real(r_double),allocatable,dimension(:)         ,intent(out) :: obs
   real(r_double),allocatable,dimension(:),optional,intent(out) :: obs_v
   real(r_double),allocatable,dimension(:)         ,intent(out) :: time, errinv_input
   character(8),allocatable,dimension(:)           ,intent(out) :: sid
   integer(i_kind)                                 ,intent(out) :: date_time

!  Declare local variables
   logical :: lexist
   integer(i_kind) :: ncid,ierr,dimid1,dimid2
   integer(i_kind) :: nobs
   integer(i_kind) :: varid1,varid2,varid3,varid4,varid5
   integer(i_kind) :: varid6,varid7,varid8,varid9,varid10
   integer(i_kind) :: varid11
   integer(i_kind) :: varid70,varid80
   integer(i_kind) :: timeLength, strlen
   integer(i_kind) :: i,j
   character(8),allocatable,dimension(:,:) :: sid_org

!  https://github.com/NCAR/rda-prepbufr-decode/blob/main/README.md      

   print*, "Read diagnostic file."
   print*, infile

!  Check if diagnostic file exists
   inquire(file=trim(infile),exist=lexist)
   if (.not.lexist) stop

!--------------------------
!  Read data
!--------------------------
   ierr =  NF90_OPEN(trim(infile),0,ncid)
   if (ierr /= nf90_noerr) call handle_err(ierr,'open')

   ierr = NF90_INQ_DIMID(ncid,'nobs',dimid1)
   if (ierr /= nf90_noerr) call handle_err(ierr,'nobs')
   ierr = nf90_inquire_dimension(ncid, dimid1, len = nobs)
   print*, 'calculate_pblh_from_raob_diag: nobs=', nobs

   ierr = NF90_INQ_DIMID(ncid,'Station_ID_maxstrlen',dimid2)
   if (ierr /= nf90_noerr) call handle_err(ierr,'Station_ID_maxstrlen')
   ierr = nf90_inquire_dimension(ncid, dimid2, len = strlen)
   print*, 'calculate_pblh_from_raob_diag: Station_ID_maxstrlen=', strlen

   allocate(lat(nobs), lon(nobs), stn_elev(nobs), pressure(nobs))
   allocate(height(nobs), prep_qc_mark(nobs), obs(nobs))
   allocate(time(nobs), errinv_input(nobs))
   allocate(sid_org(nobs, strlen), sid(nobs), obs_type(nobs))
   if (obsvar=='uv') allocate(obs_v(nobs))

!  Latitude: degree
   ierr = NF90_INQ_VARID(ncid,'Latitude',varid1)
   if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid1,lat)
!  Longitude: degree
   ierr = NF90_INQ_VARID(ncid,'Longitude',varid2)
   if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid2,lon)
!  Station_Elavation: meter
   ierr = NF90_INQ_VARID(ncid,'Station_Elevation',varid3)
   if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid3,stn_elev)
!!  Pressure: hPa
   ierr = NF90_INQ_VARID(ncid,'Pressure',varid4)
   if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid4,pressure)
!  Height: meter
   ierr = NF90_INQ_VARID(ncid,'Height',varid5)
   if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid5,height)
!  Time: Obs time - Cycle time [hr]
   ierr = NF90_INQ_VARID(ncid,'Time',varid6)
   if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid6,time)
!  Prep_QC_Mark 
   ierr = NF90_INQ_VARID(ncid,'Prep_QC_Mark',varid7)
   if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid7,prep_qc_mark)
!  Setup_QC_Mark 
!   ierr = NF90_INQ_VARID(ncid,'Setup_QC_Mark',varid70)
!   if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid70,setup_qc_mark)

!  Observation 
   if (obsvar == 'uv') then

      ierr = NF90_INQ_VARID(ncid,'u_Observation',varid8)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid8,obs)
      ierr = NF90_INQ_VARID(ncid,'v_Observation',varid80)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid80,obs_v)

   else 

      ierr = NF90_INQ_VARID(ncid,'Observation',varid8)
      if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid8,obs)

   end if

!  Errinv_Input
   ierr = NF90_INQ_VARID(ncid,'Errinv_Input',varid9)
   if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid9,errinv_input)
!  Station ID
   ierr = NF90_INQ_VARID(ncid,'Station_ID',varid10)
   if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid10,sid_org)

   do i = 1, nobs
      sid(i)=sid_org(i,1)
   end do

!  Observation Type
   ierr = NF90_INQ_VARID(ncid,'Observation_Type',varid11)
   if (ierr == nf90_noerr) ierr = NF90_GET_VAR(ncid,varid11,obs_type)

!  Read the global attribute
   ierr = NF90_INQUIRE_ATTRIBUTE(ncid, nf90_global, 'date_time')
   if (ierr == nf90_noerr) ierr = NF90_GET_ATT(ncid,nf90_global,'date_time',date_time)
!  Close the file
   ierr = NF90_CLOSE(ncid)
   if (ierr /= nf90_noerr) call handle_err(ierr,"close")

   end subroutine read_diag 

!**********************************************************************
   subroutine convert_localtime_hr(ltime_yy,ltime_mm, ltime_dd, ltime_hh, interval, ana_yy,ana_mm, ana_dd, ana_hh)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    convert_localtime_hr     convert utc time to local time for 2015 Aug and Sep
!   prgmmr: eyang        org: GMAO/NASA               date: 2023-06-26
!
! program history log:
!   2023-06-26  eyang
!
!   input argument list:
!     interval  - local time minus utc time
!     ana_yy    - year (YYYY) of utc time
!     ana_mm    - month (MM) of utc time
!     ana_dd    - day (DD) of utc time
!     ana_hh    - hour (HH) of utc time
!
!   output argument list:
!     ltime_yy    - year (YYYY) of local time
!     ltime_mm    - month (MM) of local time
!     ltime_dd    - day (DD) of local time
!     ltime_hh    - hour (HH) of local time
!
! attributes:
!   language: f90
!
!$$$

   use kinds, only: i_kind
   integer(i_kind), intent(out) :: ltime_yy, ltime_mm, ltime_dd, ltime_hh
   integer(i_kind), intent(in ) :: interval, ana_yy, ana_mm, ana_dd, ana_hh
   ltime_yy = ana_yy
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
      if (ana_mm == 12 .and. ltime_dd > 31) then
         ltime_yy = ltime_yy + 1
         ltime_mm = 1
         ltime_dd = 1
      end if

   elseif (ltime_hh < 0) then

      ltime_hh = 24 + ltime_hh 
      ltime_dd = ana_dd - 1
      if (ana_mm == 8 .and. ana_dd == 1) then
         ltime_mm = 7
         ltime_dd = 31
      end if
      if (ana_mm == 9 .and. ana_dd == 1) then
         ltime_mm = 8
         ltime_dd = 31
      end if
      if (ana_mm == 12 .and. ana_dd == 1) then
         ltime_mm = 11
         ltime_dd = 30
      end if


   end if
   end subroutine convert_localtime_hr

   subroutine convert_localtime_min(ltime_yy,ltime_mm, ltime_dd, ltime_hh, ltime_min, interval_hr_min, ana_yy,ana_mm, ana_dd, ana_hh,ana_min)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    convert_localtime_min     convert utc time to local time for 2015 Aug/Sep, 2021 Dec
!   prgmmr: eyang        org: GMAO/NASA               date: 2023-10-30
!
! program history log:
!   2023-10-30  eyang
!
!   input argument list:
!     interval  - local time minus utc time
!     ana_yy    - year (YYYY) of utc time
!     ana_mm    - month (MM) of utc time
!     ana_dd    - day (DD) of utc time
!     ana_hh    - hour (HH) of utc time
!     ana_min   - minute (MIN) of utc time
!
!   output argument list:
!     ltime_yy    - year (YYYY) of local time
!     ltime_mm    - month (MM) of local time
!     ltime_dd    - day (DD) of local time
!     ltime_hh    - hour (HH) of local time
!     ltime_min   - minute (MIN) of local time
!
! attributes:
!   language: f90
!
   use kinds, only: i_kind,r_kind
   integer(i_kind), intent(out) :: ltime_yy, ltime_mm, ltime_dd, ltime_hh, ltime_min
   integer(i_kind), intent(in ) :: ana_yy, ana_mm, ana_dd, ana_hh, ana_min
   real(r_kind), intent(in ) :: interval_hr_min
   real(r_kind) :: ana_hh_real, ltime_hh_real

   ! localtime = UTC + interval_hr_min
   ltime_yy = ana_yy
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
      if (ana_mm == 12 .and. ltime_dd > 31) then
         ltime_yy = ltime_yy + 1
         ltime_mm = 1
         ltime_dd = 1
      end if

   elseif (ltime_hh < 0) then

      ltime_hh = 24 + ltime_hh
      ltime_dd = ana_dd - 1
      if (ana_mm == 9 .and. ana_dd == 1) then
         ltime_mm = 8
         ltime_dd = 31
      end if
      if (ana_mm == 12 .and. ana_dd == 1) then
         ltime_mm = 11
         ltime_dd = 30
      end if

   end if
   end subroutine convert_localtime_min

!**********************************************************************
   subroutine convert_obstime_min(obstime_yy,obstime_mm, obstime_dd, obstime_hh,obstime_min, interval, ana_yy,ana_mm, ana_dd, ana_hh)

!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    convert_obstime_min    save obs time as yyyy, mm, dd, hh, and min format
!                                       by adding time, which is relative to analysis time, to analysis time
!                                       (for 2015 Aug and Sep, and 2021 Nov and Dec)
!
! prgmmr: eyang        org: GMAO/NASA               date: 2023-06-26
!
! program history log:
!   2023-06-26  eyang
!
!   input argument list:
!     interval  - obs time relative to analysis time [-3,3]
!     ana_yy  - year (YYYY) of analysis time
!     ana_mm    - month (MM) of analysis time
!     ana_dd    - day (DD) of analysis time
!     ana_hh    - hour (HH) of analysis time
!
!   output argument list:
!     obstime_yy  - year (YYYY) of obs time
!     obstime_mm    - month (MM) of obs time
!     obstime_dd    - day (DD) of obs time
!     obstime_hh    - hour (HH) of obs time
!     obstime_min   - min (mm) of obs time
!
! attributes:
!   language: f90
!
!$$$

   use kinds, only: r_double, i_kind
   real(r_double), intent(out) :: obstime_yy, obstime_mm, obstime_dd, obstime_hh, obstime_min
   real(r_double), intent(in ) :: interval
   integer(i_kind), intent(in ) :: ana_yy, ana_mm, ana_dd, ana_hh

   if (interval.ne.r_missing .and. .not.isnan(interval)) then

      obstime_yy   = ana_yy
      obstime_mm   = ana_mm
      obstime_dd   = ana_dd
      obstime_hh   = ana_hh + interval

      if (obstime_hh >= 24) then

         obstime_hh = obstime_hh - 24
         obstime_dd = ana_dd + 1
         if (ana_mm == 8 .and. obstime_dd > 31) then
            obstime_mm = 9
            obstime_dd = 1
         end if
         if (ana_mm == 9 .and. obstime_dd > 30) then
            obstime_mm = 10
            obstime_dd = 1
         end if
         if (ana_mm == 12 .and. obstime_dd > 31) then
            obstime_yy = obstime_yy+1
            obstime_mm = 1
            obstime_dd = 1
         end if

      elseif (obstime_hh < 0) then

         obstime_hh = 24 + obstime_hh
         obstime_dd = ana_dd - 1
         if (ana_mm == 8 .and. ana_dd == 1) then
            obstime_mm = 7
            obstime_dd = 31
         end if
         if (ana_mm == 9 .and. ana_dd == 1) then
            obstime_mm = 8
            obstime_dd = 31
         end if
         if (ana_mm == 12 .and. ana_dd == 1) then
            obstime_mm = 11
            obstime_dd = 30
         end if

      end if

      ! HH and Minutes
      obstime_min = 60*(obstime_hh - floor(obstime_hh))
      obstime_hh = floor(obstime_hh)

   end if

   end subroutine convert_obstime_min

   end program calculate_pblh_from_raob_diag
