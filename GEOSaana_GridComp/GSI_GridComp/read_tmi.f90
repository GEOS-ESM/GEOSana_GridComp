subroutine read_tmi(mype,val_tmi,ithin,rmesh,jsatid,gstime,&
     infile,lunout,obstype,nread,ndata,nodata,twind,sis,&
     mype_root,mype_sub,npe_sub,mpi_comm_sub)

!$$$  subprogram documentation block
! subprogram:    read_tmi           read TMI 1B11 bufr data
!   prgmmr: j.jin    modified from read_ssmi.f90    date: 2012-12-12
!
! abstract:  This routine reads BUFR format TMI 1B11 radiance 
!            (brightness temperature) files.  Optionally, the data 
!            are thinned to a specified resolution using simple 
!            quality control (QC) checks.
!            QC performed in this subroutine are
!             1.obs time check  |obs-anal|<time_window
!             2.remove overlap orbit
!             3.climate check  reject for tb<tbmin or tb>tbmax 
!
!            When running the gsi in regional mode, the code only
!            retains those observations that fall within the regional
!            domain
!
! program history log:
!   2012-12-12  j.jin   - read_tmi.f90, modified code of read_ssmi.f90. 
!   2013-10-21  j.jin   - add handling and call to zensun to calculate solar zenith/azimuth_ang angles.
!   2014-01-27  j.jin   - cleaup, update after TMI 1B11 bufr files' re-contructure.
!   2014-02-04  j.jin   - added obs thinning procedures:
!                         1), skip obs that far from the grid box center (dist> 0.75);      
!                         2), skip some obs at the ends of some scans.
!                         3), skip obs which show rain or thick cloud or clw/tpwc cannot be retrieved.
!   2014-02-20  j.jin   - added use_edges choices. Obs at scan edges are not read in if use_edges=F.
!                         If use_edges=T, only keep part of the obs at the edges referring the 
!                         mean number of obs between them. Edges are determined in the scaninfo file. 
!   2014-02-22  j.jin   - limit the special thinning processes made on 2014-02-04 in global models.

!   input argument list:
!     mype     - mpi task id
!     val_tmi - weighting factor applied to super obs
!     ithin    - flag to thin data
!     rmesh    - thinning mesh size (km)
!     jsatid   - satellite to read  ex. 'f15'
!     gstime   - analysis time in minutes from reference date
!     infile   - unit from which to read BUFR data
!     lunout   - unit to which to write data for further processing
!     obstype  - observation type to process
!     twind    - input group time window (hours)
!     sis      - satellite/instrument/sensor indicator
!     mype_root - "root" task for sub-communicator
!     mype_sub - mpi task id within sub-communicator
!     npe_sub  - number of data read tasks
!     mpi_comm_sub - sub-communicator for data read
!
!   output argument list:
!     nread    - number of BUFR SSM/I observations read (after eliminating orbit overlap)
!     ndata    - number of BUFR SSM/I profiles retained for further processing (thinned)
!     nodata   - number of BUFR SSM/I observations retained for further processing (thinned)
!
! attributes:
!   language: f90
!
! Note:
!   2013-10-21  j.jin   - there is not a procedure for isfcalc.
!                         (isfcalc - specifies method to determine surface fields
!                         within a FOV. When it is equal to one, integrate
!                         model fields over a FOV. When it is not equal to one, bilinearly
!                         interpolate model fields at a FOV center.)
!$$$  end documentation block
  use kinds, only: r_kind,r_double,i_kind
  use satthin, only: super_val,itxmax,makegrids,map2tgrid,destroygrids, &
      checkob,finalcheck,score_crit
  use radinfo, only: iuse_rad,jpch_rad,nusis,nuchan,nst_gsi,nstinfo,use_edges,radedge1,radedge2
  use gridmod, only: diagnostic_reg,regional,rlats,rlons,nlat,nlon,&
      tll2xy,txy2ll
  use constants, only: deg2rad,rad2deg,zero,one,two,three,four,r60inv,rearth
  use gsi_4dvar, only: l4dvar,iwinbgn,winlen
  use deter_sfc_mod, only: deter_sfc
  use gsi_nstcouplermod, only: gsi_nstcoupler_skindepth, gsi_nstcoupler_deter
  use clw_mod,  only:  retrieval_mi

  implicit none

! Declare passed variables
  character(len=*),intent(in   ) :: infile,obstype,jsatid
  character(len=*),intent(in   ) :: sis
  integer(i_kind), intent(in   ) :: mype,lunout,ithin
  integer(i_kind), intent(in   ) :: mype_root
  integer(i_kind), intent(in   ) :: mype_sub
  integer(i_kind), intent(in   ) :: npe_sub
  integer(i_kind), intent(in   ) :: mpi_comm_sub
  real(r_kind)   , intent(in   ) :: rmesh,gstime,twind
  real(r_kind)   , intent(inout) :: val_tmi
  integer(i_kind),intent(inout)  :: nread
  integer(i_kind),intent(inout)  :: ndata,nodata

! Declare local parameters
  integer(i_kind),parameter :: maxinfo=34
  integer(i_kind)           :: maxchanl

  real(r_kind),  allocatable, dimension(:) :: tbob
  integer(i_kind),allocatable,dimension(:) :: tbmin     ! different tbmin for the 9 channels (from the data document).
  real(r_double), allocatable,dimension(:) :: mirad,  & ! TBB from strtmbr
                                              fovn      ! FOVN 
  real(r_kind),parameter    :: r360=360.0_r_kind
  character(80),parameter   :: satinfo='SAID SIID OGCE GSES SACV'          !use for ufbint()
  character(80),parameter   :: hdr1b='YEAR MNTH DAYS HOUR MINU SECO ORBN'  !use for ufbint()
  integer(i_kind),parameter :: ntime=7, ninfo=5                                     !time header
  real(r_kind)              :: tbmax, satinfo_v(ninfo)
  real(r_double),dimension(ntime):: bfr1bhdr

  character(40)             :: strscan                                     !for TMI, 'TMISQ SCLAT SCLON HMSL NGQI'
  integer(i_kind),parameter :: nloc=5                                      !location dat used for ufbint()
  real(r_double),dimension(nloc) :: midat                                  !location data from 

  character(40),parameter   :: strloc='CLATH CLONH'                        !use for ufbint() 
  character(40),parameter   :: strsaza='SAZA'                              !use for ufbint() 
  real(r_double)            :: pixelloc(2),pixelsaza                       !location data
  character(40),parameter   :: strtmbr='TMBR', strfovn='FOVN'              !use for ufbrep()  
  character(8)              :: subset
  real(r_kind), parameter   :: bmiss=990_r_kind   ! miss values are 999 in tmi bufr
                                                  ! undefined value is 1.0e+11 in tmi_bufr data files.

  integer(i_kind),parameter :: n_angls=3
  real(r_double),dimension(n_angls) :: val_angls

! Declare local variables
  logical        :: tmi,assim,outside,iuse, no85GHz

  integer(i_kind):: i,k,ntest,ireadsb,ireadmg,irec,isub,next,j
  integer(i_kind):: iret,idate,nchanl
  integer(i_kind):: isflg,nreal,idomsfc
  integer(i_kind):: nmind,itx,nele,itt
  integer(i_kind):: iskip
  integer(i_kind):: lnbufr
  integer(i_kind):: ilat,ilon

  real(r_kind) :: sfcr
  real(r_kind) :: pred
  real(r_kind) :: sstime,tdiff,t4dv
  real(r_kind) :: crit1,dist1
  real(r_kind) :: timedif
  real(r_kind),allocatable,dimension(:,:):: data_all

  real(r_kind) :: disterr,disterrmax,dlon00,dlat00

  integer(i_kind) :: nscan,jc,bufsat,js,ij,npos,n, npos_bin
  integer(i_kind),dimension(5):: iobsdate
  real(r_kind):: flgch
  real(r_kind),dimension(0:3):: sfcpct
  real(r_kind),dimension(0:3):: ts
  real(r_kind),dimension(0:4):: rlndsea
  real(r_kind) :: tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10
  real(r_kind) :: zob,tref,dtw,dtc,tz_tr

  real(r_kind):: dlat,dlon,dlon_earth,dlat_earth
  real(r_kind):: sat_def_ang,sat_zen_ang        ! default and derived satellite zenith angle
  real(r_kind):: sat_scan_ang,sat_azimuth_ang, azimuth_ang

! ---- sun glint ----
  integer(i_kind):: doy,mlen(12),mday(12),mon,m
  real(r_kind)   :: time_4_sun_glint_calc,clath_sun_glint_calc,clonh_sun_glint_calc
  real(r_kind)   :: sun_zenith,sun_azimuth_ang
  data  mlen/31,28,31,30,31,30, &
             31,31,30,31,30,31/


  integer(i_kind) :: ang_nn, pos_max 
  integer(i_kind),allocatable :: pos_statis(:), npos_all(:,:)

! ---- check clw ----
  real(r_kind)    :: clw,tpwc
  integer(i_kind) :: kraintype,ierrret,nchanl2
! ---- skip some obs at the beginning and end of a scan ----
  integer(i_kind):: radedge_min,radedge_max,iscan_pos,iedge_log,j2

!**************************************************************************
! Initialize variables
  lnbufr = 15
  disterrmax=zero
  ntest=0
  iscan_pos = 8     ! id in data_all for scan positions
  iedge_log  = 32    ! id in data_all for log if obs is to be obleted beause of locating near scan edges.

  ndata  = 0
  nodata = 0
  nread  = 0
  sat_def_ang =52.8_r_kind   ! default TMI satellite zenith angle.

  ilon=3 
  ilat=4

  if(nst_gsi>0) then   
     call gsi_nstcoupler_skindepth(obstype, zob)         ! get penetration depth (zob) for the obstype
  endif
  m = 0
  do mon=1,12
     mday(mon) = m
     m = m + mlen(mon)
  end do


! Set various variables depending on type of data to be read

  tmi  = obstype  == 'tmi'

     maxchanl = 9                          ! number of channels
     nchanl = 9                            ! 9 channls
     nscan  = 208                          ! number of pixels (for high resolution)
     npos_bin = 3                          ! max number of high resolution pixels at a position 
                                           !     (for grouping the obs at a position)
     if(jsatid == 'trmm')bufsat=282        ! Satellite ID (is 282 in bufr file)
     allocate (tbmin(maxchanl))
     tbmin = (/33,66,133,80,133,133,112,70,70/) ! read from TMI (v7) HDF data document.
     tbmax = 320.0_r_kind                       ! one value for all tmi channels (see data document).
     strscan='TMISQ SCLAT SCLON HMSL NGQI' 
     ang_nn=nscan/npos_bin+1
     allocate (tbob(maxchanl), mirad(maxchanl), fovn(maxchanl))
     rlndsea(0) = zero
     rlndsea(1) = 30._r_kind
     rlndsea(2) = 30._r_kind
     rlndsea(3) = 30._r_kind
     rlndsea(4) = 100._r_kind

! If all channels of a given sensor are set to monitor or not
! assimilate mode (iuse_rad<1), reset relative weight to zero.
! We do not want such observations affecting the relative
! weighting between observations within a given thinning group.

  radedge_min = 0
  radedge_max = 1000
  assim=.false.
  search: do i=1,jpch_rad
    if (trim(nusis(i))==trim(sis)) then
        if (radedge1(i)/=-1 .and. radedge2(i)/=-1) then
           radedge_min=radedge1(i)
           radedge_max=radedge2(i)
        end if
        if (iuse_rad(i)>=0) then
           if (iuse_rad(i)>0) assim=.true.
           if (assim) exit
        endif
     endif
  end do search
  if (.not.assim) val_tmi=zero

  nchanl2=nchanl - 2    ! cha 1&2 (10 GHz) for TMI are not available for SSMI (clw/tpwc retrievals).
  no85GHz=.false.

! Make thinning grids
  call makegrids(rmesh,ithin)

! Open unit to satellite bufr file
  open(lnbufr,file=infile,form='unformatted')
  call openbf(lnbufr,'IN',lnbufr)
  call datelen(10)

! Allocate arrays to hold data
  nreal  = maxinfo + nstinfo
  nele   = nreal   + nchanl
  allocate(data_all(nele,itxmax))

!       Extract satellite id from the 1st MG.  If it is not the one we want, exit reading.
        call readmg(lnbufr, subset, iret, idate)
        rd_loop: do while (ireadsb(lnbufr)==0)
          call ufbint(lnbufr,satinfo_v,ninfo,1,iret,satinfo)
          if(nint(satinfo_v(1)) /= bufsat) then 
            write(6,*) 'READ_TMI: Bufr satellie ID SAID', nint(satinfo_v(1)), &
                       ' does not match ', bufsat
            go to 690
          endif
        enddo rd_loop
! Big loop to read data file
  next=0
  read_subset: do while(ireadmg(lnbufr,subset,idate)>=0) ! TMI scans
     next=next+1
     if(next == npe_sub)next=0
     if(next /= mype_sub)cycle
     read_loop: do while (ireadsb(lnbufr)==0)            ! TMI pixels
        call ufbrep(lnbufr,fovn,1,nchanl,iret, strfovn)
        ! npos must .LE. 90 because nstep=90 in bias_angle correction code
        !   ../../../../Applications/NCEP_Etc/NCEP_bias/main.f90
        npos = fovn(nchanl)/npos_bin + 1  ! always group scan positions according channel-9's positions. 
        if (.not. use_edges .and. &
             (npos < radedge_min .OR. npos > radedge_max ))  cycle read_loop

! ----- extract time information  
        call ufbint(lnbufr,bfr1bhdr,ntime,1,iret,hdr1b)

!       calc obs seqential time. If time is outside window, skip this obs
        iobsdate(1:5) = bfr1bhdr(1:5) !year,month,day,hour,min
        call w3fs21(iobsdate,nmind)
        t4dv=(real(nmind-iwinbgn,r_kind) + real(bfr1bhdr(6),r_kind)*r60inv)*r60inv
        if (l4dvar) then
           if (t4dv<zero .OR. t4dv>winlen) cycle read_loop
        else
           sstime=real(nmind,r_kind) + real(bfr1bhdr(6),r_kind)*r60inv
           tdiff=(sstime-gstime)*r60inv
           if(abs(tdiff) > twind)  cycle read_loop
        endif

! ----- Read header record to extract obs location information  
        call ufbint(lnbufr,midat,nloc,1,iret,strscan)
        call ufbint(lnbufr,pixelloc,2,1,iret,strloc)
        call ufbint(lnbufr,pixelsaza,1,1,iret,strsaza)

!---    Extract brightness temperature data.  Apply gross check to data. 
!       If obs fails gross check, reset to missing obs value.
        call ufbrep(lnbufr,mirad,1,nchanl,iret,strtmbr)

!          Regional case
           dlat_earth = pixelloc(1)  !deg
           dlon_earth = pixelloc(2)  !deg
           if(abs(dlat_earth)>90.0_r_kind .or. abs(dlon_earth)>r360) cycle read_loop
           if(dlon_earth< zero) dlon_earth = dlon_earth+r360
           if(dlon_earth==r360) dlon_earth = dlon_earth-r360
           dlat_earth = dlat_earth*deg2rad
           dlon_earth = dlon_earth*deg2rad

           if(regional)then
              call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
              if(diagnostic_reg) then
                 call txy2ll(dlon,dlat,dlon00,dlat00)
                 ntest=ntest+1
                 disterr=acos(sin(dlat_earth)*sin(dlat00)+cos(dlat_earth)*cos(dlat00)* &
                      (sin(dlon_earth)*sin(dlon00)+cos(dlon_earth)*cos(dlon00)))*rad2deg
                 disterrmax=max(disterrmax,disterr)
              end if

!             Check to see if in domain
              if(outside) cycle read_loop

!          Global case
           else
              dlat = dlat_earth  
              dlon = dlon_earth  
              call grdcrd1(dlat,rlats,nlat,1)
              call grdcrd1(dlon,rlons,nlon,1)
           endif
!          If available, set value of tmi zenith angle
           if (midat(1) < bmiss ) then  ! SAZA is readable if TMI dataQuality available.
              sat_zen_ang = pixelsaza
           else
              sat_zen_ang = sat_def_ang
           endif
           sat_scan_ang = asin( sin(sat_zen_ang*deg2rad)*rearth/(rearth+midat(4)) )
           sat_azimuth_ang = azimuth_ang(pixelloc(2),pixelloc(1),midat(3),midat(2))

           !  -------- Retreive Sun glint angle -----------
           clath_sun_glint_calc = pixelloc(1)
           clonh_sun_glint_calc = pixelloc(2)
           if(clonh_sun_glint_calc > 180._r_kind) clonh_sun_glint_calc = clonh_sun_glint_calc - 360.0_r_kind
           doy = mday( int(bfr1bhdr(2)) ) + int(bfr1bhdr(3))
           if ((mod( int(bfr1bhdr(1)),4)==0).and.( int(bfr1bhdr(2)) > 2))  then
              doy = doy + 1
           end if
           time_4_sun_glint_calc = bfr1bhdr(4)+bfr1bhdr(5)*r60inv+bfr1bhdr(6)*r60inv*r60inv
           call zensun(doy,time_4_sun_glint_calc,clath_sun_glint_calc,clonh_sun_glint_calc,sun_zenith,sun_azimuth_ang)
           ! output solar zenith angles are between -90 and 90
           sun_zenith = 90.-sun_zenith                                ! make sure solar zenith angles are between 0 and 180
        
!          Transfer observed brightness temperature to work array.  
!          If any temperature exceeds limits, or data_quality /= 0
!          reset observation to "bad" value.
           iskip=0
           do jc=1, nchanl 
              if(mirad(jc)<tbmin(jc) .or. mirad(jc)>tbmax &
                 .or. midat(1) > 0 .or. midat(5) > 0 ) then           ! skip data with data_quality > 0 or geoQuality > 0
                 iskip = iskip + 1
              else
                 nread=nread+1
              end if
              tbob(jc) = mirad(jc) 
           enddo
 
           ! if data for more than 7 channels are bad, skip. This removes the
           ! high frequency obs at locations where there are not low freqency obs.
           ! In other words, only use data at locations where both low and high
           ! frequcy obs are made. 
           if(tmi .and. iskip >= nchanl-2)  cycle read_loop 

           ! flgch = iskip*two   !used for thinning priority range 0-18
           ! JJJ   12/21/2012
           ! There are data for all channels 1-9 at half of the number of pixels while there are 
           ! only channels 8-9 at the other half of locations. Therefore, flgch is not a
           ! universal flag for thinning priority here. As a result, this flag
           ! is not used and flgch is set = 0.
           flgch = 0

           if (l4dvar) then
              crit1 = 0.01_r_kind+ flgch
           else
              timedif = 6.0_r_kind*abs(tdiff) ! range: 0 to 18
              crit1 = 0.01_r_kind+timedif + flgch
           endif

!          Map obs to thinning grid
           call map2tgrid(dlat_earth,dlon_earth,dist1,crit1,itx,ithin,itt,iuse,sis)  !JJ???
           if(.not. iuse)cycle read_loop

           ! if the obs is far from the grid box center, do not use it.
           if(.not. regional .and. dist1 > 0.75) cycle read_loop  
           crit1 = crit1 + 10._r_kind * float(iskip)
           call checkob(dist1,crit1,itx,iuse)
           if(.not. iuse)cycle read_loop

!          Locate the observation on the analysis grid.  Get sst and land/sea/ice
!          mask.  

!       isflg    - surface flag
!                  0 sea
!                  1 land
!                  2 sea ice
!                  3 snow
!                  4 mixed                     

           call deter_sfc(dlat,dlon,dlat_earth,dlon_earth,t4dv,isflg,idomsfc,sfcpct, &
              ts,tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10,sfcr)
           if(.not. regional .and. isflg==0 .and. &
              tmi) then
             ! retrieve tpwc and clw
             call retrieval_mi(tbob(3),nchanl2,no85GHz, &
               tpwc,clw,kraintype,ierrret )
             ! don't use obsercation for which there are thick cloudsi or rain, or clw or twc cannot be retrieved.
             if(tpwc<0 .or. clw > 3 .or. kraintype /=0_i_kind .or. ierrret > 0) cycle read_loop
           endif

           crit1 = crit1 + rlndsea(isflg)
           call checkob(dist1,crit1,itx,iuse)
           if(.not. iuse)cycle read_loop

           call finalcheck(dist1,crit1,itx,iuse)
           if(.not. iuse)cycle read_loop

!          interpolate NSST variables to Obs. location and get dtw, dtc, tz_tr
!
           if(nst_gsi>0) then
              tref  = ts(0)
              dtw   = zero
              dtc   = zero
              tz_tr = one
              if(sfcpct(0)>zero) then
                 call gsi_nstcoupler_deter(dlat_earth,dlon_earth,t4dv,zob,tref,dtw,dtc,tz_tr)
              endif
           endif

!          Transfer observation parameters to output array.  
           data_all( 1,itx) = bufsat              ! satellite id
           data_all( 2,itx) = t4dv                ! time diff between obs and anal (min)
           data_all( 3,itx) = dlon                ! grid relative longitude
           data_all( 4,itx) = dlat                ! grid relative latitude
           data_all( 5,itx) = sat_zen_ang*deg2rad ! local (satellite) zenith angle (radians)
           data_all( 6,itx) = sat_azimuth_ang     ! local (satellite) azimuth_ang angle (degrees)
           data_all( 7,itx) = sat_scan_ang        ! scan(look) angle (rad)
           data_all( 8,itx) = npos                ! scan position,  .le. 90
           data_all( 9,itx) = sun_zenith          ! solar zenith angle (deg)
           data_all(10,itx) = sun_azimuth_ang         ! solar azimuth_ang angle (deg)
           data_all(11,itx) = sfcpct(0)           ! sea percentage of
           data_all(12,itx) = sfcpct(1)           ! land percentage
           data_all(13,itx) = sfcpct(2)           ! sea ice percentage
           data_all(14,itx) = sfcpct(3)           ! snow percentage
           data_all(15,itx)= ts(0)                ! ocean skin temperature
           data_all(16,itx)= ts(1)                ! land skin temperature
           data_all(17,itx)= ts(2)                ! ice skin temperature
           data_all(18,itx)= ts(3)                ! snow skin temperature
           data_all(19,itx)= tsavg                ! average skin temperature
           data_all(20,itx)= vty                  ! vegetation type
           data_all(21,itx)= vfr                  ! vegetation fraction
           data_all(22,itx)= sty                  ! soil type
           data_all(23,itx)= stp                  ! soil temperature
           data_all(24,itx)= sm                   ! soil moisture
           data_all(25,itx)= sn                   ! snow depth
           data_all(26,itx)= zz                   ! surface height
           data_all(27,itx)= idomsfc + 0.001_r_kind ! dominate surface type
           data_all(28,itx)= sfcr                 ! surface roughness
           data_all(29,itx)= ff10                 ! ten meter wind factor
           data_all(30,itx)= dlon_earth*rad2deg   ! earth relative longitude (degrees)
           data_all(31,itx)= dlat_earth*rad2deg   ! earth relative latitude (degrees)
           data_all(iedge_log,itx) = 0            ! =0, not to be obsoleted as at scan edges
           data_all(maxinfo-1,itx)= val_tmi
           data_all(maxinfo,itx)= itt

           if(nst_gsi>0) then
              data_all(maxinfo+1,itx) = tref       ! foundation temperature
              data_all(maxinfo+2,itx) = dtw        ! dt_warm at zob
              data_all(maxinfo+3,itx) = dtc        ! dt_cool at zob
              data_all(maxinfo+4,itx) = tz_tr      ! d(Tz)/d(Tr)
           endif

           do i=1,nchanl
              data_all(i+nreal,itx)=tbob(i)
           end do

     end do read_loop
  end do read_subset
690 continue
  call closbf(lnbufr)

! If multiple tasks read input bufr file, allow each tasks to write out
! information it retained and then let single task merge files together

  call combine_radobs(mype_sub,mype_root,npe_sub,mpi_comm_sub,&
     nele,itxmax,nread,ndata,data_all,score_crit)
     write(6,*) 'READ_TMI: after combine_obs, nread,ndata is ',nread,ndata

!=========================================================================================================
  if( use_edges .and. (radedge_min > 1 .or. radedge_max < ang_nn).and. mype_sub==mype_root )then
    ! Obsolete some obs at the beginning and end positions of a scan by flagging
    !       obs at these positions with negative NPOS values.
    ! Note: This is an arbitary process. Just want to phase out part of these obs
    !       at the scan edges in the QC process (qc_ssmi, ifail_scanedge_qc=58).
    !       However, there is not a known quality issue at the edge of scans.
    !       JJJ, 2/12/2014
     pos_max=ndata  
     allocate(pos_statis(ang_nn))
     allocate(npos_all(pos_max,ang_nn))
     npos_all = 0
     pos_statis = 0
     do n=1,ndata
        i = nint(data_all(iscan_pos,n))
        pos_statis(i) = pos_statis(i) + 1
        npos_all(pos_statis(i), i) = n
     enddo

     do n=1, ndata
        i = nint(data_all(iscan_pos,n))
        if(i < radedge_min .or. i > radedge_max) then
          data_all(iedge_log,n) = 1     ! assume all at scan edges at the beginning.
        endif
     enddo
     if( radedge_min > 1 )then
       pos_max = sum(pos_statis(radedge_min : (radedge_min+1)))/2 
       do i=radedge_min-1, 1, -1
         if(pos_max==0) then
           j2=1
         else
           j2=nint(float(pos_statis(i))/pos_max)
           j2=max(1,j2)
         endif
         do j=1,pos_statis(i),j2
           n = npos_all(j,i)
           data_all(iedge_log,n)= 0     ! flag back
         enddo      
       enddo
     endif

     if( radedge_max < ang_nn )then
       pos_max = sum(pos_statis((radedge_max-1) : radedge_max))/2 
       do i=radedge_max+1,ang_nn
         if(pos_max==0) then
           j2=1
         else
           j2=nint(float(pos_statis(i))/pos_max)
           j2=max(1,j2)
         endif
         do j=1,pos_statis(i),j2
           n = npos_all(j,i)
           data_all(iedge_log,n)= 0     ! flag back
         enddo      
       enddo
     endif

     ! new pos_statis
     pos_statis=0
     do n=1,ndata
        i = nint(data_all(iscan_pos,n))
        if(data_all(iedge_log,n)>0) cycle
        pos_statis(i) = pos_statis(i) + 1
     enddo
     write(6,*) 'READ_', trim(obstype), ': after obsolete_obs near edges, ndata ', sum(pos_statis)
     !write(6,*) 'READ_', trim(obstype), ': number of observations: '
     !write(6, '(5x, 10I10)')  pos_statis
     deallocate(pos_statis, npos_all)
  endif ! use_edges, but flag part of obs at the scan edges with negative FOV values.
!=========================================================================================================

! Allow single task to check for bad obs, update superobs sum,
! and write out data to scratch file for further processing.
  if (mype_sub==mype_root.and.ndata>0) then

!    Identify "bad" observation (unreasonable brightness temperatures).
!    Update superobs sum according to observation location

     do n=1,ndata
        do i=1,nchanl
           if(data_all(i+nreal,n) > tbmin(i) .and. &
              data_all(i+nreal,n) < tbmax)nodata=nodata+1
        end do
        itt=nint(data_all(maxinfo,n))
        super_val(itt)=super_val(itt)+val_tmi
     end do
!    Write final set of "best" observations to output file
     write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
     write(lunout) ((data_all(k,n),k=1,nele),n=1,ndata)
  
  endif

! Deallocate data arrays
  deallocate(data_all)


! Deallocate satthin arrays
1000 continue
  call destroygrids

  if(diagnostic_reg .and. ntest>0 .and. mype_sub==mype_root) &
     write(6,*)'READ_TMI:  mype,ntest,disterrmax=',&
        mype,ntest,disterrmax

  deallocate(tbmin, tbob, mirad, fovn)
! End of routine
 return
end subroutine read_tmi

 
real(r_kind) function azimuth_ang(lona,lata,lonb,latb)

!$$$  subprogram documentation block
!   prgmmr: j.jin                                 date: 2013-10-22
!
! abstract:  derive azimuth_ang angle from location_1 to location_2. 
!
! program history log:
!   2013-10-22  j.jin   - debut
!   2013-12-02  j.jin   - Re-write.
!
!   input argument list:
!     lona,lata     - location array(longitude,latitude,degree) for location_1.
!     lonb,latb     - location array(longitude,latitude,degree) for location_2.
!
!   output argument list:
!     azimuth_ang  - angle between the direction of location_2 looking from 
!                    location_1 and the North direction.
!             
! attributes:
!   language: f90
!$$$
   use kinds, only: r_kind,i_kind
   use constants, only: deg2rad

   implicit none
   real(r_kind),intent(in)  :: lona,lata,lonb,latb
   real(r_kind)             :: lon1,lat1,lon2,lat2, lon_d, dfi, &
                               angles_sn(0:1), angles_we(0:1)
   integer(i_kind)          :: idir
   data angles_sn /180.0_r_kind,0.0_r_kind/
   data angles_we /270.0_r_kind,90.0_r_kind/

   lon1=lona*deg2rad
   lat1=lata*deg2rad
   lon2=lonb*deg2rad
   lat2=latb*deg2rad
   lon_d  = asin( sin(lon2)*cos(lon1) - cos(lon2)*sin(lon1) )
   dfi = cos(lat1)*tan(lat2) - sin(lat1)*cos(lon_d)
   if(lon_d == 0.0_r_kind) then
     if(lat2 /= lat1) then
         idir = int((lat2-lat1)/abs(lat2-lat1) + 1)/2
         azimuth_ang = angles_sn(idir)
     else
         azimuth_ang = 0.0_r_kind
     endif
   else
     if(dfi /= 0.0_r_kind) then
       azimuth_ang = atan(sin(lon_d)/dfi)/deg2rad
       if(lon_d > 0 .and. azimuth_ang < 0) then
         azimuth_ang = azimuth_ang + 180.0_r_kind
       else if(lon_d < 0 .and. azimuth_ang > 0) then
         azimuth_ang = azimuth_ang + 180.0_r_kind
       else if(lon_d < 0 .and. azimuth_ang < 0) then
         azimuth_ang = azimuth_ang + 360.0_r_kind
       else
       endif
     else
        idir = int(lon_d/abs(lon_d) + 1)/2
        azimuth_ang = angles_we(idir)
     endif
   endif
   return
end function azimuth_ang
