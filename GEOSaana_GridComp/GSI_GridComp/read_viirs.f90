subroutine read_sst_viirs(mype,val_viirs,ithin,rmesh,jsatid,&
     gstime,infile,lunout,obstype,nread,ndata,nodata,twind,sis, &
     mype_root,mype_sub,npe_sub,mpi_comm_sub,nobs, &
     nrec_start,dval_use)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    read_sst_viirs,   read M-Band (12, 15 and 16) VIIRS radiance data from the SST VIIRS data file
!
!   prgmmr: li, xu           org: np23                date: 2019-08-29
!
! abstract:  This routine reads BUFR format VIIRS 1b radiance (brightness
!            temperature) from the SST VIIRS files
!
!            When running the gsi in regional mode, the code only
!            retains those observations that fall within the regional
!            domain
!
! program history log:
!   2026-05-22  a.lee   - added read_viirs in GEOS system.
!
!   input argument list:
!     mype     - mpi task id
!     val_viirs- weighting factor applied to super obs
!     ithin    - flag to thin data
!     rmesh    - thinning mesh size (km)
!     jsatid   - satellite to read
!     gstime   - analysis time in minutes from reference date
!     infile   - unit from which to read BUFR data
!     lunout   - unit to which to write data for further processing
!     obstype  - observation type to process
!     twind    - input group time window (hours)
!     sis      - satellite/instrument/sensor indicator
!
!     mype_root    - "root" task for sub-communicator
!     mype_sub     - mpi task id within sub-communicator
!     npe_sub      - number of data read tasks
!     mpi_comm_sub - sub-communicator for data read
!     nrec_start   - first subset with useful information
!     dval_use     - true. if any dval weighting is used for satellite data
!
!   output argument list:
!     nread    - number of BUFR VIIRS observations read
!     ndata    - number of BUFR VIIRS profiles retained for further processing
!     nodata   - number of BUFR VIIRS observations retained for further processing
!     nobs     - array of observations on each subdomain for each processor
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,r_single,r_double,i_kind
  use satthin, only: super_val,itxmax,makegrids,map2tgrid,destroygrids, &
                     checkob, finalcheck,score_crit
  use satthin, only: radthin_time_info,tdiff2crit
  use obsmod,  only: time_window_max, ianldate
  use gridmod, only: diagnostic_reg,regional,nlat,nlon,tll2xy,txy2ll,rlats,rlons
  use constants, only: deg2rad, zero, one, two, half, rad2deg, r60inv
  use radinfo, only: iuse_rad,jpch_rad,nusis
  use gsi_4dvar, only: l4dvar,l4densvar,iwinbgn,winlen
  use deter_sfc_mod, only: deter_sfc
  use obsmod, only: bmiss
  use gsi_nstcouplermod, only: nst_gsi,nstinfo
  use gsi_nstcouplermod, only: gsi_nstcoupler_skindepth, gsi_nstcoupler_deter_viirs
  use sfcobsqc, only: get_sunangle
  use mpimod, only: npe
  implicit none


! Declare passed variables
  character(len=*),                intent(in  )  :: infile,obstype,jsatid
  character(len=20),               intent(in  )  :: sis
  integer(i_kind),                 intent(in  )  :: mype,lunout,ithin,nrec_start
  integer(i_kind),                 intent(inout) :: nread
  integer(i_kind), dimension(npe), intent(inout) :: nobs
  integer(i_kind),                 intent(inout) :: ndata,nodata
  real(r_kind),                    intent(in   ) :: rmesh,gstime,twind
  real(r_kind),                    intent(inout) :: val_viirs
  integer(i_kind),                 intent(in   ) :: mype_root
  integer(i_kind),                 intent(in   ) :: mype_sub
  integer(i_kind),                 intent(in   ) :: npe_sub
  integer(i_kind),                 intent(in   ) :: mpi_comm_sub
  logical,                         intent(in   ) :: dval_use

! Declare local parameters
!  character(6),parameter:: file_sst='SST_AN'
!  integer(i_kind),parameter:: mlat_sst = 3000
!  integer(i_kind),parameter:: mlon_sst = 5000
  real(r_kind),parameter:: r6=6.0_r_kind
  real(r_kind),parameter:: scan_start=-56.28_r_kind, scan_inc=1.26472_r_kind
  real(r_double),parameter:: r360=360.0_r_double
  real(r_kind),parameter:: tbmin=50.0_r_kind
  real(r_kind),parameter:: tbmax=550.0_r_kind

  real(r_kind),parameter :: r_earth=6370.0_r_kind, sat_hgt=829.0_r_kind
  real(r_kind),parameter :: nviirs=6304.0_r_kind,nfov=90.0_r_kind
  character(len=80),parameter ::  &
     headr='YEAR MNTH DAYS HOUR MINU SECO CLATH CLONH SAID SAZA'

! Declare local variables  
  logical outside,iuse,assim
  character(len=8) :: subset

  integer(i_kind) klon1,klatp1,klonp1,klat1
  integer(i_kind) nchanl,iret,ich_m15
  integer(i_kind) idate,maxinfo
  integer(i_kind) ilat,ilon
  integer(i_kind),dimension(5):: idate5
  integer(i_kind) nmind,isflg,idomsfc
  integer(i_kind) itx,k,i,bufsat,n
  integer(i_kind) nreal,nele,itt
!  integer(i_kind) nlat_sst,nlon_sst
  integer(i_kind) ksatid

  real(r_kind) dlon,dlat,rsc
  real(r_kind) dlon_earth,dlat_earth,sfcr
  real(r_kind) dlon_earth_deg,dlat_earth_deg
  real(r_kind) w00,w01,w10,w11,dx1,dy1
  real(r_kind) pred,crit1,tdiff,sstime,dx,dy,dist1
!  real(r_kind) dlat_sst,dlon_sst,sst_hires
  real(r_kind) t4dv, t4dv_sun

  real(r_kind),dimension(0:4):: rlndsea
  real(r_kind),dimension(0:3):: sfcpct
  real(r_kind),dimension(0:3):: ts
!  real(r_kind),dimension(mlat_sst):: rlats_sst
!  real(r_kind),dimension(mlon_sst):: rlons_sst
!  real(r_kind),allocatable,dimension(:,:):: sst_an
  real(r_kind),allocatable,dimension(:,:):: data_all
  real(r_kind) :: tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10
  real(r_kind) :: zob,tref,dtw,dtc,tz_tr
  real(r_kind) :: scan_ang,sat_zen,dfov,r01
  real(r_single) :: sol_zen
  integer(i_kind) :: scan_pos

  real(r_double), dimension(10) :: hdr
  real(r_double), dimension(2,3) :: bufrf
  integer(i_kind) lnbufr,ireadsb,ireadmg,iskip,irec,next
  integer(i_kind),allocatable,dimension(:)::nrec

  real(r_kind) cdist,disterr,disterrmax,dlon00,dlat00
  integer(i_kind) ntest

  real(r_kind)    :: ptime,timeinflat,crit0
  integer(i_kind) :: ithin_time,n_tbin,it_mesh

  character(len=10) :: str
  integer(i_kind) :: idate_sun

!**************************************************************************

! Start routine here.  Set constants.  Initialize variables
  maxinfo = 31
  lnbufr = 10
  disterrmax=zero
  ntest=0
  ndata  = 0
  nodata  = 0
  nread   = 0
  nchanl = 5
  ich_m15    = 2 
  r01 = 0.01_r_kind

  dfov = (nviirs - one)/nfov

  ilon=3
  ilat=4

  if(nst_gsi>0) then
     call gsi_nstcoupler_skindepth(obstype, zob)         ! get penetration depth (zob) for the obstype
  endif

  rlndsea(0) = zero
  rlndsea(1) = 30._r_kind
  rlndsea(2) = 20._r_kind
  rlndsea(3) = 30._r_kind
  rlndsea(4) = 30._r_kind


  if (jsatid == 'npp') then
     bufsat = 224
  elseif (jsatid == 'n20' .or. jsatid == 'j1') then
     bufsat = 225
  elseif (jsatid == 'n21' .or. jsatid == 'j2') then
     bufsat = 226
  else
     write(*,*) 'READ_SST_VIIRS: Unrecognized value for jsatid '//jsatid//':RETURNING'
     return
  end if


! If all channels of a given sensor are set to monitor or not
! assimilate mode (iuse_rad<1), reset relative weight to zero.
! We do not want such observations affecting the relative
! weighting between observations within a given thinning group.

  assim=.false.
  search: do i=1,jpch_rad
     if ((nusis(i)==sis) .and. (iuse_rad(i)>0)) then
        assim=.true.
        exit search
     endif
  end do search
  if (.not.assim) val_viirs=zero


  call radthin_time_info(obstype, jsatid, sis, ptime, ithin_time)
  if( ptime > 0.0_r_kind) then
     n_tbin=nint(2*time_window_max/ptime)
  else
     n_tbin=1
  endif
! Make thinning grids
  call makegrids(rmesh,ithin)

! Read hi-res sst analysis
!     allocate(sst_an(mlat_sst, mlon_sst))
!     call rdgrbsst(file_sst,mlat_sst,mlon_sst,&
!     sst_an,rlats_sst,rlons_sst,nlat_sst,nlon_sst)

! Allocate arrays to hold all data for given satellite
  if(dval_use) maxinfo = maxinfo + 2
  nreal = maxinfo + nstinfo
  nele  = nreal   + nchanl
  allocate(data_all(nele,itxmax),nrec(itxmax))

  open(lnbufr,file=trim(infile),form='unformatted')         ! open bufr data file

! Associate the tables file with the message file, and identify the 
! latter to BUFRLIB software
  call openbf (lnbufr,'IN',lnbufr)

  next=0
  nrec=999999
  irec=0
! Read BUFR VIIRS 1b data
  read_msg: do while (ireadmg(lnbufr,subset,idate) >= 0)
     irec=irec+1
     if(irec < nrec_start) cycle read_msg
     next=next+1
     if(next == npe_sub)next=0
     if(next /= mype_sub)cycle
     read_loop: do while (ireadsb(lnbufr) == 0)
            call ufbint(lnbufr,hdr,10,1,iret,headr)
            call ufbrep(lnbufr,bufrf,2,3,iret,'CHNM TMBR')
        if(iret <= 0) cycle read_loop
        ksatid  = nint(hdr(9))                ! Extract satellite id from bufr file
        if(ksatid /= bufsat) cycle read_loop  ! If this sat is not the one we want, read next record
 
        iskip = 0
            do k=1,nchanl-2
               if(bufrf(2,k) < zero .or. bufrf(2,k) > tbmax) then
               iskip=iskip+1
           end if
        end do
        if(iskip >= nchanl)cycle read_loop

!       Extract date information.  If time outside window, skip this obs
        idate5(1) = nint(hdr(1))    !year
        idate5(2) = nint(hdr(2))    !month
        idate5(3) = nint(hdr(3))    !day
        idate5(4) = nint(hdr(4))    !hour
        idate5(5) = nint(hdr(5))    !minute
        rsc       = hdr(6)          !second in real
        call w3fs21(idate5,nmind)
        write(str,'(I4.4,I2.2,I2.2,I2.2)') idate5(1), idate5(2), idate5(3), idate5(4)
        read(str,*) idate_sun
        t4dv=(real((nmind-iwinbgn),r_kind) + rsc*r60inv)*r60inv
        t4dv_sun=(real(idate5(5),r_kind) + rsc*r60inv)*r60inv
        sstime=real(nmind,r_kind) + rsc*r60inv
        tdiff=(sstime-gstime)*r60inv

        if (l4dvar.or.l4densvar) then
           if (t4dv<zero .OR. t4dv>winlen) cycle read_loop
        else
           if (abs(tdiff) > twind) cycle read_loop
        endif

!       Convert obs location to radians
        if (abs(hdr(7))>90.0_r_double .or. abs(hdr(8))>r360) cycle read_loop
        if (hdr(8)==r360) hdr(8)=hdr(8)-r360
        if (hdr(8)< zero) hdr(8)=hdr(8)+r360

        dlon_earth_deg = hdr(8)
        dlat_earth_deg = hdr(7)

        dlon_earth = hdr(8)*deg2rad   !convert degrees to radians
        dlat_earth = hdr(7)*deg2rad
        sat_zen =  hdr(10)

!       Regional case
        if(regional)then
           call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)
           if(diagnostic_reg) then
              call txy2ll(dlon,dlat,dlon00,dlat00)
              ntest=ntest+1
              cdist=sin(dlat_earth)*sin(dlat00)+cos(dlat_earth)*cos(dlat00)* &
                   (sin(dlon_earth)*sin(dlon00)+cos(dlon_earth)*cos(dlon00))
              cdist=max(-one,min(cdist,one))
              disterr=acos(cdist)*rad2deg
              disterrmax=max(disterrmax,disterr)
           end if
           if(outside) cycle read_loop

!       Global case: Get relative lat/lon in analysis grids
        else
           dlat = dlat_earth
           dlon = dlon_earth
           call grdcrd1(dlat,rlats,nlat,1)
           call grdcrd1(dlon,rlons,nlon,1)
        endif

        nread = nread + 1

        call deter_sfc(dlat,dlon,dlat_earth,dlon_earth,t4dv,isflg,idomsfc,sfcpct, &
           ts,tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10,sfcr)
        if(abs(sfcpct(0) - 1.0) > zero)  cycle read_loop

        if(nst_gsi>0) then
           tref  = ts(0)
           dtw   = zero
           dtc   = zero
           tz_tr = one
           if(sfcpct(0)>zero) then
              !call gsi_nstcoupler_deter(dlat_earth,dlon_earth,t4dv,zob,tref,dtw,dtc,tz_tr)
              call gsi_nstcoupler_deter_viirs(dlat_earth,dlon_earth,t4dv,zob,tref,dtw,dtc,tz_tr)
           if(tref < 1e-9)  cycle read_loop
           endif
        endif

        crit0 = 0.01_r_kind
        timeinflat=6.0_r_kind
        call tdiff2crit(tdiff,ptime,ithin_time,timeinflat,crit0,crit1,it_mesh)
        call map2tgrid(dlat_earth,dlon_earth,dist1,crit1,itx,ithin,itt,iuse,sis,it_mesh=it_mesh)
        if(.not. iuse)cycle read_loop

!       Interpolate hi-res sst analysis to observation location
!          dlat_sst = dlat_earth
!          dlon_sst = dlon_earth
!          call grdcrd1(dlat_sst,rlats_sst,nlat_sst,1)
!          call grdcrd1(dlon_sst,rlons_sst,nlon_sst,1)
!
!          klon1=int(dlon_sst); klat1=int(dlat_sst)
!          dx  =dlon_sst-klon1; dy  =dlat_sst-klat1
!          dx1 =one-dx;         dy1 =one-dy
!          w00=dx1*dy1; w10=dx1*dy; w01=dx*dy1; w11=dx*dy
!
!          klat1=min(max(1,klat1),nlat_sst); klon1=min(max(0,klon1),nlon_sst)
!          if(klon1==0) klon1=nlon_sst
!          klatp1=min(nlat_sst,klat1+1); klonp1=klon1+1
!          if(klonp1==nlon_sst+1) klonp1=1
!
!           sst_hires=w00*sst_an(klat1,klon1 ) + w10*sst_an(klatp1,klon1 ) + &
!                     w01*sst_an(klat1,klonp1) + w11*sst_an(klatp1,klonp1)
!
             !####################################################### Only process the obs. at the location where hsst(sst_an,
             !in this case?)  is in the current processed range             
!       if ( sst_hires < zero ) then
!          print*,' sst_hires,klat1,klon1 : ',sst_hires,klat1,klon1
!       endif

!     "Score" observation.   We use this information to id "best" obs.

!     Locate the observation on the analysis grid.  Get sst and land/sea/ice
!     mask.  

!     isflg    - surface flag
!                0 sea
!                1 land
!                2 sea ice
!                3 snow
!                4 mixed                          
!        call deter_sfc(dlat,dlon,dlat_earth,dlon_earth,t4dv,isflg,idomsfc,sfcpct, &
!           ts,tsavg,vty,vfr,sty,stp,sm,sn,zz,ff10,sfcr)
!        if(isflg /= zero)  print*, 'really? in viirs?' 
!        if(isflg /= zero)  cycle read_loop
        !if(isflg /= zero .or. abs(sfcpct(0)) - 1.0) .lt. 0.0001)  cycle read_loop

        call checkob(dist1,crit1,itx,iuse)
        if(.not. iuse)cycle read_loop


!       Set common predictor parameters

!       Compute "score" for observation.  All scores>=0.0.  Lowest score is "best"

        pred   = (600.0_r_kind - bufrf(2,ich_m15)) * r01

        crit1=crit1+rlndsea(isflg)
        crit1 = crit1+pred  
        call finalcheck(dist1,crit1,itx,iuse)

        if(.not. iuse)cycle read_loop
!
!          calculate scan angle
!
           scan_ang = asin(r_earth/(r_earth+sat_hgt)*sin(sat_zen*deg2rad))/deg2rad
!
!          get scan position (1 to 90)
!
           scan_pos = 1+(scan_ang-scan_start)/scan_inc
           if ( scan_pos > nfov ) scan_pos = nfov
!
!          calculate solar zenith angle

!
        call get_sunangle(idate_sun,t4dv_sun,dlon_earth,dlat_earth,sol_zen)
        !call get_sunangle(ianldate,t4dv,dlon_earth,dlat_earth,sol_zen)

!       interpolate NSST variables to Obs. location and get dtw, dtc, tz_tr
!
!        if(nst_gsi>0) then
!           tref  = ts(0)
!           dtw   = zero
!           dtc   = zero
!           tz_tr = one
!           if(sfcpct(0)>zero) then
!              call gsi_nstcoupler_deter(dlat_earth,dlon_earth,t4dv,zob,tref,dtw,dtc,tz_tr)
!              !call gsi_nstcoupler_deter_viirs(dlat_earth,dlon_earth,t4dv,zob,tref,dtw,dtc,tz_tr)
!           !if(tref == zero)  cycle read_loop 
!           endif
!        endif
           
!       Transfer information to work array
        data_all(1, itx) = hdr(9)                 ! satellite id 
        data_all(2, itx) = t4dv                   ! time (hours)
        data_all(3, itx) = dlon                   ! grid relative longitude
        data_all(4, itx) = dlat                   ! grid relative latitude
        data_all(5, itx) = sat_zen*deg2rad        ! satellite zenith angle (radians)
        data_all(6, itx) = bmiss                  ! satellite azimuth angle
        data_all(7, itx) = scan_ang*deg2rad       ! scan angle
        data_all(8, itx) = real(scan_pos)         ! scan position
        data_all(9, itx) = 90.0_r_kind-sol_zen    ! solar zenith angle (degree)
        data_all(10,itx) = bmiss                  ! solar azimuth angle (radians)
        data_all(11,itx) = sfcpct(0)              ! sea percentage of
        data_all(12,itx) = sfcpct(1)              ! land percentage
        data_all(13,itx) = sfcpct(2)              ! sea ice percentage
        data_all(14,itx) = sfcpct(3)              ! snow percentage
        data_all(15,itx) = ts(0)                  ! ocean skin temperature (from surface file: sst_full)
        data_all(16,itx) = ts(1)                  ! land skin temperature (from surface file: sst_full)
        data_all(17,itx) = ts(2)                  ! ice skin temperature (from surface file: sst_full)
        data_all(18,itx) = ts(3)                  ! snow skin temperature (from surface file: sst_full)
        data_all(19,itx) = tsavg                  ! average skin temperature
        data_all(20,itx) = vty                    ! vegetation type
        data_all(21,itx) = vfr                    ! vegetation fraction
        data_all(22,itx) = sty                    ! soil type
        data_all(23,itx) = stp                    ! soil temperature
        data_all(24,itx) = sm                     ! soil moisture
        data_all(25,itx) = sn                     ! snow depth
        data_all(26,itx) = zz                     ! surface height
        data_all(27,itx) = idomsfc + 0.001_r_kind ! dominate surface type
        data_all(28,itx) = sfcr                   ! surface roughness
        data_all(29,itx) = ff10                   ! ten meter wind factor
        data_all(30,itx) = dlon_earth_deg         ! earth relative longitude (degrees)
        data_all(31,itx) = dlat_earth_deg         ! earth relative latitude (degrees)
        if(dval_use)then
           data_all(32,itx) = val_viirs              ! weighting factor applied to super obs
           data_all(33,itx) = itt                    !
        end if

        if(nst_gsi>0) then
           data_all(maxinfo+1,itx) = tref            ! foundation temperature
           data_all(maxinfo+2,itx) = dtw             ! dt_warm at zob
           data_all(maxinfo+3,itx) = dtc             ! dt_cool at zob
           data_all(maxinfo+4,itx) = tz_tr           ! d(Tz)/d(Tr)
        endif

           data_all(1+nreal,itx) = bufrf(2,1)
           data_all(2+nreal,itx) = bmiss
           data_all(3+nreal,itx) = bmiss
           data_all(4+nreal,itx) = bufrf(2,2)
           data_all(5+nreal,itx) = bufrf(2,3)

           nrec(itx)=irec

!    End of satellite read block

     enddo read_loop
  enddo read_msg
  call closbf(lnbufr)

  call combine_radobs(mype_sub,mype_root,npe_sub,mpi_comm_sub,&
     nele,itxmax,nread,ndata,data_all,score_crit,nrec)

! Allow single task to check for bad obs, update superobs sum,
! and write out data to scratch file for further processing.
  if (mype_sub==mype_root.and.ndata>0) then
   do n=1,ndata
      do k=1,nchanl
         if(data_all(k+nreal,n) > tbmin .and. &
            data_all(k+nreal,n) < tbmax) nodata=nodata+1
      end do
   end do
   if(dval_use .and. assim)then
      do n=1,ndata
         itt=nint(data_all(33,n))
         super_val(itt)=super_val(itt)+val_viirs
      end do
   end if
 
!  Write retained data to local file
   call count_obs(ndata,nele,ilat,ilon,data_all,nobs)
   write(lunout) obstype,sis,nreal,nchanl,ilat,ilon
   write(lunout) ((data_all(k,n),k=1,nele),n=1,ndata)
  endif
! write(6,*) 'READ_VIIRS:  mype, total number of obs info, nread,ndata : ',mype, nread,ndata

! Deallocate local arrays
  deallocate(data_all,nrec)
!  deallocate(sst_an)

! Deallocate arrays
  call destroygrids

  if(diagnostic_reg.and.ntest>0 .and. mype_sub==mype_root) &
     write(6,*)'READ_VIIRS-m:  ',&
     'mype,ntest,disterrmax=',mype,ntest,disterrmax

! End of routine
  return
end subroutine read_sst_viirs
