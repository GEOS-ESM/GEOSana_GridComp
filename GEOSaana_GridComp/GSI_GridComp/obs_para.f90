subroutine obs_para(ndata,mype)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    obs_para    assign and distribute observations to subdomains
!   prgmmr: weiyu yang       org: np20                date: 1998-05-27
!
! abstract: Based on observation location and domain decomposition, assign
!           and distribute observations to subdomains
!
! program history log:
!   1998-05-27  weiyu yang
!   1999-08-24  derber, j., treadon, r., yang, w., first frozen mpp version
!   2004-06-15  treadon, reformat documenation
!   2004-07-23  derber - modify to include conventional sst
!   2004-07-28  treadon - add only to module use, add intent in/out
!   2004-11-19  derber - modify to eliminate file and change to use logical 
!                        rather than weights
!   2005-06-14  wu      - add OMI total ozone
!   2005-09-08  derber - simplify data set handling
!   2005-10-17  treadon - fix bug in defnition of mype_diaghdr
!   2006-02-03  derber  - modify for new obs control and obs count
!   2006-07-28  derber  - eliminate counting of data and clean up
!   2008-04-29  safford - rm unused vars and uses
!   2008-06-30  derber  - optimize calculation of criterion
!   2008-09-08  lueken  - merged ed's changes into q1fy09 code
!   2009-04-21  derber  - reformulate to remove communication
!   2008-05-10  meunier - handle for lagrangian data
!   2013-01-26  parrish - attempt fix for bug flagged by WCOSS debug compiler.
!                            Replace 
!                              "call dislag(.....,nobs_s)"
!                            with
!                              "call dislag(.....,nobs_s(mm1))"
!                           nobs_s is an array in current subroutine, but is a
!                           scalar inside subroutine dislag.
!   2014-10-03  carley  - add creation mpi subcommunicator needed for 
!                          buddy check QC to distinguish among pe subdomains 
!                          with and without obs (only for t obs and twodvar_regional 
!                          at the moment) 
!   2015-10-28  guo     - added ioid_s(:) to support split-observer GSI with a
!                         control vector grid different from the guess state
!                         grid.
!   2023-11-22  eyang   - add mpi_allreduce for inflating ens spread (coefficient array)
!   2023-11-29  eyang   - add option to decide to inflate ens spread near pblh or not
!   2023-12-20  eyang   - add inflating ens spread near pblrf only over ocean (idomsfc=0)
!   2024-03-01  zhu     - handles pblrf as pblrf + pblkh (secondary minimum)
!
!   input argument list:
!     ndata(*,1)- number of prefiles retained for further processing
!     ndata(*,2)- number of observations read
!     ndata(*,3)- number of observations keep after read
!     mype     - mpi task number
!     ipoint   - pointer in array containing information about all obs type to process
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: i_kind
  use constants, only: zero
  use jfunc, only: factqmin,factqmax
  use mpimod, only: npe,mpi_itype,mpi_comm_world,ierror
  use obsmod, only: obs_setup,dtype,mype_diaghdr,ndat,nsat1
  use obsmod, only: obsfile_all,dplat,nobs_sub,obs_sub_comm 
  use gridmod, only: twodvar_regional 
  use qcmod, only: buddycheck_t,buddydiag_save 
  use gsi_io, only: verbose, print_obs_para

! mpi_allreduce
  use hybrid_ensemble_parameters, only: ges_coef_inf_ens_grid_global_t
  use hybrid_ensemble_parameters, only: ges_coef_inf_ens_grid_global_q
  use hybrid_ensemble_parameters, only: ges_coef_inf_ens_grid_global_t1
  use hybrid_ensemble_parameters, only: ges_coef_inf_ens_grid_global_q1
  use hybrid_ensemble_parameters, only: ges_vlocal_ens_grid_global
  use hybrid_ensemble_parameters, only: ges_vlocal_ens_grid_global1
  use hybrid_ensemble_parameters, only: InfEnsprdPBLH ! if ture, inflate ensemble spread near pblh
  use hybrid_ensemble_parameters, only: grd_ens

  use mpimod, only: ierror,mpi_comm_world,mpi_rtype4,mpi_sum,mpi_max

  implicit none

! Declare passed variables
  integer(i_kind)                  ,intent(in   ) :: mype
  integer(i_kind),dimension(ndat,3),intent(in   ) :: ndata

! Declare local variables
  integer(i_kind) lunout,is,ii
  integer(i_kind) mm1
  integer(i_kind) ndatax_all,ikey_yes,ikey_no,newprocs,newrank 
  integer(i_kind),dimension(npe):: ikey,icolor 
  integer(i_kind) nobs_s
  logical print_verbose
 
  integer(i_kind) :: i,j,k
  print_verbose=.false.
  if(verbose)print_verbose=.true.
!
!****************************************************************
! Begin obs_para here
!
! Distribute observations as a function of pe number.
  nsat1=0
  obs_sub_comm=0 
  mype_diaghdr = -999
  mm1=mype+1
  ndatax_all=0
  lunout=55

!
! Set number of obs on each pe for each data type
  open(lunout,file=obs_setup,form='unformatted')
  rewind lunout
  !print*, 'obs_para L98: ndat=',ndat ! ndat=103
  do is=1,ndat

     if(dtype(is) /= ' ' .and. ndata(is,1) > 0)then
!        if(dtype(is)=='pblri') then
!           print*, 'opara L102: dtype(is)=',dtype(is)
!           print*, 'opara L103: mm1,is=',mm1,is ! is=14
!           print*, 'opara L104: ndata(is,1)=',ndata(is,1) ! 624
!        end if

        ndatax_all=ndatax_all + ndata(is,1)

        if (dtype(is)=='lag') then    ! lagrangian data
           call dislag(ndata(is,1),mm1,lunout,obsfile_all(is),dtype(is),&
                nobs_s) 
        nsat1(is)= nobs_sub(mm1,is)
        else 
           obproc:do ii=1,npe
             if(nobs_sub(ii,is) > 0)then
               mype_diaghdr(is) = ii-1
               exit obproc
             end if
           end do obproc
                
           if(nobs_sub(mm1,is) > zero)then    ! classical observations
!              if(dtype(is)=='pblri')then
!                 print*, 'opara L121: mm1,is=',mm1,is
!                 print*, 'opara L122: nobs_sub(mm1,is)=',nobs_sub(mm1,is)
!                 print*, 'opara L123: obsfile_all(is)=',obsfile_all(is)
!                 print*, 'opara L124: ndata(is,1)=',ndata(is,1)
!                 print*, 'opara L125: lunout=',lunout
!              end if
              call disobs(ndata(is,1),nobs_sub(mm1,is),mm1,lunout, &
                  obsfile_all(is),dtype(is))
           end if

        end if
        nsat1(is)= nobs_sub(mm1,is)
        if(mm1 == 1 .and. (print_verbose .or. print_obs_para))then
           write(6,1000)dtype(is),dplat(is),(nobs_sub(ii,is),ii=1,npe)
1000       format('OBS_PARA: ',2A10,8I10,/,(10X,10I10))
        end if

        ! Simple logic to organize which tasks do and do not have obs. 
        !  Needed for buddy check QC.   
        if (twodvar_regional .and. dtype(is) == 't' .and.  buddycheck_t) then 
           ikey_yes=0 
           ikey_no=0 
           ikey=0 
           do ii=1,npe 
              if (nobs_sub(ii,is)>0) then 
                 icolor(ii)=1 
                 ikey(ii)=ikey_yes 
                 ikey_yes=ikey_yes+1 
              else 
                 icolor(ii)=2 
                 ikey(ii)=ikey_no 
                 ikey_no=ikey_no+1 
              end if 
           end do 
 
           ! With organized colors and keys, now create the new MPI communicator 
           !   which only talks to pe's who have obs on their subdomains.  This is 
           !   needed for MPI communication within the setup* routines (e.g. a buddy check). 
               
           call mpi_comm_split(mpi_comm_world,icolor(mm1),ikey(mm1),obs_sub_comm(is),ierror)   
           CALL MPI_COMM_SIZE(obs_sub_comm(is), newprocs, ierror) 
           CALL MPI_COMM_RANK(obs_sub_comm(is), newrank, ierror) 
           if (buddydiag_save) write(6,'(A,I3,I10,A,I20,A,I3,A,I3)') 'obs_para: mype/myobs=',& 
                              mype,nobs_sub(mm1,is),'newcomm=',obs_sub_comm(is),'newprocs=', & 
                              newprocs,'newrank=',newrank            
        end if 
        
     end if

  end do
  close(lunout)


  !-------------------
  ! When inflating ensemble spread near pblh
  ! - mpi_allreduce -
  !   Collect and assign max values ens_grid_global array for coef_inf_ens_grid and vlocal_ens_grid
  if (InfEnsprdPBLH) then
     call mpi_allreduce(ges_coef_inf_ens_grid_global_t,ges_coef_inf_ens_grid_global_t1,size(ges_coef_inf_ens_grid_global_t1),&
          mpi_rtype4,mpi_max,mpi_comm_world,ierror)
     call mpi_allreduce(ges_coef_inf_ens_grid_global_q,ges_coef_inf_ens_grid_global_q1,size(ges_coef_inf_ens_grid_global_q1),&
          mpi_rtype4,mpi_max,mpi_comm_world,ierror)
     call mpi_allreduce(ges_vlocal_ens_grid_global,ges_vlocal_ens_grid_global1,size(ges_vlocal_ens_grid_global1),&
          mpi_rtype4,mpi_max,mpi_comm_world,ierror)
  end if
  !-------------------

! If there are no obs available, turn off moisture constraint.  
! If the user still wants the moisture constraint active when no obs are
! present, comment out the block of code below.
  if (ndatax_all == 0) then
     factqmin=zero
     factqmax=zero
     if (mype==0) write(6,*)'OBS_PARA:  ***WARNING*** no observations to be  ',&
          ' assimilated. reset factqmin,factqmax=',factqmin,factqmax
  endif


  return
end subroutine obs_para

subroutine disobs(ndata,nobs,mm1,lunout,obsfile,obstypeall)

!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    disobs  distribute observations into each pe subdomain
!   prgmmr: weiyu yang                                date: 1998-04-06
!
! abstract: distribute observations into each pe subdomain
!
! program history log:
!   1998-04-06  weiyu yang
!   1999-08-24  derber, j., treadon, r., yang, w., first frozen mpp version
!   2004-06-15  treadon, reformat documenation
!   2004-07-28  treadon - add only to module use, add intent in/out
!   2004-11-19  derber - change to eliminate additional file and use logical
!                        rather than weights
!   2005-10-17  treadon - remove conversion to grid-relative cooridnates
!   2006-04-06  middlecoff - changed lunin from 15 to 11
!   2008-04-29  safford - rm unused vars and uses
!   2008-09-08  lueken  - merged ed's changes into q1fy09 code
!   2009-04-21  derber  - reformulate to remove communication
!   2023-09-26  eyang   - add code to define array regarding PBLH location to modify hybrid vertical localization length scale near PBLH obs
!   2023-11-13  eyang   - assign pblhobs array on ens grid directly instead of regridding
!   2023-11-28  eyang   - consider different inflation coefficient for t and rh (t: 5, rh: 10) to have similar ratio of max/min ens spread after inflation below 600 hPa (model levels from 1 to 21). For example, before inflation: ratio of T and q (median)= 2.4, 5.3 vs after inflation: ratio of T and q = 8.0, 16.2, which is 3 times ratio for T and q after inflation (inflating 5 and 10 for t and q ens pert., respectively).
!   2023-11-28  eyang   - add option inflating ens spread near pblri or q_gridient_pblh
!   2023-11-29  eyang   - add option to decide to inflate ens spread near pblh or not
!   2023-12-20  eyang   - add inflating ens spread near pblrf only over ocean (idomsfc=0)
!   2024-01-07  eyang   - consider overlapping grid points regarding inflating ens spread
!
!   input argument list:
!     nn_obs   - number of an observation data
!     ndata    - number of observations
!     lunin    - unit from which to read all obs of a given type
!     lunout   - unit to which to write subdomain specific observations
!     mm1      - mpi task number + 1
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
  use gridmod, only: periodic_s,nlon,nlat,jlon1,ilat1,istart,jstart

!--------------------------------------------------
! For inflating ens sprd and modifying hybrid vertical localization length scale near PBLH obs,
! 1) Save pblh obs grid information (global domain in ens grid)
! 2) Modify hybrid vertical localization length scale near PBLH obs
! In this code, we will do 1)
  use guess_grids, only: hrdifsig,nfldsig
  use gridmod, only: regional,rlats,rlons,use_sp_eqspace
  use gridmod, only: grd_a,lat1,lon1,lat2,lon2
  use guess_grids, only: ges_geopi
  use general_sub2grid_mod, only: general_suba2sube
  use hybrid_ensemble_parameters, only: grd_a1, grd_e1
! same as grd_anl/grd_ens, but with communication set up for a single 3d grid
  use gridmod, only: nsig
  use hybrid_ensemble_parameters, only: ges_coef_inf_ens_grid_global_t
  use hybrid_ensemble_parameters, only: ges_coef_inf_ens_grid_global_q
  use hybrid_ensemble_parameters, only: ges_vlocal_ens_grid_global
  use hybrid_ensemble_parameters, only: grd_ens, grd_anl,p_e2a
  use hybrid_ensemble_parameters, only: inf_ensprd_pblri_type ! 1: pblri, 2: q-gradient pblh
  use hybrid_ensemble_parameters, only: InfEnsprdPBLH ! if ture, inflate ensemble spread near pblh
  use hybrid_ensemble_parameters, only: inf_coef_t,inf_coef_q ! coef 10 and 22.5 for Gaussian function inflate en_sprd 5 and 10 tiems
  use hybrid_ensemble_parameters, only: inf_ensprd_pblri ! if true, inflate ens sprd near pblri
  use hybrid_ensemble_parameters, only: inf_ensprd_pblrf ! if true, inflate ens sprd near pblrf (i.e., pblrf_low)
  use hybrid_ensemble_parameters, only: inf_ensprd_pblkh ! if true, inflate ens sprd near pblkh (i.e., pblrf_high)
  use gsi_bundlemod, only: gsi_grid
  use gsi_bundlemod, only: gsi_bundle

  use constants, only: pi ! for Gaussian dist.
  use constants, only: deg2rad,rad2deg,zero ! Find index of obs on ens grid
  use gsi_enscouplermod, only: gsi_enscoupler_localization_grid

!  use mpimod, only: ierror,mpi_comm_world,mpi_rtype,mpi_sum
!  use mpl_allreducemod, only: mpl_allreduce ! to sum ens_grid_global for all processors
!  use hybrid_ensemble_parameters, only: grid_coef_inf_ens_grid_global1
!--------------------------------------------------

  implicit none

! Declare passed variables
  integer(i_kind)               ,intent(in   ) :: ndata,lunout,mm1,nobs
  character(len=*)              ,intent(in   ) :: obsfile
  character(len=*)              ,intent(in   ) :: obstypeall

! Declare local variables
  integer(i_kind) lon,lat,lat_data,lon_data,n,k,lunin
  integer(i_kind) jj,nreal,nchanl,nn_obs
  integer(i_kind) ndata_s
  integer(i_kind),dimension(mm1):: ibe,ibw,ibn,ibs
  logical,allocatable,dimension(:):: luse_s
  integer(i_kind),allocatable,dimension(:):: ioid_s
        ! IOID represents Initial-Obs-ID, which is the serial index of given
        ! observation record of a given observation source, before they are
        ! distributed according to the guess-grid partition.  This ID is saved
        ! for all future needs of obs. redistribution, such that the sorted
        ! sequences can be reproduced.
  real(r_kind),allocatable,dimension(:,:):: obs_data,data1_s
  character(10):: obstype
  character(20):: isis

!--------------------------------------------------
! Vertical localization
  type(gsi_grid)   :: grid_ens,grid_anl
  type(gsi_bundle) :: work_ens,work_anl
  real(r_kind),allocatable :: tmp_ens(:,:,:,:),tmp_anl(:,:,:,:)
! Declare external calls for code analysis
  external:: tintrp2a1
  external:: tintrp2a11
  external:: tintrp2a11_indx
  external:: tintrp2a11_indx_ens_grd
  integer(i_kind) iz,izp,j,t
  integer(i_kind) ix,ixp,iy,iyp,it,itp
  integer(i_kind) ix_ens_global,ixp_ens_global,iy_ens_global,iyp_ens_global
  real(r_kind) dlat,dlon,dtime
  real(r_kind) dlat_ens,dlon_ens
  real(r_kind) lat_deg,lon_deg
  real(r_kind) dlat_earth_deg,dlon_earth_deg
  real(r_kind) dlat_earth,dlon_earth ! radians
  integer(i_kind) mype
  real(r_kind) pblh_obs
  real(r_kind),allocatable,dimension(:) :: geopiges
  real(r_kind) rlats_ens_local(grd_ens%nlat)
  real(r_kind) rlons_ens_local(grd_ens%nlon)

! Inflation coefficient based on Gaussian Distribution
  real(r_kind) sd2 ! variance = square of standard deviation
  real(r_kind) coeff_t ! 10.0 for t (multiplied by pdf)
  real(r_kind) coeff_q ! 22.5 for q (rh) (multiplied by pdf)
  real(r_kind) lev_diff ! 0.0, 1.0, 2.0 : model level - model level closest to pblh
  real(r_kind),dimension(2,2,2) :: inf_coef0_t, inf_coef1_t, inf_coef2_t ! inflation coeff. ( = 1 + pdf * coeff)
  real(r_kind),dimension(2,2,2) :: inf_coef0_q, inf_coef1_q, inf_coef2_q ! inflation coeff. ( = 1 + pdf * coeff)

! Declare local parameters
  real(r_kind),parameter:: r360 = 360.0_r_kind
!--------------------------------------------------

! Read and write header

  do k=1,mm1

!    ibw,ibe,ibs,ibn west,east,south and north boundaries of total region
!    (mm1=0)  0, 36, 0, 9 -> 36,72, 0,9 -> 72,108, 0,9
!    (mm1=17) 0, 36, 9,18 -> 36,72, 9,18-> 72,108, 9,18
     ibw(k)=jstart(k)-1
     ibe(k)=jstart(k)+jlon1(k)-1
     ibs(k)=istart(k)-1
     ibn(k)=istart(k)+ilat1(k)-1

  end do

  lunin=11
  !obsfile = obs_input.xx (obstype)
  open(lunin,file=trim(obsfile),form='unformatted')
  read(lunin)obstype,isis,nreal,nchanl,lat_data,lon_data
  if(trim(obstype) /=trim(obstypeall)) &
        write(6,*)'DISOBS:  ***ERROR***   obstype,obstypeall=',trim(obstype),trim(obstypeall)

  print*, 'disobs L310: obstype=',trim(obstype) ! pblri?
  nn_obs = nreal + nchanl

  allocate(obs_data(nn_obs,ndata))
!  Read in all observations of a given type along with subdomain flags
  read(lunin) obs_data ! cdata_all
  close(lunin)

  ndata_s=0
  allocate(data1_s(nn_obs,nobs),luse_s(nobs),ioid_s(nobs))

! Loop over all observations.  Locate each observation with respect
! to subdomains.
  obs_loop: do n=1,ndata
     lat=obs_data(lat_data,n)
     lat=min(max(1,lat),nlat)

     if(lat>=ibs(mm1).and.lat<=ibn(mm1))then
        lon=obs_data(lon_data,n)
        lon=min(max(0,lon),nlon)
        if((lon >= ibw(mm1).and. lon <=ibe(mm1))  .or.  &
           (lon == 0   .and. ibe(mm1) >=nlon) .or.  &
           (lon == nlon    .and. ibw(mm1) <=1) .or. periodic_s(mm1)) then
           ndata_s=ndata_s+1
           ioid_s(ndata_s)=n
           luse_s(ndata_s)= .true.
           do jj= 1,nn_obs
              data1_s(jj,ndata_s) = obs_data(jj,n)
           end do
           !----------------------------------------
           ! Inflating ens spread and Vertical localzaiton near PBLH obs
           !----------------------------------------
           ! - pblri_raobpblh and pblrf_gnssro
           !   (wrong)   pblri and ikx=5 (nkx=120), cdata_all(8,n)=ikx <-- ikx changes
           !   (correct) pblri and kx=120, cdata_all(17,n)=kx
           if (InfEnsprdPBLH) then
              !obs_data(17,n): index of ob type (kx or nkx)
              !if ( (pblri and 120 raob_pblh) or (pblrf and over water) )
              !   and usage (obs_data(11,n))==0

              !if ( ( (trim(obstype)=='pblri' .and. obs_data(8,n)==5 .and. inf_ensprd_pblri) .or. &
              if ( ( (trim(obstype)=='pblri' .and. obs_data(17,n)==120 .and. inf_ensprd_pblri) .or. &
                     (trim(obstype)=='pblrf' .and. obs_data(15,n)==0 .and. inf_ensprd_pblrf) .or. &
                     (trim(obstype)=='pblkh' .and. obs_data(15,n)==0 .and. obs_data(6,n)==870 .and. inf_ensprd_pblkh) ) &
                      .and. obs_data(11,n)==0) then
                 !----------------------------------------------------
                 ! 1) Obtain grid index (subdomain) near pblri observation location and time
                 !    - on analysis grid (for it, itp, iz)
                 !    - on ensemble grid (for ix_ens_global, iy_ens_global)
                 !--------------------------
                 !subroutine tintrp2a11(f,g,dx,dy,obstime,gridtime, &
                 !     mype,nflds)
                 !subroutine tintrp2a11_indx(dx,dy,obstime,gridtime, &
                 !     mype,nflds,ix,ixp,iy,iyp,itime,itimep)
                 !abstract: same as tintrp2a11 but for horizontal grid indexes (anal) surrounding
                 !           an observation point
                 !     ix,iy,ixp,iyp - horizontal grid indexes on anal grid
                 !     itime,itimep  - time grid indexes
                 !---------------------------------------------------
                 !    1a) ANALYSIS GRID -> to get it, itp (index of model time w.r.t. obs time) for inflating ens sprd coef and
                 !                         iz (index of model vertical level w.r.t. pblh obs) for inflating and vertical localization
                 !---------------------------------------------------
                 ! dlat, dlon <-- full analysis grid
                 ! ix, ixp <-- subdomain analysis grid
                 ! it, itp <-- we will use this info
                 ! itime=7, obs_data(itime,n): index of obs time in data array
                 mype=mm1-1
                 dtime=obs_data(7,n) ! index of obs time in data array
                 dlat=obs_data(lat_data,n)
                 dlon=obs_data(lon_data,n)
!                 print*, 'disobs 120 raob_pblh L362: obs_data(8,n),mm1=',obs_data(8,n),mm1
!                 print*, 'disobs L363: index: dlat,dlon,dtime,mm1=',dlat,dlon,dtime,mm1
!                 print*, 'disobs L364: lat_earth_deg,lon_earth_deg,mm1=',obs_data(13,n),obs_data(12,n),mm1
                 !print*, 'disobs L365: jstart(mm1),mm1=',jstart(mm1),mm1
                 !print*, 'disobs L366: istart(mm1),mm1=',istart(mm1),mm1
!                 print*, 'disobs L367_mype1: rlats(1),rlats(2),rlats(90)=',rlats(1),rlats(2),rlats(90)
!                 print*, 'disobs L367_mype1: rlons(1),rlons(2),rlons(90)=',rlons(1),rlons(2),rlons(90)
                 call tintrp2a11_indx(dlat,dlon,dtime,hrdifsig,&
                                      mype,nfldsig,ix,ixp,iy,iyp,it,itp)

                 print*, "disobs L367: anal_grid_subdomain, here we will use only it and itp info for coef_inf array: it,itp, ix,iy,mm1=",it,itp,ix,iy,mm1
                 print*, "disobs L368: ------1a) anal_grid_subdomain ---- it and itp info=",it,itp
                 !---------------------------------------------------
                 !    1b) ENSEMBLE GRID -> to get ix_ens_global(index of model lat/lon grid (4 points) w.r.t. obs location)
                 !                         for inflating and vertical localization
                 !---------------------------------------------------
                 call gsi_enscoupler_localization_grid (rlats_ens_local,rlons_ens_local)
                 !print*, 'disobs ens grid:L374: rlats_ens_local,rlons_ens_local=',rlats_ens_local,rlons_ens_local
                 ! global domain: rlats_ens_local, rlons_ens_local
                 ! rlats_ens_local = -1.57, ..., +1.57 (radian)
                 ! rlons_ens_local = 0, ..., +6.261368 (radian)

                 ! Find grid point near lat_earth, lon_earth on ens grid (global)
                 lon_deg=obs_data(12,n) ! deg
                 lat_deg=obs_data(13,n) ! deg
                 if(lon_deg>= r360)lon_deg=lon_deg-r360
                 if(lon_deg < zero)lon_deg=lon_deg+r360
                 dlon_earth_deg=lon_deg
                 dlat_earth_deg=lat_deg
                 dlon_earth=lon_deg*deg2rad ! radians
                 dlat_earth=lat_deg*deg2rad ! radians

                 dlon_ens = dlon_earth
                 dlat_ens = dlat_earth
                 call grdcrd1(dlon_ens,rlons_ens_local,grd_ens%nlon,1)
                 call grdcrd1(dlat_ens,rlats_ens_local,grd_ens%nlat,1)
                 ! dlon_ens, dlat_ens: point converted to grid units (grid relative lon,lat) (global ens grid)
!                 print*, 'disobs ens L375: index in global ens grid: dlat_ens,dlon_ens,mm1=', dlat_ens,dlon_ens,mm1
!                 print*, 'disobs ens L376: rlats_ens_local(int(dlat_ens))*rad2deg,rlons_ens_local(int(dlon_ens))*rad2deg,mm1=', rlats_ens_local(int(dlat_ens))*rad2deg,rlons_ens_local(int(dlon_ens))*rad2deg,mm1
!                 print*, 'disobs ens L377: global ens grid degree lat,lon,mm1=', rlats_ens_local(int(dlat_ens))*rad2deg,rlons_ens_local(int(dlon_ens))*rad2deg,mm1
                 !ix_ens_global,iy_ens_global - horizontal grid indexes on global ens grid
                 ix_ens_global =int(dlat_ens)
                 ixp_ens_global=int(dlat_ens)+1
                 iy_ens_global =int(dlon_ens)
                 iyp_ens_global=int(dlon_ens)+1
                 print*, 'disobs ens L378: ----- 1b) global ens grid ----- mm1,ix_ens_global,ixp,iy,iyp=', mm1,ix_ens_global,ixp_ens_global,iy_ens_global,iyp_ens_global
                 !call tintrp2a11_indx_ens_grd(dlat_ens,dlon_ens,dtime,hrdifsig,&
                 !                     mype,nfldsig,ix_ens,ixp_ens,iy_ens,iyp_ens,it,itp)
                 !print*, 'disobs ens L377: mm1,ix_ens,ixp_ens,iy_ens,iyp_ens=',mm1,ix_ens,ixp_ens,iy_ens,iyp_ens

                 !-------------------------------------------------------
                 ! 2) Obtain the correnponding model level to PBLH
                 !    geop_hgti, ges_geopi: (public) guess geopot height at level interfaces
                 !    load_geop_hgt is called in read_guess.F90
                 !-------------------------------------------------------
                 !    2a) Interpolate to get ges_geopi at obs location/time -> output: geopiges(nsig+1)
                 allocate(geopiges(nsig+1))
                 print*, 'disobs L381: dlat,dlon,dtime,mm1=',dlat,dlon,dtime,mm1
                 print*, 'disobs L382: mm1,ges_geopi(4,4,10,1)=',mm1,ges_geopi(4,4,10,1)
                 call tintrp2a1(ges_geopi,geopiges,dlat,dlon,dtime,hrdifsig,&
                               (nsig+1),mype,nfldsig)
                 print*,"disobs L383: mm1=",mm1,",geopiges",geopiges

                 !----------------------------------------------------------
                 !    2b) Find model level near obs pblh
                 !        depending on pblh type (1: pblri, 2: q-gradient pblh)
                 !----------------------------------------------------------

                 !    2b_1) convert log(obs pblh) to obs pblh
                 !if (data(ipblri,i) == -9999.0) then
                 !if (obs_data(5,n) == -999.0 .or. obs_data(5,n) == -9999.0) then
                 if (obs_data(5,n) < 0 ) then
                    print*,"disobs L390:obs_data(5,n)=",obs_data(5,n)
                    pblh_obs=-9999.0
                 else
                    !----------------------------------------------------------
                    ! inflate ens spread near pblh 
                    ! 1) for pblri -> q-gradient pblh (or pblri)
                    ! 1) for pblrf -> pblhlow
                    ! 1) for pblkh -> pblhhgh
                    !----------------------------------------------------------
                    if (trim(obstype)=='pblri') then ! pblri
                    
                       if (inf_ensprd_pblri_type==1) then ! pblri
                          !pblh_obs=exp(obs_data(5,n)) 
                          pblh_obs=obs_data(5,n) 
                       else if (inf_ensprd_pblri_type==2) then ! q-gradient pblh
                          pblh_obs=obs_data(16,n) ! [meters] 
                       end if
                       
                    else if (trim(obstype)=='pblrf') then ! pblrf (gnssro pblhlow)
                    
                       !pblh_obs=exp(obs_data(5,n)) 
                       pblh_obs=obs_data(5,n) 

                    else if (trim(obstype)=='pblkh') then ! gnssro pblhhgh
                       pblh_obs=obs_data(5,n)

                    end if

                 end if

!                 print*,"disobs L395: mm1,obs_data(5,n)=",mm1,obs_data(5,n)
!                 if ( trim(obstype)=='pblri' .and. inf_ensprd_pblri) then
!                    print*,"disobs L396: mm1,inf_ensprd_pblri_type for pblri (1=pblri, 2=qgrad_pblh)=",inf_ensprd_pblri_type
!                 else if ( trim(obstype)=='pblrf' .and. inf_ensprd_pblrf) then
!                    print*,"disobs L396: mm1, inflate ens spread near pblrf"
!                 end if
!                 print*,"disobs L397: mm1,obstype,obs_data(5,n),pblh_obs=",mm1,trim(obstype),obs_data(5,n),pblh_obs

                 !    2b_2) find model level near obs pblh (model geopotential height (geopiges) is interpolated to obs location/time)
                 do k = 1, nsig ! bottom to top
                    if (pblh_obs>=geopiges(k) .and. pblh_obs<geopiges(k+1)) then
                       print*, "dobs L401:mm1=",mm1,",geopiges(k),(k+1)=",geopiges(k),geopiges(k+1)
                       print*, "dobs L402:obstype=",obstype
                       print*, "mm1,obs_data(5,n)=",mm1,obs_data(5,n)
                       print*, "mm1,pblh_obs=",mm1,pblh_obs,",k=",k
                       iz=k
                       izp=k+1
                       exit
                    end if
                 end do

                 if (pblh_obs < 0) then ! missing value
                    iz=-9999
                    izp=-9999
                 end if

                 deallocate(geopiges)

                 !===============================================================
                 ! 3) (Iflating ens spread and Vertical localization) 
                 !    Save grid information of pblh obs on ens global grid
                 !    to modify ens spread and hybrid vertical localization length scale near PBLH obs
                 !===============================================================

                 !----------------------------------------------------
                 !   3a) ensemble inflation coefficient (x,y,z,t)
                 !     (ix_ens_global, ixp_ens_global, iy_ens_global, iyp_ens_global, iz, izp, it, itp)
                 !     -> initialized value is 1.0
                 !----------------------------------------------------
                 !    Assign coefficient for inflating ens spread near obs pblh
                 !    at the adjacent 3 mid-layer model levels
                 !    at the adjacent 5 mid-layer model levels
                 !    ---------------  4
                 !    - - - - - - - -
                 !    ---------------  3
                 !    - - * - - - - -
                 !    ---------------  2
                 !    - - - - - - - -
                 !    ---------------  1
                 !    if obs is between interface 2 and 3, then we will assign coefficient for
                 !    mid-layer model level 1, 2, and 3 (interface 1,2,3,4)
                 !    ges_coef_inf(lat2,lon2,nsig,nfldsig)
                 !    if iz is lowest or izp is highest model, then consider 2 levels
                 !    ---------------------------
                 !    3a_1) for three adjacent vertical model levels
                 !    ---------------------------
                 !     if (iz == 1) then
                 !        ges_coef_inf(ix:ixp, iy:iyp, iz:izp, it:itp) = 20.0 ! 5.0
                 !     else if (iz == nsig) then
                 !        ges_coef_inf(ix:ixp, iy:iyp, (iz-1):iz, it:itp) = 20.0 !5.0
                 !     else
                 !        ges_coef_inf(ix:ixp, iy:iyp, (iz-1):izp, it:itp) = 20.0 !5.0
                 !     end if
                 !    ---------------------------
                 !    3a_2) for five adjacent vertical model levels
                 !    ---------------------------
                 !     if (iz == 1) then ! 3 levels
                 !        ges_coef_inf(ix:ixp, iy:iyp, iz:(iz+2), it:itp) = 20.0 ! 5.0
                 !     else if (iz == 2) then ! 4 levels
                 !        ges_coef_inf(ix:ixp, iy:iyp, (iz-1):(iz+2), it:itp) = 20.0 ! 5.0
                 !     else if (iz == nsig) then ! 3 levels
                 !        ges_coef_inf(ix:ixp, iy:iyp, (iz-2):iz, it:itp) = 20.0 !5.0
                 !     else if (iz == (nsig-1)) then ! 4 levels
                 !        ges_coef_inf(ix:ixp, iy:iyp, (iz-2):(iz+1), it:itp) = 20.0 !5.0
                 !     else ! 5 levels
                 !        ges_coef_inf(ix:ixp, iy:iyp, (iz-2):(iz+2), it:itp) = 20.0 !5.0
                 !     end if

                 !    ---------------------------
                 !    3a_3) for five adjacent vertical model levels (Gaussian dist.)
                 !    ---------------------------
                 !!!!!!!!!!!!!!!!!!!!
                 sd2=1.0_r_kind ! variance = square of standard deviation 
                 coeff_t=inf_coef_t ! coefficient for 5 times inflation with sd=1 for T
                 coeff_q=inf_coef_q ! coefficient for 10 times inflation with sd=1 for q
                 !coeff_t=10.0_r_kind ! coefficient for 5 times inflation with sd=1 for T
                 !coeff_q=22.5_r_kind ! coefficient for 10 times inflation with sd=1 for q
                 !lev_diff = 0.0, 1.0, 2.0 : model level - model level closest to pblh
                 ! 0) at model level closest to PBLH
                 lev_diff = 0.0_r_kind
                 !inf_coef0-2: inflation coefficient (1+pdf*coef)
                 inf_coef0_t = 1.0_r_kind + coeff_t*1.0_r_kind/( sqrt(2.0_r_kind*pi*sd2) ) * exp( -0.5_r_kind*(lev_diff)**2.0_r_kind/sd2)
                 inf_coef0_q = 1.0_r_kind + coeff_q*1.0_r_kind/( sqrt(2.0_r_kind*pi*sd2) ) * exp( -0.5_r_kind*(lev_diff)**2.0_r_kind/sd2)
                 !inf_coef0 ~ 5 
                 !inf_coef0 ~ 5 or 10.0 or 30
 
                 ! 1) +/- 1 from model level closest to PBLH
                 lev_diff = 1.0_r_kind
                 inf_coef1_t = 1.0_r_kind + coeff_t*1.0_r_kind/( sqrt(2.0_r_kind*pi*sd2) ) * exp( -0.5_r_kind*(lev_diff)**2.0_r_kind/sd2)
                 inf_coef1_q = 1.0_r_kind + coeff_q*1.0_r_kind/( sqrt(2.0_r_kind*pi*sd2) ) * exp( -0.5_r_kind*(lev_diff)**2.0_r_kind/sd2)
                 !inf_coef1 ~ 3.45
                 !inf_coef1 ~ 2.8 or 5.0 or 14
 
                 ! 2) +/- 2 from model level closest to PBLH
                 lev_diff = 2.0_r_kind
                 inf_coef2_t = 1.0_r_kind + coeff_t*1.0_r_kind/( sqrt(2.0_r_kind*pi*sd2) ) * exp( -0.5_r_kind*(lev_diff)**2.0_r_kind/sd2)
                 inf_coef2_q = 1.0_r_kind + coeff_q*1.0_r_kind/( sqrt(2.0_r_kind*pi*sd2) ) * exp( -0.5_r_kind*(lev_diff)**2.0_r_kind/sd2)
                 !inf_coef2 ~ 1.55
                 !inf_coef2 ~ 1.17 or 1.5 or 2
 
                 if (mm1==216) print*, 'disobs L555, inf_coef_t=',inf_coef0_t,inf_coef1_t,inf_coef2_t,',inf_coef_q=',inf_coef0_q,inf_coef1_q,inf_coef2_q
                 !print*, 'disobs L556, inf_ens_sprd: inf_coef1=',inf_coef1
                 !print*, 'disobs L557, inf_ens_sprd: inf_coef2=',inf_coef2
!                 print*, 'disobs L558, inf_ens: mm1,ix_ens,iy,iz,it=',mm1,ix_ens_global,iy_ens_global,iz,it
!                 print*, 'disobs L559, inf_ens: mm1=',mm1
                 !grid_coef_inf_ens_grid_global(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz, it:itp) = 1.0_r_kind

                 if (iz > 0) then ! 5 levels (to exclude missing value (iz=-9999))

                    ! t
                    !ges_coef_inf_ens_grid_global_t(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz+2, it:itp) = inf_coef2_t !1.55
                    !ges_coef_inf_ens_grid_global_t(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz+1, it:itp) = inf_coef1_t !3.45
                    !ges_coef_inf_ens_grid_global_t(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz, it:itp) = inf_coef0_t   !5
                    !ges_coef_inf_ens_grid_global_t(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz-1, it:itp) = inf_coef1_t !3.45
                    !ges_coef_inf_ens_grid_global_t(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz-2, it:itp) = inf_coef2_t !1.55
                    ges_coef_inf_ens_grid_global_t(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz+2, it:itp) = & 
                    max(ges_coef_inf_ens_grid_global_t(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz+2, it:itp), inf_coef2_t) !1.55
                    ges_coef_inf_ens_grid_global_t(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz+1, it:itp) = & 
                    max(ges_coef_inf_ens_grid_global_t(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz+1, it:itp), inf_coef1_t) !3.45
                    ges_coef_inf_ens_grid_global_t(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz, it:itp) = &
                    max(ges_coef_inf_ens_grid_global_t(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz, it:itp), inf_coef0_t)   !5
                    ges_coef_inf_ens_grid_global_t(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz-1, it:itp) = &
                    max(ges_coef_inf_ens_grid_global_t(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz-1, it:itp), inf_coef1_t) !3.45
                    ges_coef_inf_ens_grid_global_t(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz-2, it:itp) = &
                    max(ges_coef_inf_ens_grid_global_t(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz-2, it:itp), inf_coef2_t) !1.55

                    ! q
                    !ges_coef_inf_ens_grid_global_q(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz+2, it:itp) = inf_coef2_q !2.23
                    !ges_coef_inf_ens_grid_global_q(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz+1, it:itp) = inf_coef1_q !6.45
                    !ges_coef_inf_ens_grid_global_q(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz, it:itp) = inf_coef0_q   !10
                    !ges_coef_inf_ens_grid_global_q(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz-1, it:itp) = inf_coef1_q !6.45
                    !ges_coef_inf_ens_grid_global_q(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz-2, it:itp) = inf_coef2_q !2.23
                    ges_coef_inf_ens_grid_global_q(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz+2, it:itp) = &
                    max(ges_coef_inf_ens_grid_global_q(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz+2, it:itp), inf_coef2_q) !2.23
                    ges_coef_inf_ens_grid_global_q(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz+1, it:itp) = &
                    max(ges_coef_inf_ens_grid_global_q(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz+1, it:itp), inf_coef1_q) !6.45
                    ges_coef_inf_ens_grid_global_q(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz, it:itp) = &
                    max(ges_coef_inf_ens_grid_global_q(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz, it:itp), inf_coef0_q)   !10
                    ges_coef_inf_ens_grid_global_q(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz-1, it:itp) = &
                    max(ges_coef_inf_ens_grid_global_q(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz-1, it:itp), inf_coef1_q) !6.45
                    ges_coef_inf_ens_grid_global_q(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz-2, it:itp) = &
                    max(ges_coef_inf_ens_grid_global_q(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz-2, it:itp), inf_coef2_q) !2.23

                 end if
 
                 !----------------------------------------------------
                 !   3b) vertical localzation (x,y,z) no time info required
                 !     (ix_ens_global, ixp_ens_global, iy_ens_global, iyp_ens_global, iz)
                 !     -> initialized value is 0.0
                 !----------------------------------------------------
                 !ges_vlocal(ix:ixp, iy:iyp, iz,1) = 1.0_r_kind
                 ! 4 Surrounding ens grid points at one level (near pblh)
                 if (iz > 0) then ! only for non-missing value (missing obs -> iz = -9999)
                    ges_vlocal_ens_grid_global(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz) = 1.0_r_kind
                    !print*, 'dobs vlocal L592: mm1,iz,izp=',mm1,iz,izp
                    !print*, 'dobs vlocal L593: mm1,ges_vlocal_ens_grid_global(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz:izp)=',mm1,ges_vlocal_ens_grid_global(ix_ens_global:ixp_ens_global, iy_ens_global:iyp_ens_global, iz:izp)
                 end if

              end if
           end if
           !------
           prec_loop: do k=1,mm1-1
              if(lat>=ibs(k).and.lat<=ibn(k)) then
                 if((lon >= ibw(k).and. lon <=ibe(k))  .or.  &
                    (lon == 0     .and. ibe(k) >=nlon) .or.  &
                    (lon == nlon  .and. ibw(k) <=1) .or. periodic_s(k)) then
                    luse_s(ndata_s)= .false.
                    exit prec_loop
                 end if
              end if
           end do prec_loop
           if(nobs == ndata_s) exit obs_loop
         end if
     end if
  end do obs_loop

  deallocate(obs_data)

! Write observations for given task to output file
  write(lunout) obstypeall,isis,nreal,nchanl
  write(lunout) data1_s,luse_s,ioid_s
  deallocate(data1_s,luse_s,ioid_s)

  return
end subroutine disobs
subroutine count_obs(ndata,nn_obs,lat_data,lon_data,obs_data,nobs_s)

!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    count_obs  counts number of observations on each subdomain
!   prgmmr: derber                                    date: 2014-12-16
!
! abstract: count observations in each pe subdomain
!
! program history log:
!   2014-12-16  derber
!
!   input argument list:
!     ndata    - number of observations
!     nn_obs   - number of an observation data for each ob 
!     lat_data - location of lattitude
!     lon_data - location of longitude
!     obs_data - observation information array
!
!   output argument list:
!     nobs_s   - number of observations on all subdomains
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
  use gridmod, only: periodic_s,nlon,nlat,jlon1,ilat1,istart,jstart
  use mpimod, only: npe
  implicit none

! Declare passed variables
  integer(i_kind)               ,intent(in   ) :: ndata,lat_data,lon_data
  integer(i_kind)               ,intent(in   ) :: nn_obs
  integer(i_kind),dimension(npe),intent(inout) :: nobs_s
  real(r_kind),dimension(nn_obs,ndata),intent(in) :: obs_data

! Declare local variables
  integer(i_kind) lon,lat,n,k
  integer(i_kind),dimension(npe):: ibe,ibw,ibn,ibs

! Read and write header

  do k=1,npe

!    ibw,ibe,ibs,ibn west,east,south and north boundaries of total region
     ibw(k)=jstart(k)-1
     ibe(k)=jstart(k)+jlon1(k)-1
     ibs(k)=istart(k)-1
     ibn(k)=istart(k)+ilat1(k)-1

  end do


  nobs_s=0

! Loop over all observations.  Locate each observation with respect
! to subdomains.
  do n=1,ndata
     lat=obs_data(lat_data,n)
     lat=min(max(1,lat),nlat)

     do k=1,npe
        if(lat>=ibs(k).and.lat<=ibn(k)) then
           lon=obs_data(lon_data,n)
           lon=min(max(0,lon),nlon)
           if((lon >= ibw(k).and. lon <=ibe(k))  .or.  &
              (lon == 0   .and. ibe(k) >=nlon) .or.  &
              (lon == nlon    .and. ibw(k) <=1) .or. periodic_s(k)) then
                 nobs_s(k)=nobs_s(k)+1
           end if
        end if
     end do
  end do 
  return
end subroutine count_obs

! ------------------------------------------------------------------------
subroutine dislag(ndata,mm1,lunout,obsfile,obstypeall,ndata_s)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    dislag  distribute lagrangian observations into each pe 
!                subdomain
!   prgmmr: lmeunier                                  date: 2009-03-12
!
! abstract: distribute lagrangian observations into each pe subdomain
!           (based on disobs). All observations of one balloon are
!           are associated with a given processor
!
! program history log:
!   2009-03-12  lmeunier
!
!   input argument list:
!     ndata_s  - number of observations in each pe sub-domain
!     ndata    - number of observations
!     obsfile  - unit from which to read all obs of a given type
!     lunout   - unit to which to write subdomain specific observations
!     mm1      - mpi task number + 1
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
  use lag_fields, only: orig_lag_num
  implicit none

! Declare passed variables
  integer(i_kind),intent(in   ) :: ndata,lunout,mm1
  integer(i_kind),intent(inout) :: ndata_s
  character(14)  ,intent(in   ) :: obsfile
  character(10)  ,intent(in   ) :: obstypeall

! Declare local variables
  integer(i_kind) num_data,n,lunin,num
  integer(i_kind) jj,nreal,nchanl,nn_obs,ndatax
  logical,allocatable,dimension(:):: luse,luse_s,luse_x
  real(r_kind),allocatable,dimension(:,:):: obs_data,data1_s
  character(10):: obstype
  character(20):: isis

  ndata_s=0

  lunin=11
  open(lunin,file=obsfile,form='unformatted')
  read(lunin)obstype,isis,nreal,nchanl,num_data
  if(obstype /=obstypeall) &
        write(6,*)'DISLAG:  ***ERROR***   obstype,obstypeall=',obstype,obstypeall

  nn_obs = nreal + nchanl

  allocate(obs_data(nn_obs,ndata))
! Read in all observations of a given type along with subdomain flags
  read(lunin) obs_data
  close(lunin)

  allocate(luse(ndata),luse_x(ndata))
  luse=.false.
  luse_x=.false.

! Loop over all observations.  Locate each observation with respect
! to subdomains.
  use: do n=1,ndata

!    Does the observation belong to the subdomain for this task?
     num=int(obs_data(num_data,n),i_kind)  
     if ((mm1-1)==orig_lag_num(num,2)) then
        ndata_s=ndata_s+1
        luse(n)=.true.
        luse_x(n)=.true.  ! Never on another subdomain
     end if

  end do use

  if(ndata_s > 0)then
     allocate(data1_s(nn_obs,ndata_s),luse_s(ndata_s))
     ndatax=0
     do n=1,ndata

        if(luse(n))then

           ndatax=ndatax+1
           luse_s(ndatax)=luse_x(n)
 
           do jj= 1,nn_obs
              data1_s(jj,ndatax) = obs_data(jj,n)
           end do

        end if

     end do

!    Write observations for given task to output file
     write(lunout) obstypeall,isis,nreal,nchanl
     write(lunout) data1_s,luse_s
     deallocate(data1_s,luse_s)
  endif
  deallocate(obs_data,luse,luse_x)

  return
end subroutine dislag
