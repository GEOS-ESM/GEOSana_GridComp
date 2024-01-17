module tgas_setup
  implicit none
  private
  public :: setup
  interface setup; module procedure setuptgas; end interface

contains
subroutine setuptgas(obsLL, odiagLL, lunin, mype, stats_tgas, nchanl, nreal,   &
                     nobs, obstype, isis, is, ltgasdiagsave, init_pass)
!$$$  subprogram documentation block
!                .      .    .
! subprogram:    setuptgas --- Compute rhs of oi for trace gas obs
!    prgrmmr:    weir             org: gmao                date: 2014-04-14
!
!   abstract:    For trace gas observations, this routine
!                  a) reads obs assigned to given mpi task (geographic region),
!                  b) simulates obs from guess,
!                  c) applies some quality control to obs,
!                  d) loads weight and innovation arrays used in minimization
!                  e) collects statistics for runtime diagnostic output
!                  f) writes additional diagnostic information to output file
!
! program history log:
!   2014-03-30   weir  - initial code
!   2014-07-24   weir  - changed to interpolation of prior and averaging
!                        kernel to model levels instead of interpolation of
!                        guess to averaging kernel levels
!   2014-11-05   weir  - fixed some bugs with the diag output
!   2014-11-14   weir  - fixed a bug with how pressure level at surface was
!                        computed for mopitt; changed how extrapolation is
!                        done at the top of the atmosphere in zintrp, this
!                        needs more work
!   2015-07-01   weir  - removed all the changes with levels, etc. above, and
!                        made everything follow ozone more closely
!
!
!   input argument list:
!     obsLL         - a linked-list of obsNode holding observations
!     odiagLL       - a linked-list of obs_diag holding obsdiags
!     lunin         - unit from which to read observations
!     mype          - mpi task id
!     nchanl        - the number of elements of the obs vector; 1 for column
!                     obs; nlays for layer obs; nlays+1 for interface obs
!     nreal         - number of pieces of non-data info (location, time, etc.)
!                     per obs
!     nobs          - number of observations
!     isis          - sensor/instrument/satellite id
!     is            - integer counter for number of obs types to process
!     obstype       - type of tgas obs
!     ltgasdiagsave - switch on diagnostic output (.false. = no output)
!     init_pass     - state of setup processing
!     stats_tgas    - sums for various statistics as a function of data input
!
!   output argument list:
!     obsLL         - a linked-list of obsNode holding observations
!     odiagLL       - a linked-list of obs_diag holding obsdiags
!     stats_tgas    - sums for various statistics as a function of data input
!
!   notes about function of use variables:
!     luse  - processor use flag, set by obspara
!                      = .false., obs not assigned to this processor
!                      = .true.,  obs     assigned to this processor
!     qcuse - gross and qc check use flag, set by this routine
!                      = 0, obs does not pass gross and qc checks
!                      = 1, obs does     pass gross and qc checks
!     iouse - channel use flag, set by tgasinfo
!                      = -2, do not use
!                      = -1, monitor if diagnostics produced
!                      =  0, monitor and use in QC only
!                      =  1, use data with complete quality control
!     (the below are used by rad, but not implemented here)
!                      =  2, use data with no airmass bias correction
!                      =  3, use data with no angle dependent bias correction
!                      =  4, use data with no bias correction
!
!$$$ end documentation block

! Declare module variables
  use m_obsdiagNode, only: obs_diag
  use m_obsdiagNode, only: obs_diags
  use m_obsdiagNode, only: obsdiagLList_nextNode
  use m_obsdiagNode, only: obsdiagNode_set
  use m_obsdiagNode, only: obsdiagNode_get
  use m_obsdiagNode, only: obsdiagNode_assert

  use m_obsNode,  only: obsNode
  use m_tgasNode, only: tgasNode
  use m_tgasNode, only: tgasNode_appendto
  use m_obsLList, only: obsLList

  use nc_diag_write_mod, only: nc_diag_init, nc_diag_header, nc_diag_metadata, &
                               nc_diag_write, nc_diag_data2d
  use nc_diag_read_mod,  only: nc_diag_read_init, nc_diag_read_get_dim,        &
                               nc_diag_read_close

  use mpeu_util,         only: die
  use kinds,             only: r_kind, r_single, i_kind
  use constants,         only: zero, half, one, two, tiny_r_kind, huge_r_kind, &
                               cg_term, wgtlim, varlen => max_varname_length
  use constants,         only: vv_to_molec
  use obsmod,            only: dplat, nobskeep, mype_diaghdr, dirname,         &
                               time_offset, ianldate, luse_obsdiag,            &
                               lobsdiag_allocated, lobsdiagsave, netcdf_diag,  &
                               binary_diag, wrtgeovals
  use gsi_4dvar,         only: nobs_bins, mn_obsbin
  use gridmod,           only: get_ij, nsig
  use guess_grids,       only: nfldsig, ntguessig, hrdifsig, ges_pe => ges_prsi
  use guess_grids,       only: no2_priorc, so2_priorc, pbltopl
  use gsi_bundlemod,     only: gsi_bundlegetpointer
  use gsi_chemguess_mod, only: gsi_chemguess_get, gsi_chemguess_bundle
  use tgasinfo,          only: jpch_tgas, nusis_tgas, iuse_tgas, gross_tgas,   &
                               error_tgas, b_tgas, pg_tgas, bias_tgas,         &
                               vname_tgas, vunit_tgas, ntgas, tgnames,         &
                               nreact, reactname, tgas_minsigobs,              &
                               tgas_maxsigobs, tgas_sigobsscal
  use jfunc,             only: jiter, last, miter
  use m_dtime,           only: dtime_setup, dtime_check

  implicit none

! Declare passed variables
  type(obsLList),  target, dimension(:), intent(inout) :: obsLL     ! a L-List array of tgasNode
  type(obs_diags), target, dimension(:), intent(inout) :: odiagLL   ! a L-List array of obs_diag

  integer(i_kind),   intent(in   ) :: lunin, mype, nchanl, nreal, nobs, is
  character(len=20), intent(in   ) :: isis
  character(len=*),  intent(in   ) :: obstype
  logical,           intent(in   ) :: ltgasdiagsave, init_pass
  real(r_kind),      intent(inout) :: stats_tgas(15,jpch_tgas)

! Declare local parameters  
  integer(i_kind),  parameter :: iint  = 1
  integer(i_kind),  parameter :: ireal = 3
  integer(i_kind),  parameter :: irmin = 6
  real(r_kind),     parameter :: tol_error = 1.e4_r_kind
  character(len=*), parameter :: myname    = 'setuptgas'
  logical,          parameter :: ldebug    = .false.
  logical,          parameter :: lcenter   = .false.
  real(r_kind),     parameter :: rmiss     = -9999.9_r_kind

! Declare local variables  
  real(r_kind) :: obs, bias, omg, soundn, grdtime
  real(r_kind) :: grdlat, grdlon, deglat, deglon
  real(r_kind) :: cg_tgas, wgross, wgt, arg, exp_arg, term
  real(r_kind) :: errorinv, fact, delz

  real(r_kind),    dimension(nchanl) :: obsinv, varinv, raterr2
  real(r_kind),    dimension(nchanl) :: pchanl, gross, sclstd, sclmod
  real(r_kind),    dimension(nchanl) :: priorobs, gesobs, uncert, error
  real(r_kind),    dimension(nchanl) :: ratio_errors
  real(r_kind),    dimension(nsig+1) :: pemod, zemod, peuse
  real(r_kind),    dimension(nsig)   :: geslmod, qvmod, dpmod, plmod, pluse
  integer(i_kind), dimension(nchanl) :: isbad

  real(r_kind),    dimension(nreal+nchanl,nobs) :: tgasdata

  real(r_single),  dimension(ireal,nobs) :: diagbuf

  integer(i_kind), allocatable, dimension(:)     :: itgas
  real(r_kind),    allocatable, dimension(:)     :: peavg, grdpe, priorpro
  real(r_kind),    allocatable, dimension(:)     :: peges, dpges
  real(r_kind),    allocatable, dimension(:)     :: gespro, grdpch
  real(r_kind),    allocatable, dimension(:,:)   :: avgker, avgwgt
  real(r_single),  allocatable, dimension(:,:,:) :: rdiagbuf
  real(r_kind),    allocatable, dimension(:,:,:,:)   :: ges, ges_ze, ges_qv
  real(r_kind),    allocatable, dimension(:,:,:,:,:) :: ges_tgas

  integer(i_kind) :: nchon, navg, nedge
  integer(i_kind) :: jj, ibin, itg, iz, izp
  integer(i_kind) :: jc, idout, jdiag, irdim1, ier, istatus
  integer(i_kind) :: ikeep, mm1, fsign
  integer(i_kind) :: inum, itime, iglat, iglon, idlat, idlon, iread
  integer(i_kind) :: i, j, k, l, n

  integer(i_kind), dimension(iint,nobs) :: idiagbuf
  integer(i_kind), dimension(nchanl)    :: ipos, iouse, qcuse

  character(len=12)     :: string
  character(len=10)     :: filex
  character(len=128)    :: diag_tgas_file
  character(len=varlen) :: vname, vunit
  logical               :: tgasdiagexist

  logical :: lmaybepassive, lmaybeany, ldebugob
  logical :: in_curbin, in_anybin, useak

  logical, dimension(nobs) :: luse

  integer(i_kind), dimension(nobs) :: ioid  ! initial (pre-distribution) observation ID

  ! stuff for NO2/SO2 (call it 'DOAS' here)
  logical                                        :: isdoas
  real(r_kind),    allocatable, dimension(:,:)   :: tmpwgt
  real(r_kind),    allocatable, dimension(:)     :: qvavg
  real(r_kind),    allocatable, dimension(:)     :: delps
  real(r_kind),    dimension(nsig)               :: priortgas
  real(r_kind)                                   :: scd_tot, scd_trop, amf, amftrop, pbl, iuncert
  integer(i_kind)                                :: ipbl
  real(r_kind)                                   :: lscal, surfscal
  real(r_kind)                                   :: minsigobs, maxsigobs, sigobsscal
  logical                                        :: ispos
  real(r_kind), dimension(:,:,:), pointer        :: rank3 => NULL()

  type(tgasNode),  pointer :: my_head
  type(obs_diag),  pointer :: my_diag
  type(obs_diags), pointer :: my_diagLL
  type(obsLList),  pointer, dimension(:) :: tgashead
  tgashead => obsLL(:)

! Check if this is an averaging-kernel obs type
  isdoas = ( trim(obstype) == 'omno2' .or. trim(obstype) == 'mindsno2' .or. &
             trim(obstype) == 'omso2' .or. trim(obstype) == 'nmso2'    .or. &
             trim(obstype) == 's5pno2' )
  useak  = ( trim(obstype) == 'tgav'  .or. trim(obstype) == 'tgaz'     .or. &
             trim(obstype) == 'acos'  .or. isdoas )

! Determine level variables navg and nedge from nreal and nchanl
  call init_levs_(useak, nreal, nchanl, navg, nedge)

! If required guess fields are available, extract from bundle
  call init_vars_

  mm1 = mype + 1


!$$$ -- 1. ALLOCATE AND INITIALIZE VARIABLES --
!===============================================================================
  do k = 1,nchanl
     ipos(k)   = 0
     iouse(k)  = -2
     gross(k)  = huge_r_kind
     sclstd(k) = huge_r_kind
     sclmod(k) = 1.e+00_r_kind
  end do
  vname = 'null'
  vunit = 'null'

! Locate data for satellite in tgasinfo arrays, and remove channels if
! necessary
  lmaybepassive = .false.
  jc = 0
  if (ldebug) print *, 'jpch_tgas = ', jpch_tgas, 'nchanl = ', nchanl
  do j = 1,jpch_tgas
     if (trim(isis) == trim(nusis_tgas(j))) then
        if (ldebug) print *, 'j = ', j, 'isis = ', isis

        jc = jc + 1
        if (nchanl < jc) then
           write(0,*) myname // ': *** ERROR *** in channel number'
           write(0,*) myname // ': jc    = ', jc, ', nchanl = ', nchanl
           call die(myname)
        end if

        ipos(jc)   = j
        iouse(jc)  = iuse_tgas(j)
        gross(jc)  = gross_tgas(j)
        sclstd(jc) = error_tgas(j)

!       The constituent and its units are restricted to be the same for all
!       channels for averaging kernel type data (this restriction could be
!       relaxed, but would complicate the code).  Sample measurements can
!       have mixed constituents and types
        if (useak .and. 1 < jc) then
           if (trim(vname_tgas(j)) /= vname) then
              write(0,*) myname // ': all channels of averaging kernel ' //    &
                         'type obs must be of the same constituent'
              call die(myname)
           end if

           if (trim(vunit_tgas(j)) /= vunit) then
              write(0,*) myname // ': all channels of averaging kernel ' //    &
                         'type obs must have the same units'
              call die(myname)
           end if
        end if

        vname = trim(vname_tgas(j))
        vunit = trim(vunit_tgas(j))

!       Set every time (inefficient) to keep code in same place
        if (vunit == 'ppd') sclmod(jc) = 1.e+01_r_kind
        if (vunit == 'ppc') sclmod(jc) = 1.e+02_r_kind
        if (vunit == 'ppk') sclmod(jc) = 1.e+03_r_kind
        if (vunit == 'ppm') sclmod(jc) = 1.e+06_r_kind
        if (vunit == 'ppb') sclmod(jc) = 1.e+09_r_kind
        if (vunit == 'ppt') sclmod(jc) = 1.e+12_r_kind

        if (-1 < iouse(jc)) lmaybepassive = .true.
        if (-2 < iouse(jc)) lmaybeany     = .true.

        if (ldebug) then
           print *, 'iouse  = ', iouse(jc)
           print *, 'gross  = ', gross(jc)
           print *, 'sclstd = ', sclstd(jc), 'sclmod = ', sclmod(jc)
           print *, 'ltgasdiagsave = ', ltgasdiagsave
        end if
     end if
  end do
  nchon = jc

! If there's no data for assimilation or diagnostic output, return
  if (.not. lmaybeany .or. nchon == 0) then
     if (mype == 0) write(0,*) myname // ': no channels found for ', isis
     if (0  < nobs) read(lunin)

     call final_vars_
     return
  end if

! Handle error conditions
  if (nchon < nchanl) then
!    Unsure how to properly implement this, so reporting as an error
!    Instead, turn off channels with the iuse flag
!    if (mype == 0) write(0,*) myname // ': level number reduced for ',        &
!                              obstype, ' ', nchanl,' --> ', nchon
     write(0,*) myname // ': *** ERROR: Incompatible channel numbers ***'
     write(0,*) myname // ': nchon = ', nchon, ', nchanl = ', nchanl
     call die(myname)
  end if

! If requested, save data for diagnostic ouput
  if (ltgasdiagsave) then
     idout  = 0
     irdim1 = irmin
     if (lobsdiagsave) irdim1 = irdim1 + 4*miter + 1
     allocate(rdiagbuf(irdim1,nchanl,nobs))
     if(netcdf_diag) call init_netcdf_diag_
  end if

  if (ldebug) then
     print *, 'nreal     = ', nreal,     'nobs    = ', nobs
     print *, 'nchanl    = ', nchanl,    'navg    = ', navg
     print *, 'ntguessig = ', ntguessig, 'nfldsig = ', nfldsig
     print *, 'hrdifsig  = '
     print *, hrdifsig
     if (useak) print *, 'nedge     = ', nedge
  end if

! Allocate local variables
  allocate(itgas(nchanl), gespro(navg))
  allocate(avgker(nchanl,navg), avgwgt(navg,nsig))

  if (useak) then
     allocate(peavg(nedge), grdpe(nedge), priorpro(navg))
     allocate(peges(nedge), dpges(nedge-1))
  else
     allocate(grdpch(nchanl))
  end if

  ! NO2/SO2 stuff 
  if( isdoas ) then
     allocate(tmpwgt(navg,nsig),qvavg(navg))
     allocate(delps(navg))
  endif

  call dtime_setup


!$$$ -- 2. READ OBSERVATIONAL DATA --
!$$$       2A. ALL DATA TYPES
!===============================================================================
  read(lunin) tgasdata, luse, ioid ! luse and ioid come from obs_para

! Index information for data array (see read_tgas subroutine)
  inum  = 1   ! index of sounding id
  itime = 2   ! index of analysis relative obs time
  iglon = 3   ! index of grid relative obs location (x)
  iglat = 4   ! index of grid relative obs location (y)
  idlon = 5   ! index of earth relative longitude (degrees)
  idlat = 6   ! index of earth relative latitude (degrees)

  do i = 1,nobs
     soundn  = tgasdata(inum,i)
     grdtime = tgasdata(itime,i)
     grdlat  = tgasdata(iglat,i)
     grdlon  = tgasdata(iglon,i)
     deglat  = tgasdata(idlat,i)
     deglon  = tgasdata(idlon,i)

     ldebugob = (ldebug .and. i <= 3)

!    If debugging, allow centering of obs in space and time to eliminate
!    interpolation
     if (ldebug .and. lcenter) then
        if (mype == 0) write(0,*) myname // ': *** WARNING: centering obs ' // &
                                  'in time and moving to nearest grid '     // &
                                  'point ***'
        grdtime = 3._r_kind
        grdlat  = real(int(grdlat))
        grdlon  = real(int(grdlon))
     end if

!    Check which time bin of the current obs (assumes obs were screened to
!    fit in time window in read routine; largely perfunctory in standard
!    setups as far as I can tell)
     call dtime_check(grdtime, in_curbin, in_anybin)
     if (.not. in_anybin) cycle

if (in_curbin) then
     if (ltgasdiagsave .and. luse(i)) then
        idout = idout + 1
        idiagbuf(1,idout) = mype                            ! mpi task number
        diagbuf(1,idout)  = tgasdata(idlat,i)               ! lat (degree)
        diagbuf(2,idout)  = tgasdata(idlon,i)               ! lon (degree)
        diagbuf(3,idout)  = tgasdata(itime,i) - time_offset ! time (hours rel. to ana.)
     end if

     if (nobskeep > 0) then
        write(0,*) myname // ': nobskeep ', nobskeep
        call die(myname)
     end if

!    Interpolate guess edge pressure (ges_pe) to obs time and place (pemod)
     call tintrp2a1(ges_pe, pemod, grdlat, grdlon, grdtime, hrdifsig,          &
                    nsig+1, mype, nfldsig)

!    Convert pemod units from kPa to hPa
     pemod = pemod * 10._r_kind
     dpmod =  abs(pemod(2:nsig+1) - pemod(1:nsig))
     plmod = 0.5*(pemod(2:nsig+1) + pemod(1:nsig))

!    The same for edge geopotential heights (ges_ze) and water vapor (ges_qv)
     call tintrp2a1(ges_ze, zemod, grdlat, grdlon, grdtime, hrdifsig,          &
                    nsig+1, mype, nfldsig)
     call tintrp2a1(ges_qv, qvmod, grdlat, grdlon, grdtime, hrdifsig,          &
                    nsig,   mype, nfldsig)

!    Hack to use geometric height (m) instead of pressure for
!    tgez, tgaz, and NOAA ObsPack data
     peuse = pemod
     if (obstype == 'tgez' .or. obstype == 'tgaz') peuse = zemod
     if (obstype == 'tgop'                       ) peuse = zemod
     pluse = 0.5*(peuse(2:nsig+1) + peuse(1:nsig))

     if (ldebugob) then
        print *, 'soundn     = ', soundn
        print *, 'grdlat     = ', grdlat, 'grdlon     = ', grdlon
        print *, 'grdtime    = ', grdtime
        print *, 'pemod      = '
        print *, pemod
        print *, 'peuse      = '
        print *, peuse
     end if

!    Read pressure, qcflag, and uncertainty for each channel
     iread = idlat
     do k = 1,nchanl
        iread = iread + 1
        pchanl(k) = tgasdata(iread,i)
     end do
     do k = 1,nchanl
        iread = iread + 1
        isbad(k) = tgasdata(iread,i)
     end do
     do j = 1,nchanl
        iread = iread + 1
        uncert(j) = tgasdata(iread,i)
     end do

!$$$       2B. ADDITIONAL DATA FOR AVERAGING KERNEL RETRIEVALS
!===============================================================================
     if (useak) then
!       1. peavg - edge pressures/altitudes of averaging kernel levels
        do k = 1,nedge
           iread = iread + 1
           peavg(k) = tgasdata(iread,i)
        end do
!       2. avgker - averaging kernel (times pwf for column obs)
        do k = 1,navg
           do j = 1,nchanl
              iread = iread + 1
              avgker(j,k) = tgasdata(iread,i)
           end do
        end do
!       3. priorpro - a priori profile of observed variable
        do k = 1,navg
           iread = iread + 1
           priorpro(k) = tgasdata(iread,i)
        end do
!       4. priorobs - a priori value of observed variable
        do j = 1,nchanl
           iread = iread + 1
           priorobs(j) = tgasdata(iread,i)
        end do
     end if


!$$$ -- 3. COMPUTE GUESS VALUES --
!$$$       3A. AVERAGING KERNEL RETRIEVALS: ACOS GOSAT/OCO, MOPITT, IASI, ...
!===============================================================================
     if (useak) then
!       Map averaging kernel edge pressures/altitudes (peavg) to grid relative
!       vertical coordinate (grdpe)
        grdpe = peavg
        fsign = sign(1, int(peuse(2) - peuse(1)))
        do k = 1,nedge
           grdpe(k) = min(grdpe(k), maxval(peuse))
           grdpe(k) = max(grdpe(k), minval(peuse))
           call grdcrd1(grdpe(k), peuse, nsig+1, fsign)
        end do

!       Fill itgas (for int/stp routines) and gesobs (interpolate guess values
!       to obs time, location, and pressure, then average over obs layers)
        vname = trim(vname_tgas(ipos(1)))
        call gsi_chemguess_get('var::'//vname, itg, ier)
        itgas(:) = itg

        if (ier /= 0) then
           write(0,*) myname // ': ', vname, ' not found in chem bundle'
           call die(myname)
        end if

        ges = ges_tgas(:,:,:,:,itg)

!       Transform when input units are log-based
        vunit = trim(vunit_tgas(ipos(1)))
        if (vunit == 'ln' ) ges =   log(ges)
        if (vunit == 'log') ges =   log(ges)
        if (vunit == 'l10') ges = log10(ges)

!       Interpolate in time & horizontal space
        call tintrp2a1(ges, geslmod, grdlat, grdlon, grdtime, hrdifsig,        &
                       nsig, mype, nfldsig)

!       Average onto averaging kernel levels:
        ! for reactive trace gases, recalculate AMF using scattering weights and background trace gas profiles
        if ( isdoas ) then

           ! Get guess fields on obs levels (averaging kernel levels) 
           call zavgtgas_(pemod, geslmod, grdpe, gespro, avgwgt)

           ! get background (= a-priori) NO2/SO2 on model levels
           priortgas(:) = 0.0
           if ( obstype=='omso2' .or. obstype=='nmso2' ) then
              call tintrp2a1(so2_priorc, priortgas, grdlat, grdlon, grdtime, hrdifsig, &
                             nsig, mype, nfldsig)
           else
              call tintrp2a1(no2_priorc, priortgas, grdlat, grdlon, grdtime, hrdifsig, &
                             nsig, mype, nfldsig)
           endif
     
           ! map q onto obs levels
           call zavgtgas_(pemod, qvmod, grdpe, qvavg, tmpwgt)

           ! map prior tgas onto obs levels
           priorpro(:) = 0.0
           call zavgtgas_(pemod, priortgas, grdpe, priorpro, tmpwgt)

           ! convert both the a-priori and the current guess profile from mol/mol to 1e15 molec cm-2
           ! also accumulate guess observation, which is sum of all partial columns of the a-priori profile
           do k = 1,navg
              surfscal = max(0.0,min(1.0,(pemod(1)-peavg(k+1))/(peavg(k)-peavg(k+1))))
              ! conversion factor to go from v/v dry to molec cm-3. 
              ! Multiply by 1-Q for conversion of kg/kg dry to kg/kg total.
              lscal = surfscal * (peavg(k)-peavg(k+1)) * vv_to_molec * (1.-qvavg(k))
              gespro(k) = gespro(k) * lscal
              priorpro(k) = priorpro(k) * lscal
              avgwgt(k,:) = avgwgt(k,:) * lscal
           enddo

           ! for amftrop calculation: need to know pbl level at this location:          
           call tintrp2a11(pbltopl, pbl, grdlat, grdlon, grdtime, hrdifsig, mype, nfldsig)
           ipbl = nint(pbl)

           ! calculate AMF and actual averaging kernel using scattering weights and a-priori columns 
           amf     = sum(avgker(1,:)*gespro(:)) / sum(gespro)
           amftrop = sum(avgker(1,1:ipbl)*gespro(1:ipbl)) / sum(gespro(1:ipbl))
           avgker(1,:) = avgker(1,:) / amf
           ! cap averaging kernels at 1.0
           ! cakelle2, 20240116: now allow Ak's > 1.0 for doas type assimilation
           !where(avgker(1,:)>1.0) avgker(1,:) = 1.0

           ! a-priori observation
           priorobs(1) = sum(priorpro)

!       a. Mixing ratios in mol/mol, log, ...
        elseif (vunit /= 'molm2') then
           call zavgtgas_(pemod, geslmod, grdpe, gespro, avgwgt)

!       b. Partial columns in mol/m2 (fix hard-coded constants below)
        else
!          Undry profile because pressures are total-air
           geslmod = geslmod * (1.0 - qvmod)
           call zavgtgas_(pemod, geslmod, grdpe, gespro, avgwgt)

           do k = 1,navg
              avgwgt(k,:) = avgwgt(k,:) * (1.0 - qvmod)
           end do

!          Compute total-air subcolumns (dpges) and convert to mol/m2
           do k = 1,nedge
!             Copying what's in tintrp31
              iz  = grdpe(k)
              iz  = max(1, min(iz, nsig+1))
              izp = min(iz+1, nsig+1)

              delz = grdpe(k) - float(iz)
              peges(k) = pemod(iz)*(1.0 - delz) + pemod(izp)*delz
           end do
           dpges = abs(peges(2:nedge) - peges(1:nedge-1))
!          Should be MAPL_GRAV and MAPL_MOLMW or something close
           dpges = 1.0E+5/(9.80665*28.965) * dpges

!          For now assuming navg = nedge - 1
           gespro = gespro * dpges

           do j = 1,nsig
              avgwgt(:,j) = avgwgt(:,j) * dpges
           end do

           if (ldebugob) then
              print *, 'dpmod      = '
              print *, dpmod
              print *, 'qvmod      = '
              print *, qvmod
              print *, 'geslmod    = '
              print *, geslmod
              print *, 'peges      = '
              print *, peges
              print *, 'dpges      = '
              print *, dpges
           end if
        end if

!       Convert model units to prior units (obs treated below)
!       (This is applied inverted instead of to the model.  Otherwise, would
!       have to account for the scaling in gradients, etc., which is harder
!       to maintain, IMHO)
        priorobs = priorobs / sclmod(1)
        priorpro = priorpro / sclmod(1)

!       Create guess corresponding to obs
        if ( isdoas ) then
           gesobs = matmul(avgker, gespro)
        else
           gesobs = priorobs + matmul(avgker, gespro - priorpro)
        endif

!$$$       B. SAMPLE MEASUREMENTS: EZ, OBSPACK, MLS, ACE-FTS, ...
!===============================================================================
     else
!       Map obs pressure (pchanl) to grid relative vertical coordinate (grdpch)
        grdpch = pchanl
        fsign  = sign(1, int(pluse(2) - pluse(1)))
        do k = 1,nchanl
           grdpch(k) = min(grdpch(k), maxval(pluse))
           grdpch(k) = max(grdpch(k), minval(pluse))
           call grdcrd1(grdpch(k), pluse, nsig, fsign)
        end do

!       Fill itgas (for int/stp routines) and gesobs (interpolate guess values
!       to obs time, location, and pressure)
        do k = 1,nchanl
           vname = trim(vname_tgas(ipos(k)))
           call gsi_chemguess_get('var::'//vname, itg, ier)
           itgas(k) = itg

           if (ier /= 0) then
              write(0,*) myname // ': ', vname, ' not found in chem bundle'
              call die(myname)
           end if

           ges = ges_tgas(:,:,:,:,itg)

!          Transform log-based obs
           vunit = trim(vunit_tgas(ipos(k)))
           if (vunit == 'ln' ) ges =   log(ges)
           if (vunit == 'log') ges =   log(ges)
           if (vunit == 'l10') ges = log10(ges)

!          Interpolate to time and space of obs
           call tintrp31(ges, gesobs(k), grdlat, grdlon, grdpch(k), grdtime,   &
                         hrdifsig, mype, nfldsig)
        end do

!       Hack to make int/stp routines think this has an averaging kernel
        avgker = zero
        avgwgt = zero
        do k = 1,nchanl
           avgker(k,k) = one

!          Copying what's in tintrp31
           iz = grdpch(k)
           iz = max(1, min(iz, nsig))  
           delz = grdpch(k) - float(iz)
           delz = max(zero, min(delz, one))

           izp = min(iz+1, nsig)
           avgwgt(k,iz)  = one - delz
           avgwgt(k,izp) = avgwgt(k,izp) + delz
        end do
        gespro = gesobs
     end if

!    Account for log-space of background in averaging kernel
     do k = 1,nchanl
        vunit = trim(vunit_tgas(ipos(k)))
        if (vunit == 'ln' .or. vunit == 'log') then
           avgker(k,:) = avgker(k,:) / exp(gespro)
        else if (vunit == 'l10') then
           avgker(k,:) = avgker(k,:) / (log(10._r_kind) * 10._r_kind**gespro)
        end if
     end do

!    Print diagnostic output if debugging
     if (ldebugob) then
        print *,    'pchanl     = '
        print *,    pchanl
        print *,    'isbad      = '
        print *,    isbad
        print *,    'uncert     = '
        print *,    uncert

        if (useak) then
           print *, 'peavg      = '
           print *, peavg
           print *, 'grdpe      = '
           print *, grdpe
           print *, 'avgker     = '
           print *, avgker
           print *, 'priorpro   = '
           print *, priorpro
           print *, 'gespro     = '
           print *, gespro
           print *, 'priorobs   = '
           print *, priorobs
        end if
     end if


!$$$ -- 4.  CALCULATE INNOVATIONS, PERFORM GROSS CHECKS, AND ACCUMUALTE STATS
!===============================================================================
     do k = 1,nchanl
        j = ipos(k)

!       Get observation
        obs  = tgasdata(nreal+k,i)

        ! local copy of uncertainty. Probably not needed
        iuncert = uncert(k)

        ! for DOAS style obs, convert SCD to VCD using AMF calculated above
        if ( isdoas ) then
           scd_tot   = obs
           obs = scd_tot / amf
           ! get obs uncertainty 
           iuncert = obs * ( iuncert / scd_tot )

           if ( trim(obstype) == 'omso2' .or. trim(obstype)=='nmso2' ) then
              minsigobs  = tgas_minsigobs(2)
              maxsigobs  = tgas_maxsigobs(2)
              sigobsscal = tgas_sigobsscal(2)
           else
              minsigobs  = tgas_minsigobs(1)
              maxsigobs  = tgas_maxsigobs(1)
              sigobsscal = tgas_sigobsscal(1)
           endif
           iuncert = min(max(iuncert,minsigobs),maxsigobs)*sigobsscal
           !if ( sza(k) > 50.0 ) uncert(k) = uncert(k) * 1.+max(0.0,((sza(k)-50.0)/(tgas_szamax-sza(k))))
        endif

!       Subtract any specified bias from obs
        bias = bias_tgas(j)
        obs  = obs - bias

!       Convert model units to prior units (applied as inverse)
        obs = obs / sclmod(k)
        iuncert = iuncert / sclmod(k)

!       Compute observation innovation and error standard deviation
        obsinv(k) = obs - gesobs(k)
        error(k)  = sclstd(k) * iuncert 

!       Toss the obs that fail qc and gross checks
!       and set inverse obs error variance (varinv) and
!       square of ratio of final obs error to original obs error (raterr2)
        if (isbad(k) == 0 .and. abs(obsinv(k)/error(k)) <= gross(k)) then
           qcuse(k)   = 1
           varinv(k)  = one/(error(k)**2)
           raterr2(k) = one
        else
           qcuse(k)   = 0
           varinv(k)  = zero
           raterr2(k) = zero
           if (luse(i)) stats_tgas(2,j) = stats_tgas(2,j) + one                 ! # obs tossed
        end if

        if (ldebugob) then
           print *, 'gesobs     = ', gesobs(k),  'obs        = ', obs
           print *, 'o-g        = ', obsinv(k)
           print *, 'raterr2    = ', raterr2(k), 'varinv     = ', varinv(k)
           print *, '---'
        end if

!       Accumulate numbers for statistics
        if (qcuse(k) == 1 .and. luse(i)) then
           if (-1 < iouse(k) .or. (-1 == iouse(k) .and. ltgasdiagsave)) then
              omg = obsinv(k)

              stats_tgas(1,j) = stats_tgas(1,j) + one                           ! # obs used
              stats_tgas(3,j) = stats_tgas(3,j) + omg                           ! (o-g)
              stats_tgas(4,j) = stats_tgas(4,j) + omg*omg                       ! (o-g)**2
              stats_tgas(5,j) = stats_tgas(5,j) + omg*omg*varinv(k)*raterr2(k)  ! penalty
              stats_tgas(6,j) = stats_tgas(6,j) + obs                           ! obs

!             Compute hemispheric totals
              if (ltgasdiagsave) then
                 if (23 < deglat) then
                    stats_tgas(10,j) = stats_tgas(10,j) + one
                    stats_tgas(11,j) = stats_tgas(11,j) + omg
                    stats_tgas(12,j) = stats_tgas(12,j) + omg*omg*varinv(k)*raterr2(k)
                 else if (deglat < -23) then
                    stats_tgas(13,j) = stats_tgas(13,j) + one
                    stats_tgas(14,j) = stats_tgas(14,j) + omg
                    stats_tgas(15,j) = stats_tgas(15,j) + omg*omg*varinv(k)*raterr2(k)
                 end if
              end if

!             I think this needs to come out of the luse control and
!             way up to the top if raterr2 is going to depend on it
              exp_arg  = -half*varinv(k)*omg**2
              errorinv = sqrt(varinv(k))
              if (      tiny_r_kind < errorinv .and. tiny_r_kind < pg_tgas(j)  &
                  .and. tiny_r_kind < one - pg_tgas(j)) then
                 cg_tgas = b_tgas(j)*errorinv
                 wgross  = cg_term * pg_tgas(j)/((one - pg_tgas(j))*cg_tgas)
                 arg     = exp(exp_arg)
                 term    = log(arg + wgross) - log(one + wgross)
                 wgt     = one - wgross/(arg + wgross)
              else
                 term    = exp_arg
                 wgt     = one
              end if
!             raterr2(k) = term/exp_arg

              stats_tgas(8,j) = stats_tgas(8,j) - two*raterr2(k)*term           ! qc penalty
              if (wgt < wgtlim) stats_tgas(9,j) = stats_tgas(9,j) + one         ! hueristic obs count
           end if

           stats_tgas(7,j) = stats_tgas(7,j) + one                              ! # obs kept
        end if

!       Optionally save data for diagnostics
        if (ltgasdiagsave .and. luse(i)) then
           rdiagbuf(1,k,idout) = obs                         ! obs
           rdiagbuf(2,k,idout) = obsinv(k)                   ! obs - ges
           rdiagbuf(3,k,idout) = sqrt(varinv(k)*raterr2(k))  ! inverse obs error std dev
           rdiagbuf(4,k,idout) = pchanl(k)                   ! (nominal) obs pressure (hPa)
           rdiagbuf(5,k,idout) = pemod(1)                    ! model surface pressure (hPa)
           rdiagbuf(6,k,idout) = soundn                      ! sounding number

           ! undefined: nlevs
           if (netcdf_diag) then
              call nc_diag_metadata("TopLevelPressure",sngl(pemod(nsig+1)) )
              call nc_diag_metadata("BottomLevelPressure",sngl(pemod(1)) )
              call nc_diag_metadata("MPI_Task_Number", mype                      )
              call nc_diag_metadata("Latitude",        sngl(tgasdata(idlat,i))       )
              call nc_diag_metadata("Longitude",       sngl(tgasdata(idlon,i))       )
              call nc_diag_metadata("Time",            sngl(tgasdata(itime,i)-time_offset) )
              call nc_diag_metadata("Analysis_Use_Flag",      iouse(k)           )
              call nc_diag_metadata("Observation",     sngl(obs))
              call nc_diag_metadata("Inverse_Observation_Error",    sngl(errorinv))
              call nc_diag_metadata("Input_Observation_Error", sngl(error(k)))
              call nc_diag_metadata("Obs_Minus_Forecast_adjusted",  sngl(obsinv(k)))
              call nc_diag_metadata("Obs_Minus_Forecast_unadjusted",sngl(obsinv(k)))
              call nc_diag_metadata("Forecast_unadjusted", sngl(gesobs(k)))
              call nc_diag_metadata("Forecast_adjusted",sngl(gesobs(k)))
              if ( isdoas ) then
                 call nc_diag_metadata("Slant_Column_Density",      sngl(scd_tot) )
                 call nc_diag_metadata("SCD_error",                 sngl(uncert(k)) )
                 call nc_diag_metadata("Air_Mass_Factor"     ,      sngl(amf) )
                 call nc_diag_metadata("Air_Mass_Factor_Trop",      sngl(amftrop) )
                 call nc_diag_metadata("PBL_level_index"     ,      sngl(ipbl))
                 call nc_diag_metadata("A_Priori_VCD"        ,      sngl(priorobs(1)) )
                 call nc_diag_metadata("Solar_Zenith_Angle"  ,      sngl(rmiss) )
                 call nc_diag_metadata("Scan_Position"       ,      sngl(rmiss) )
              else
                 call nc_diag_metadata("Slant_Column_Density",      sngl(rmiss) )
                 call nc_diag_metadata("SCD_error",                 sngl(rmiss) )
                 call nc_diag_metadata("Air_Mass_Factor"     ,      sngl(rmiss) )
                 call nc_diag_metadata("Air_Mass_Factor_Trop",      sngl(rmiss) )
                 call nc_diag_metadata("PBL_level_index"     ,      sngl(rmiss) )
                 call nc_diag_metadata("A_Priori_VCD"        ,      sngl(rmiss) )
                 call nc_diag_metadata("Solar_Zenith_Angle",        sngl(rmiss) )
                 call nc_diag_metadata("Scan_Position",             sngl(rmiss) )
              endif 
              if (wrtgeovals) then
                 call nc_diag_data2d("mole_fraction_of_species", sngl(geslmod))
                 call nc_diag_data2d("air_pressure_levels",sngl(plmod))
                 call nc_diag_data2d("a_priori_profile",sngl(priorpro))
                 call nc_diag_data2d("guess_profile",sngl(gespro))
                 call nc_diag_data2d("averaging_kernel",sngl(avgker(1,:)))
              endif
           endif
        end if

!       If not assimilating this observation, reset qcuse and variances
        if (iouse(k) < 1) then
           qcuse(k)   = 0
           varinv(k)  = zero
           raterr2(k) = zero
        end if

!       Hack to pre-condition obs and avoid some GSI badness
        if (qcuse(k) == 1) then
           obsinv(k) = obsinv(k) / error(k)
           varinv(k) = one
           do j = 1,navg
              avgker(k,j) = avgker(k,j) / error(k)
           end do
        end if
     end do

!    Check all information for obs.  If there is at least one piece of
!    information that passed quality control, use this observation.
     ikeep = maxval(qcuse)
end if ! (in_curbin)

     if (ldebug) then
        print *, 'lmaybepassive = ', lmaybepassive
        print *, 'ikeep = ', ikeep
     end if

!    In principle, we want ALL obs in the diagnostics structure, but for
!    passive obs (monitoring), it is difficult to do if ltgasdiagsave
!    is not on in the first outer loop. For now we use lmaybepassive
     if (lmaybepassive) then
!       Link observation to appropriate observation bin
        if (nobs_bins > 1) then
           ibin = nint(grdtime*60/mn_obsbin) + 1
        else
           ibin = 1
        end if
        if (ibin < 1 .or. nobs_bins < ibin) then
           write(0,*) myname // ': ', mype, ' error nobs_bins, ibin = ',       &
                      nobs_bins, ibin
           call die(myname)
        end if

        if (ldebug) print *, 'ibin  = ', ibin, 'nobs_bins = ', nobs_bins

        if (luse_obsdiag) my_diagLL => odiagLL(ibin)

if (in_curbin) then
!       Process obs that have at least one piece of information that passed qc
!       checks
        if (.not. last .and. ikeep == 1) then
           allocate(my_head)
           call tgasNode_appendto(my_head,tgashead(ibin))

           my_head%idv  = is
           my_head%iob  = ioid(i)
           my_head%elat = grdlat
           my_head%elon = grdlon

           my_head%luse    = luse(i)
           my_head%time    = grdtime
           my_head%nchanl  = nchanl
           my_head%navg    = navg
           my_head%obstype = obstype

           allocate(my_head%res(nchanl),       my_head%err2(nchanl),           &
                    my_head%raterr2(nchanl),   my_head%ipos(nchanl),           &
                    my_head%itgas(nchanl),     my_head%avgker(nchanl,navg),    &
                    my_head%avgwgt(navg,nsig), stat=istatus)
           if (istatus /= 0) then
              write(0,*) myname // ': allocation error for tgas pointer, ' //  &
                         'istatus = ', istatus
              call die(myname)
           end if
           if (luse_obsdiag) then
              allocate(my_head%diags(nchanl), stat=istatus)
              if (istatus /= 0) then
                 write(0,*) myname // ': allocation error for tgas ' //        &
                            'pointer, istatus = ', istatus
                 call die(myname)
              end if
           end if

!          Set (i,j) indices of guess gridpoint that bound obs location
           call get_ij(mm1, grdlat, grdlon, my_head%ij, my_head%wij)

           if (ldebug) then
              print *, 'ij  = ', my_head%ij
              print *, 'wij = ', my_head%wij
           end if

           my_head%res     = obsinv
           my_head%err2    = varinv
           my_head%raterr2 = raterr2
           my_head%ipos    = ipos
           my_head%itgas   = itgas
           my_head%avgker  = avgker
           my_head%avgwgt  = avgwgt

           my_head => null()
        end if ! (.not. last .and. ikeep == 1)
end if ! (in_curbin)

!       Link obs to diagnostics structure
        do k = 1,nchanl
           if (luse_obsdiag) then
              my_diag => obsdiagLList_nextNode(my_diagLL,                      &
                 create = .not. lobsdiag_allocated, idv = is, iob = ioid(i),   &
                 ich = k, elat = grdlat, elon = grdlon,                        &
                 luse = luse(i) .and. (qcuse(k) == 1), miter = miter)
           end if

           if (.not. associated(my_diag)) then
              write(0,*) myname // ': obsdiagLList_nextNode(), create = ',     &
                         .not. lobsdiag_allocated
              call die(myname)
           end if

if (in_curbin) then
           if (luse_obsdiag) then
              call obsdiagNode_set(my_diag, wgtjo = varinv(k)*raterr2(k),      &
                 jiter = jiter, muse = (ikeep==1), nldepart = obsinv(k))
           end if

           if (.not. last .and. ikeep == 1) then
              my_head => tailNode_typecast_(tgashead(ibin))
              if (.not. associated(my_head)) then
                 write(0,*) myname // ': unexpected, associated(my_head) = ',  &
                            associated(my_head)
                 call die(myname)
              end if

              if (luse_obsdiag) then
                 call obsdiagNode_assert(my_diag, my_head%idv, my_head%iob,    &
                    k, myname, 'my_diag:my_head')
                 my_head%diags(k)%ptr => my_diag
              end if
              my_head => null()
           end if

           if (ltgasdiagsave .and. lobsdiagsave .and. luse(i)) then
              associate(odiag => my_diag)
                 jdiag = irmin
                 do jj = 1,miter
                    jdiag = jdiag + 1
                    if (odiag%muse(jj)) then
                       rdiagbuf(jdiag,k,idout) =  one
                    else
                       rdiagbuf(jdiag,k,idout) = -one
                    end if
                 end do
                 do jj = 1,miter+1
                    jdiag = jdiag + 1
                    rdiagbuf(jdiag,k,idout) = odiag%nldepart(jj)
                 end do
                 do jj = 1,miter
                    jdiag = jdiag + 1
                    rdiagbuf(jdiag,k,idout) = odiag%tldepart(jj)
                 end do
                 do jj = 1,miter
                    jdiag = jdiag + 1
                    rdiagbuf(jdiag,k,idout) = odiag%obssen(jj)
                 end do
              end associate ! odiag
           end if
end if ! (in_curbin)
        end do ! (k = 1,nchanl)

     else ! (.not. lmaybepassive)
if (in_curbin) then
        if (ltgasdiagsave .and. lobsdiagsave .and. luse(i)) then
           rdiagbuf(irmin+1:irdim1,1:nchanl,idout) = zero
        end if
end if ! (in_curbin)
 
     end if ! (lmaybepassive)
  end do ! (i = 1,nobs)

! If requested, write to diagnostic file
  if (ltgasdiagsave) then

     if (netcdf_diag) call nc_diag_write

     if (binary_diag .and. idout>0 ) then
        write(string,100) jiter
100     format('_',i2.2)
        diag_tgas_file =  trim(dirname) // trim(isis) // string
        if (init_pass) then
           open(4, file=diag_tgas_file, form='unformatted', status='unknown',     &
                position='rewind')
           if (mype == mype_diaghdr(is)) then
              write(6,*) trim(myname), ': write header record for ', isis, iint,  &
                         ireal, idout, ' to file ', trim(diag_tgas_file), ' ',    &
                         ianldate
              write(4) isis, dplat(is), obstype, jiter, nchanl, ianldate, iint,   &
                       ireal, irdim1, irmin
              write(4) real(pchanl,r_single), real(gross,r_single),               &
                       real(sclstd,r_single), iouse
           end if
        else
           open(4, file=diag_tgas_file, form='unformatted', status='old',         &
                position='append')
        end if
        write(4) idout
        write(4) idiagbuf(:,1:idout), diagbuf(:,1:idout), rdiagbuf(:,:,1:idout)
        close(4)
     endif
  end if

! Clean up
  deallocate(itgas, gespro, avgker, avgwgt)
  if (useak) then
     deallocate(peavg, grdpe, priorpro)
     deallocate(peges, dpges)
  else
     deallocate(grdpch)
  end if
  ! doas stuff
  if(allocated(tmpwgt    )) deallocate(tmpwgt)
  if(allocated(qvavg     )) deallocate(qvavg)
  if(allocated(delps     )) deallocate(delps)

  if (ltgasdiagsave) deallocate(rdiagbuf)

! End of routine
  call final_vars_
  return

  contains

  function tailNode_typecast_(oll) result(ptr_)
! Cast the tailNode of oll to an obsNode, as in
!      ptr_ => typecast_(tailNode_(oll))
     use m_tgasNode, only: tgasNode, typecast_ => tgasNode_typecast
     use m_obsLList, only: obsLList, tailNode_ => obsLList_tailNode
     use m_obsNode,  only: obsNode

     implicit none

     type(tgasNode), pointer :: ptr_
     type(obsLList), target, intent(in) :: oll
     class(obsNode), pointer :: inode_

     inode_ => tailNode_(oll)
     ptr_   => typecast_(inode_)
  end function tailNode_typecast_

  subroutine init_vars_
     use guess_grids,      only: geop_hgti
     use gsi_metguess_mod, only: gsi_metguess_bundle

     character(len=*), parameter :: myname_ = myname // '::init_vars_'

     real(r_kind), dimension(:,:  ), pointer :: rank2 => NULL()
     real(r_kind), dimension(:,:,:), pointer :: rank3 => NULL()

     character(len=varlen) :: vname

     integer(i_kind) :: i, k, ier

     if (size(gsi_chemguess_bundle) /= nfldsig) then
        write(0,*) myname_ // ': inconsistent vector sizes (nfldsig, ' //      &
                   'size(gsi_chemguess_bundle) = ', nfldsig,                   &
                   size(gsi_chemguess_bundle)
        call die(myname_)
     end if

!    Get dimensions and allocate tgas variables
     if (0 < ntgas) then
        vname = trim(tgnames(1))
        call gsi_bundlegetpointer(gsi_chemguess_bundle(1), vname, rank3, ier)

        if (ldebug) print *, 'size(ges_tgas) = ', size(rank3,1),               &
                             size(rank3,2), size(rank3,3), nfldsig, ntgas

        allocate(     ges(size(rank3,1),size(rank3,2),size(rank3,3),nfldsig))
        allocate(ges_tgas(size(rank3,1),size(rank3,2),size(rank3,3),nfldsig,   &
                          ntgas))
     end if

!    If required guess vars available, extract from bundle
     do n = 1,ntgas
        vname = trim(tgnames(n))

        do i = 1,nfldsig
           call gsi_bundlegetpointer(gsi_chemguess_bundle(i), vname, rank3, ier)
           ges_tgas(:,:,:,i,n) = rank3

           if (ier /= 0) then
              write(0,*) myname_ // ': ', vname, ' not found in chem bundle'
              call die(myname_)
           end if
        end do
     end do

!    You should recreate the correct ges_ze here from z and tv.  At
!    present, it could change if a switch is changed in guess_grids, and it
!    doesn't match the definition in the GCM (fixme)
     if (allocated(ges_ze)) then
        write(0,*) myname_ // ': ges_ze already incorrectly allocated'
        call die(myname_)
     end if

     allocate(ges_ze(size(rank3,1),size(rank3,2),size(rank3,3)+1,nfldsig))

     do i = 1,nfldsig
        call gsi_bundlegetpointer(gsi_metguess_bundle(i), 'z', rank2, ier)

        if (ier /= 0) then
           write(0,*) myname // ': z not found in met bundle'
           call die(myname_)
        end if

        do k = 1,nsig+1
           ges_ze(:,:,k,i) = rank2 + geop_hgti(:,:,k,i)
        end do
     end do

!    Allocate and extract water vapor from bundle
     if (allocated(ges_qv)) then
        write(0,*) myname_ // ': ges_qv already incorrectly allocated'
        call die(myname_)
     end if

     allocate(ges_qv(size(rank3,1),size(rank3,2),size(rank3,3),nfldsig))

     do i = 1,nfldsig
        call gsi_bundlegetpointer(gsi_metguess_bundle(i), 'q', rank3, ier)
        ges_qv(:,:,:,i) = rank3

        if (ier /= 0) then
           write(0,*) myname // ': q not found in met bundle'
           call die(myname_)
        end if
     end do
  end subroutine init_vars_

  subroutine final_vars_
     if (allocated(ges))      deallocate(ges)
     if (allocated(ges_tgas)) deallocate(ges_tgas)
     if (allocated(ges_ze))   deallocate(ges_ze)
     if (allocated(ges_qv))   deallocate(ges_qv)
  end subroutine final_vars_

  subroutine init_levs_(useak, nreal, nchanl, navg, nedge)
     character(len=*), parameter    :: myname_ = myname // '::init_levs_'

     logical,         intent(in   ) :: useak
     integer(i_kind), intent(in   ) :: nreal, nchanl
     integer(i_kind), intent(  out) :: navg,  nedge

     integer(i_kind) :: nn

!    If not using averaging kernels, just return nchanl for navg and a value
!    for nedge that will throw an error
     if (.not. useak) then
        navg  = nchanl
        nedge = 0
        return
     end if

!    We can determine navg and nedge from nreal and nchanl as follows
!    ---
!    nreal = 6 + 3*nchanl + nedge + (navg + navg*nchanl + nchanl)
!    Let nn = nreal - 6 - 4*nchanl.  Then
!        nn = nedge + navg + navg*nchanl
!
!    There are two options:
!    1) nedge = navg
!       nn    = navg*(2 + nchanl)
!       navg  = nn/(2 + nchanl)
!    2) nedge = navg + 1,
!       nn    = navg*(2 + nchanl) + 1
!       navg  = (nn - 1)/(2 + nchanl)
!
!    However, only one can be true because if both nn/(2 + nchanl) and
!    (nn - 1)/(2 + nchanl) are integers, then their difference, 1/(2 + nchanl)
!    is also an integer, which is impossible since nchanl > 0.

!    First assume option 1, if untrue fill with option 2
     nn    = nreal - 6 - 4*nchanl
     navg  = nn/(2 + nchanl)
     nedge = navg
     if (navg*(2 + nchanl) /= nn) then
        navg  = (nn - 1)/(2 + nchanl)
        nedge = navg + 1
     end if

!    Make sure we got it right
     if (nn /= nedge + navg + navg*nchanl) then
        write(0,*) myname_ // ': *** ERROR: cannot determine navg and nedge ***'
        write(0,*) myname_ // ':     nreal = ', nreal, ', nchanl = ', nchanl
        call die(myname_)
     end if
  end subroutine init_levs_

  subroutine zavgtgas_(pe, f, ke, g, w)
     logical,          parameter :: ldebug_ = .false.
     character(len=*), parameter :: myname_ = myname // '::zavgtgas_'

     real(r_kind), intent(in   ) :: pe(nsig+1), f(nsig), ke(nedge)
     real(r_kind), intent(  out) :: g(navg), w(navg,nsig)

     real(r_kind)    :: avg(nedge-1), wgt(nedge-1,nsig)

     real(r_kind)    :: dz1, dz2, delp, delz, dlay
     integer(i_kind) :: nlays, iz1, iz2, j, l
     logical         :: top2bot

     nlays = nedge - 1
     avg   = zero
     wgt   = zero

!    Determine if levels are ordered from top to bottom (top2bot)
     top2bot = ke(1) < ke(nedge)

     do j = 1,nlays
        if (top2bot) then
          dz1 = ke(j+1)
          dz2 = ke(j)
        else
          dz1 = ke(j)
          dz2 = ke(j+1)
        end if

        iz1 = min(int(dz1), nsig)
        iz2 = min(int(dz2), nsig)

!       Obs layer falls inside single model layer
        if (iz1 == iz2) then
           avg(j)     = f(iz1)
           wgt(j,iz1) = one
           cycle
        end if

!       Obs layer covers multiple model layers
        dlay = zero
        do l = iz1,iz2,-1
!          Weight is 1 unless a boundary value
           delz = one
           if (l == iz1) delz =         dz1 - iz1
           if (l == iz2) delz = delz - (dz2 - iz2)

           delp = pe(l) - pe(l+1)

           dlay     = dlay   +      delp*delz
           avg(j)   = avg(j) + f(l)*delp*delz
           wgt(j,l) =               delp*delz
        end do

        avg(j)   = avg(j)   / dlay
        wgt(j,:) = wgt(j,:) / dlay

        if (ldebug_) then
           print *, 'iz1  = ', iz1,  'iz2  = ', iz2
           print *, 'dlay = ', dlay, 'avg  = ', avg(j)
        end if
     end do

!    Is it a layer profile or edge profile?
     if (nedge == navg+1) then
!       Copy layer data to outputs
        g = avg
        w = wgt
     else if (nedge == navg) then
!       Interpolate layer data to edge data
        g(1) = avg(1)
        do j = 2,nlays
           g(j) = 0.5*(avg(j-1) + avg(j))
        end do
        g(nedge) = avg(nlays)

        do l = 1,nsig
           w(1,l) = wgt(1,l)
           do j = 2,nlays
              w(j,l) = 0.5*(wgt(j-1,l) + wgt(j,l))
           end do
           w(nedge,l) = wgt(nlays,l)
        end do
     else
!       Throw an error if nedge & navg don't match up
        write(0,*) myname_ // ': *** ERROR: Inconsistent averaging levels ' // &
                                     'and edge numbers ***'
        write(0,*) myname_ // ': *** nedge must be navg or navg+1'
        write(0,*) myname_ // ':     nedge = ', nedge, ', navg   = ', navg
        call die(myname_)
     end if
  end subroutine zavgtgas_
                         
  subroutine init_netcdf_diag_
  character(len=80) string
  integer(i_kind) ncd_fileid,ncd_nobs 
  logical append_diag
  logical,parameter::verbose=.true.
           
     write(string,900) jiter
900  format('_',i2.2,'.nc4')
     filex=obstype
     diag_tgas_file = trim(dirname) // trim(filex) // '_' // trim(dplat(is)) // (string)
        
     inquire(file=diag_tgas_file, exist=append_diag)
        
     if (append_diag) then
        call nc_diag_read_init(diag_tgas_file,ncd_fileid)
        ncd_nobs = nc_diag_read_get_dim(ncd_fileid,'nobs')
        call nc_diag_read_close(diag_tgas_file)
  
        if (ncd_nobs > 0) then
           if(verbose) print *,'file ' // trim(diag_tgas_file) // ' exists.  Appending.  nobs,mype=',ncd_nobs,mype
        else
           if(verbose) print *,'file ' // trim(diag_tgas_file) // ' exists but contains no obs.  Not appending. nobs,mype=',ncd_nobs,mype
           append_diag = .false. ! if there are no obs in existing file, then do not try to append
        endif
     end if
     
     call nc_diag_init(diag_tgas_file, append=append_diag)

     if (.not. append_diag) then ! don't write headers on append - the module will break?
        call nc_diag_header("date_time",ianldate )
        call nc_diag_header("Satellite_Sensor", isis)
        call nc_diag_header("Satellite", dplat(is))
        call nc_diag_header("Observation_type", obstype)
     endif

  end subroutine init_netcdf_diag_
  subroutine contents_binary_diag_
  end subroutine contents_binary_diag_
  subroutine contents_netcdf_diag_
! Observation class
  character(7),parameter     :: obsclass = '  tgas'
! contents interleafed above should be moved here (RTodling)
  end subroutine contents_netcdf_diag_

end subroutine setuptgas
end module tgas_setup
