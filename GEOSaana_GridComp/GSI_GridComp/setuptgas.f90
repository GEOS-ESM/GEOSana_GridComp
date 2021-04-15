subroutine setuptgas(lunin, mype, stats_tgas, nchanl, nreal, nobs, obstype,    &
                     isis, is, ltgasdiagsave, init_pass, last_pass)
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
!     last_pass     - state of setup processing
!     stats_tgas    - sums for various statistics as a function of level
!
!   output argument list:
!     stats_tgas    - sums for various statistics as a function of level
!
!$$$ end documentation block

! Declare module variables
  use mpeu_util,         only: die, perr
  use kinds,             only: r_kind, r_single, i_kind
  use constants,         only: zero, half, one, two, tiny_r_kind, huge_r_kind, &
                               cg_term, wgtlim, h300
  use obsmod,            only: tgashead, tgastail, i_tgas_ob_type,             &
                               dplat, nobskeep, mype_diaghdr, dirname,         &
                               time_offset, ianldate, obsdiags,                &
                               lobsdiag_allocated, lobsdiagsave,               &
                               tgas_ob_type, obs_diag, nlmopitt, nlacos,       &
                               nlflask
  use gsi_4dvar,         only: nobs_bins, hr_obsbin
  use gridmod,           only: get_ij, nsig
  use guess_grids,       only: nfldsig, ges_prsi, ntguessig, hrdifsig
  use gsi_bundlemod,     only: gsi_bundlegetpointer
  use gsi_chemguess_mod, only: gsi_chemguess_get, gsi_chemguess_bundle
  use tgasinfo,          only: jpch_tgas, pchanl_tgas, error_tgas, gross_tgas, &
                               nusis_tgas, iuse_tgas, b_tgas, pg_tgas,         &
                               ihave_co, ihave_co2
  use jfunc,             only: jiter, last, miter
  use m_dtime,           only: dtime_setup, dtime_check, dtime_show

  implicit none

! Declare passed variables
  integer(i_kind),  intent(in   ) :: lunin, mype, nchanl, nreal, nobs, is
  character(20),    intent(in   ) :: isis
  character(len=*), intent(in   ) :: obstype
  logical,          intent(in   ) :: ltgasdiagsave, init_pass, last_pass
  real(r_kind),     intent(inout) :: stats_tgas(9,jpch_tgas)

! Declare local parameters  
  integer(i_kind),  parameter :: iint  = 1
  integer(i_kind),  parameter :: ireal = 3
  integer(i_kind),  parameter :: irmin = 6
  real(r_kind),     parameter :: r10   = 10._r_kind
  real(r_kind),     parameter :: r19   = 19._r_kind
  real(r_kind),     parameter :: tol_error = 1.e4_r_kind
  character(len=*), parameter :: myname  = 'setuptgas'
  logical,          parameter :: ldebug  = .false.
  logical,          parameter :: lcenter = .false.

! Declare local variables  
  real(r_kind) :: obs, omg, rat_err2, grdtime, grdlat, grdlon, psurf
  real(r_kind) :: cg_tgas, wgross, wnotgross, wgt, arg, exp_arg, term
  real(r_kind) :: errorinv, fact

  real(r_kind), dimension(nchanl) :: obsinv, varinv3, ratio_errors
  real(r_kind), dimension(nchanl) :: pchanl, error, gross, tnoise
  real(r_kind), dimension(nchanl) :: priorobs, gesobs, uncert
  real(r_kind), dimension(nsig)   :: colmod, co2lmod
  real(r_kind), dimension(nsig+1) :: prsimod

  real(r_kind), dimension(nreal+nchanl,nobs) :: tgasdata

  real(r_single), dimension(ireal,nobs) :: diagbuf

  real(r_kind),   allocatable, dimension(:)       :: prsiobs, grdprsi
  real(r_kind),   allocatable, dimension(:)       :: priorpro, gespro
  real(r_kind),   allocatable, dimension(:,:)     :: avgker, avgwgt
  real(r_kind),   allocatable, dimension(:,:,:,:) :: ges_co, ges_co2
  real(r_single), allocatable, dimension(:,:,:)   :: rdiagbuf

  integer(i_kind) :: i, nchon, npro, nlays, jj, iextra, ibin, ifld
  integer(i_kind) :: k, l, kcut, j, jc, id, jdiag, irdim1, ier, istatus, ioff
  integer(i_kind) :: itoss, ikeep, nkeep, mm1
  integer(i_kind) :: isd, itime, iglat, iglon, idlat, idlon, iqcf, ipsurf

  integer(i_kind), dimension(iint,nobs) :: idiagbuf
  integer(i_kind), dimension(nchanl)    :: ipos, iouse
  integer(i_kind), dimension(nobs_bins) :: n_alloc
  integer(i_kind), dimension(nobs_bins) :: m_alloc

  character(12)  :: string
  character(10)  :: filex
  character(128) :: diag_tgas_file
  character(256) :: ctype

  logical :: l_may_be_passive, lproceed, in_curbin, in_anybin, lqcflag
  logical :: lacos, lmopitt, lflask

  logical, dimension(nobs) :: luse

  type(tgas_ob_type), pointer :: my_head
  type(obs_diag),     pointer :: my_diag

! Shortcut to not have to redo test many times
  lacos   = (trim(obstype) == 'acos')
  lmopitt = (trim(obstype) == 'mopitt')
  lflask  = (trim(obstype) == 'flask')

! Initialize observation dependent variables
  if (lacos) then
     npro  = nlacos
     nlays = nlacos - 1

     if (nchanl /= 1) then
        write(6,*) trim(myname), ': ***ERROR*** in channel number'
        write(6,*) trim(myname), ': nchanl = ', nchanl
        write(6,*) trim(myname), ' ***STOP***'
        call stop2(70)
     end if

  else if (lmopitt) then
     npro  = nlmopitt
     nlays = nlmopitt

     if (nchanl /= npro) then
        write(6,*) trim(myname), ': ***ERROR*** in channel number'
        write(6,*) trim(myname), ': nchanl = ', nchanl
        write(6,*) trim(myname), ' ***STOP***'
        call stop2(69)
     end if

  else if (lflask) then
     nlays = 0
     npro  = nlflask

!    if (nchanl /= npro) then
     if (nchanl /= 2) then
        write(6,*) trim(myname), ': ***ERROR*** in channel number'
        write(6,*) trim(myname), ': nchanl = ', nchanl
        write(6,*) trim(myname), ' ***STOP***'
        call stop2(69)
     end if

! Return if unsupported obstype
  else
     return

  end if

! Allocate local variables
  allocate(prsiobs(nlays+1), grdprsi(nlays+1))
  allocate(priorpro(npro), gespro(npro))
  allocate(avgker(nchanl,npro), avgwgt(npro,nsig))

! If required guess fields are unavaiable, simply return
  call check_vars_(lproceed)
  if (.not. lproceed) return

! If required guess fields are available, extract from bundle
  call init_vars_

  n_alloc(:) = 0
  m_alloc(:) = 0

  mm1 = mype + 1

! Initialize arrays
  do k = 1,nchanl
     ipos(k)   = 0
     pchanl(k) = zero
     gross(k)  = huge_r_kind
     tnoise(k) = huge_r_kind
     iouse(k)  = -2
  end do

! Locate data for satellite in tgasinfo arrays, and remove channels if
! necessary
  itoss = 1
  l_may_be_passive = .false.
  jc = 0
  if (ldebug) print *, 'jpch_tgas = ', jpch_tgas, 'nchanl = ', nchanl
  do j = 1,jpch_tgas
     if (ldebug) print *, 'isis = ', isis, 'nusis_tgas(j) = ', nusis_tgas(j)
     if (isis == nusis_tgas(j)) then
        jc = jc + 1
        if (nchanl < jc) then
           write(6,*) trim(myname), ': ***ERROR*** in channel number'
           write(6,*) trim(myname), ': jc    = ', jc, ', nchanl = ', nchanl
           write(6,*) trim(myname), ' ***STOP***'
           call stop2(71)
        end if

        ipos(jc)   = j
        pchanl(jc) = pchanl_tgas(j)
        gross(jc)  = min(r10*gross_tgas(j), h300)
        tnoise(jc) = error_tgas(j)
        iouse(jc)  = iuse_tgas(j)

        if (-1 < iouse(jc)) l_may_be_passive = .true.
        if (tnoise(jc) < tol_error) itoss = 0

        if (ldebug) then
           print *, 'pchanl = ', pchanl(jc)
           print *, 'gross  = ', gross(jc)
           print *, 'tnoise = ', tnoise(jc)
           print *, 'iouse  = ', iouse(jc)
           print *, 'ltgasdiagsave = ', ltgasdiagsave
        end if
     end if
  end do
  nchon = jc

! Handle error conditions
  if (nchon < nchanl) then
!    Unsure how to properly implement this, so reporting as an error
!    Instead, turn off channels with the iuse flag
!    if (mype == 0) write(6,*) trim(myname), ': level number reduced for ',    &
!                              obstype, ' ', nchanl,' --> ', nchon
     write(6,*) trim(myname), ': ***ERROR*** in channel number'
     write(6,*) trim(myname), ': nchon = ', nchon, ', nchanl = ', nchanl
     write(6,*) trim(myname), ' ***STOP***'
     call stop2(72)
  end if
  if (nchon == 0) then
     if (mype == 0) write(6,*) trim(myname), ': no channels found for ', isis
     if (nobs >  0) read(lunin) 
     goto 135
  end if
  if (itoss == 1) then
     if (mype == 0) write(6,*) trim(myname), ': all obs variances > 1.e4, ' // &
                               'do not use data from satellite ', isis
     if (nobs >  0) read(lunin)
     goto 135
  end if

! Initialize variables used in processing
  nkeep = 0
  do i = 1,nobs
     ikeep = 0
     do k = 1,nchanl
        if (0 < iouse(k) .or. ltgasdiagsave) ikeep = 1
     end do
     nkeep = nkeep + ikeep
  end do

! If none of the data will be assimilated and don't need diagnostics, return to
! calling program
  if (nkeep == 0) then
     call final_vars_
     return
  end if

! If requested, save data for diagnostic ouput
  if (ltgasdiagsave) then
     id     = 0
     irdim1 = irmin
     if (lobsdiagsave) irdim1 = irdim1 + 4*miter + 1
     allocate(rdiagbuf(irdim1,nchanl,nobs))
  end if

! Read and transform data
  read(lunin) tgasdata, luse

! Index information for data array (see read_tgas subroutine)
  isd    = 1   ! index of satellite (unused in ozone routines, fixme)
  itime  = 2   ! index of analysis relative obs time
  iglon  = 3   ! index of grid relative obs location (x)
  iglat  = 4   ! index of grid relative obs location (y)
  idlon  = 5   ! index of earth relative longitude (degrees)
  idlat  = 6   ! index of earth relative latitude (degrees)
  iqcf   = 7   ! index of obs quality control flag
  ipsurf = 8   ! index of surface pressure
      
  if (ldebug) then
     print *, 'ntguessig = ', ntguessig, 'nfldsig = ', nfldsig
     print *, 'hrdifsig  = '
     print *, hrdifsig
  end if

  call dtime_setup

  do i = 1,nobs
     grdtime = tgasdata(itime,i)
     grdlat  = tgasdata(iglat,i)
     grdlon  = tgasdata(iglon,i)
     lqcflag = (abs(tgasdata(iqcf,i) - one) < 1.e-6_r_kind)                    ! Paranoid casting to logical
     psurf   = tgasdata(ipsurf,i)

!    If debugging, allow centering of obs in space and time to eliminate
!    interpolation
     if (ldebug .and. lcenter) then
       if (mype == 0) write(6,*) trim(myname), ': *** WARNING: centering ',    &
                                 'obs in time and moving to nearest grid ',    &
                                 'point ***'
        grdtime = 3._r_kind
        grdlat  = real(int(grdlat))
        grdlon  = real(int(grdlon))
     end if

     call dtime_check(grdtime, in_curbin, in_anybin)
     if (.not. in_anybin) cycle

if (in_curbin) then
     if (nobskeep > 0) then
        write(6,*) trim(myname), ': nobskeep ', nobskeep
        call stop2(259)
     end if

!    Interpolate guess interface pressure (ges_prsi) to obs time and
!    location (prsimod)
     call tintrp2a1(ges_prsi, prsimod, grdlat, grdlon, grdtime, hrdifsig,      &
                    nsig+1, mype, nfldsig)

!    Read and construct obs variables:
!       1. priorobs = prior for observed variable
!       2. avgker   = averaging kernel (times pres. wght. fcn. for column obs)
!       3. gespro   = guess/model profile
!       4. priorpro = prior for   profile
!       5. uncert   = standard deviation of obs error
!    For observations of multiple species, the "profile" can contain
!    observations of each of those species; see flask obstype.  This
!    functionality is still being developed.
     if (lacos) then
        prsiobs(1) = psurf*1.e-4_r_kind
        do k = 1,nlays
           prsiobs(k+1) = psurf - real(nlays-k, r_kind)/r19*psurf
        end do
        pchanl = psurf

!       Map obs pressure (prsiobs) to grid relative vertical coordinate (grdprsi)
        grdprsi = prsiobs/r10
        do k = 1,nlays+1
           if (prsimod(1) < grdprsi(k)) grdprsi(k) = prsimod(1)
           call grdcrd1(grdprsi(k), prsimod, nsig+1, -1)
        end do

!       Interpolate guess values to obs time and location, then average over
!       obs layers
        call tintrp2a1(ges_co2, co2lmod, grdlat, grdlon, grdtime, hrdifsig,    &
                       nsig, mype, nfldsig)
        call zavgtgas_(co2lmod, .false., gespro, avgwgt)

!       Read a priori values, averaging kernel, and uncertainty
        do k = 1,npro
           priorpro(k) = tgasdata(k+8,i)
        end do
        do k = 1,npro
           avgker(1,k) = tgasdata(k+8+npro,i)
        end do

        priorobs(1) = tgasdata(1+8+2*npro,i)
        uncert(1)   = tgasdata(2+8+2*npro,i)

     else if (lmopitt) then
        prsiobs(1) = psurf
        do k = 2,nlays
           prsiobs(k) = 1000._r_kind - real(k-1, r_kind)*100._r_kind
           prsiobs(k) = min(prsiobs(k), psurf)
        end do
!       Not sure if this is the best choice (fixme)
        prsiobs(nlays+1) = 26._r_kind
        pchanl = prsiobs(1:nlays)

!       Map obs pressure (prsiobs) to grid relative vertical coordinate (grdprsi)
        grdprsi = prsiobs/r10
        do k = 1,nlays+1
           if (prsimod(1) < grdprsi(k)) grdprsi(k) = prsimod(1)
           call grdcrd1(grdprsi(k), prsimod, nsig+1, -1)
        end do

!       Interpolate guess values to obs time and location, then average over
!       obs layers
        call tintrp2a1(ges_co, colmod, grdlat, grdlon, grdtime, hrdifsig,      &
                       nsig, mype, nfldsig)
        call zavgtgas_(colmod, .true., gespro, avgwgt)

!       Read a priori values, averaging kernel, and uncertainty
        do k = 1,npro
           priorpro(k) = tgasdata(k+8,i)
        end do
        do k = 1,npro
           do j = 1,npro
              avgker(j,k) = tgasdata(k+(j-1)*npro+8+npro,i)
           end do
        end do

        do k = 1,npro
!          uncert(k) = sqrt(tgasdata(k+8+npro+npro**2,i))
!          For now, just use ones because product has errors (fixme)
           uncert(k) = one
        end do

!       Prevent nans when we take the log
        do kcut = 1,npro
           if (tiny_r_kind < priorpro(kcut)) exit
        end do
        priorpro(1:kcut-1) = log10(priorpro(kcut))
        do k = kcut,npro
           priorpro(k) = log10(priorpro(k))
        end do

        priorobs = priorpro

     else if (lflask) then
!       This will need to be read/computed from data (fixme)
        psurf      = 1000._r_kind
        prsiobs(1) = psurf
        pchanl     = psurf

!       Map obs pressure (prsiobs) to grid relative vertical coordinate (grdprsi)
        grdprsi = prsiobs/r10
        do k = 1,nlays+1
           if (prsimod(1) < grdprsi(k)) grdprsi(k) = prsimod(1)
           call grdcrd1(grdprsi(k), prsimod, nsig+1, -1)
        end do

!       Interpolate guess values to obs time and location (called gespro)
        call tintrp31(ges_co,  gespro(1), grdlat, grdlon, grdprsi(1), grdtime, &
                      hrdifsig, mype, nfldsig)
        call tintrp31(ges_co2, gespro(2), grdlat, grdlon, grdprsi(1), grdtime, &
                      hrdifsig, mype, nfldsig)

        priorpro = (/ zero, zero /)
!       uncert   = (/ tgasdata(1+8,i), tgasdata(2+8,i) /)
!       For now, just use ones because still testing (fixme)
        uncert   = (/ one, one /)

        avgker   = reshape((/ one, zero, zero, one /), (/ 2, 2 /))
        priorobs = priorpro

     end if

!    Create guess corresponding to obs
     gesobs = priorobs + matmul(avgker, (gespro - priorpro))

     if (ldebug) then
        print *, 'grdlat     = ', grdlat, 'grdlon     = ', grdlon
        print *, 'grdtime    = ', grdtime
        print *, 'qcflag     = ', lqcflag
        print *, 'psurf      = ', psurf
        print *, 'prsimod*10 = '
        print *, prsimod*r10
        print *, 'prsiobs    = '
        print *, prsiobs
        print *, 'grdprsi    = '
        print *, grdprsi
        print *, 'avgker     = '
        print *, avgker
        print *, 'priorpro   = '
        print *, priorpro
        print *, 'gespro     = '
        print *, gespro
        print *, 'uncert     = '
        print *, uncert
        print *, 'priorobs   = '
        print *, priorobs
     end if

     if (ltgasdiagsave .and. luse(i)) then
        id = id + 1
        idiagbuf(1,id) = mype                            ! mpi task number
        diagbuf(1,id)  = tgasdata(idlat,i)               ! lat (degree)
        diagbuf(2,id)  = tgasdata(idlon,i)               ! lon (degree)
        diagbuf(3,id)  = tgasdata(itime,i) - time_offset ! time (hours rel. to ana.)
     end if

!    Calculate innovations, perform gross checks, and accumualte statistics
     do k = 1,nchanl
        j    = ipos(k)
        ioff = nreal + k

        obs = tgasdata(ioff,i)

!       Apply log transform to co data
        if (lmopitt) then
           if (obs < tiny_r_kind) then
              obs = gesobs(k)
           else
              obs = log10(obs)
           end if

        else if (lflask .and. k == 1) then
           if (obs < tiny_r_kind) then
              obs = gesobs(k)
           else
              obs = log10(obs)
           end if

        end if

!       Compute innovation and load obs error into local array
        obsinv(k) = obs - gesobs(k)
        error(k)  = tnoise(k) * uncert(k)

        if (ldebug) print *, 'gesobs    = ', gesobs(k), 'obs       = ', obs

!       Set inverse obs error squared and ratio_errors
        if (error(k) < tol_error) then
           varinv3(k)      = one/(error(k)**2)
           ratio_errors(k) = one
        else
           varinv3(k)      = zero
           ratio_errors(k) = zero
        end if

!       Toss the obs not recommended by the data provider, and perform gross
!       check
        if (lqcflag .or. gross(k) < abs(obsinv(k))) then
           varinv3(k)      = zero
           ratio_errors(k) = zero
           if (luse(i)) stats_tgas(2,j) = stats_tgas(2,j) + one                ! # obs tossed
        end if

!       Accumulate numbers for statistics
        rat_err2 = ratio_errors(k)**2

        if (luse(i)) then
           if (     (iouse(k) == -1 .and. ltgasdiagsave)                       &
               .or. (tiny_r_kind < varinv3(k))         ) then
              omg = obsinv(k)

              if (ldebug) print *, trim(myname), ': o-g = ', omg

              stats_tgas(1,j) = stats_tgas(1,j) + one                          ! # obs
              stats_tgas(3,j) = stats_tgas(3,j) + omg                          ! (o-g)
              stats_tgas(4,j) = stats_tgas(4,j) + omg*omg                      ! (o-g)**2
              stats_tgas(5,j) = stats_tgas(5,j) + omg*omg*varinv3(k)*rat_err2  ! penalty
              stats_tgas(6,j) = stats_tgas(6,j) + obs                          ! obs

              exp_arg  = -half*varinv3(k)*omg**2
              errorinv = sqrt(varinv3(k))
              if (tiny_r_kind < pg_tgas(j) .and. tiny_r_kind < errorinv) then
                 arg       = exp(exp_arg)
                 wnotgross = one - pg_tgas(j)
                 cg_tgas   = b_tgas(j)*errorinv
                 wgross    = cg_term*pg_tgas(j)/(cg_tgas*wnotgross)
                 term      = log((arg + wgross)/(one + wgross))
                 wgt       = one - wgross/(arg + wgross)
              else
                 term      = exp_arg
                 wgt       = one
              end if

              stats_tgas(8,j) = stats_tgas(8,j) - two*rat_err2*term
              if (wgt < wgtlim) stats_tgas(9,j) = stats_tgas(9,j) + one
           end if

!          Optionally save data for diagnostics
!          (Currently only using first 5 components, but keeping the first
!          dimension of rdiagbuf at 6 to preserve compatibility with ozone
!          diagnostics)
           if (ltgasdiagsave) then
              rdiagbuf(1,k,id) = obs                    ! obs
              rdiagbuf(2,k,id) = obsinv(k)              ! obs-ges
              rdiagbuf(3,k,id) = varinv3(k)*rat_err2    ! inverse (obs error)**2
              rdiagbuf(4,k,id) = pchanl(k)              ! obs pressure
              rdiagbuf(5,k,id) = prsimod(1)*r10         ! mod surface pressure
              rdiagbuf(6,k,id) = zero                   !
           end if

           if (tiny_r_kind < rat_err2*varinv3(k)) then
              stats_tgas(7,j) = stats_tgas(7,j) + one
           end if
        end if

!       If not assimilating this observation, reset inverse variance to zero
        if (iouse(k) < 1) then
           varinv3(k) = zero
           ratio_errors(k) = zero
           rat_err2 = zero
        end if
     end do

!    Check all information for obs.  If there is at least one piece of
!    information that passed quality control, use this observation.
     ikeep = 0
     do k = 1,nchanl
        if (ldebug) then
           print *, 'ratio_errors = ', ratio_errors(k)
           print *, 'varinv3      = ', varinv3(k)
           print *, 'ratio_errors**2*varinv3 = ',                              &
                    ratio_errors(k)**2*varinv3(k)
        end if
        if (tiny_r_kind < (ratio_errors(k)**2)*varinv3(k)) ikeep = 1
     end do
end if ! (in_curbin)

     if (ldebug) print *, 'ikeep = ', ikeep

!    In principle, we want ALL obs in the diagnostics structure, but for
!    passive obs (monitoring), it is difficult to do if ltgasdiagsave
!    is not on in the first outer loop. For now we use l_may_be_passive
     if (ldebug) print *, 'l_may_be_passive = ', l_may_be_passive
     if (l_may_be_passive) then
!       Link observation to appropriate observation bin
        if (nobs_bins > 1) then
           ibin = nint(grdtime/hr_obsbin) + 1
        else
           ibin = 1
        end if
        if (ibin < 1 .or. nobs_bins < ibin) then
           write(6,*) mype, ' error nobs_bins, ibin = ', nobs_bins, ibin
        end if

if (in_curbin) then
!       Process obs that have at least one piece of information that passed qc
!       checks
        if (.not. last .and. ikeep == 1) then
           if (.not. associated(tgashead(ibin)%head)) then
              allocate(tgashead(ibin)%head, stat=istatus)
              if (istatus /= 0) write(6,*) trim(myname), ': failure to ' //    &
                                           'write tgashead'
              tgastail(ibin)%head => tgashead(ibin)%head
           else
              allocate(tgastail(ibin)%head%llpoint, stat=istatus)
              if (istatus /= 0) write(6,*) trim(myname), ': failure to ' //    &
                                           'write tgastail%llpoint'
              tgastail(ibin)%head => tgastail(ibin)%head%llpoint
           end if

           m_alloc(ibin) =  m_alloc(ibin) + 1
           my_head       => tgastail(ibin)%head
           my_head%idv   =  is
           my_head%iob   =  i

!          Increment data counter and save information used in inner loop
!          minimization (int* and stp* routines)
           tgastail(ibin)%head%luse    = luse(i)
           tgastail(ibin)%head%time    = grdtime
           tgastail(ibin)%head%nchanl  = nchanl
           tgastail(ibin)%head%npro    = npro
           tgastail(ibin)%head%obstype = obstype

           allocate(tgastail(ibin)%head%res(nchanl),         &
                    tgastail(ibin)%head%diags(nchanl),       &
                    tgastail(ibin)%head%err2(nchanl),        &
                    tgastail(ibin)%head%raterr2(nchanl),     &
                    tgastail(ibin)%head%ipos(nchanl),        &
                    tgastail(ibin)%head%avgker(nchanl,npro), &
                    tgastail(ibin)%head%avgwgt(npro,nsig),   &
                    stat=istatus)

           if (istatus /= 0) then
              write(6,*) trim(myname), ': allocation error for tgas ' //       &
                         'pointer, istatus = ', istatus
           end if

!          Set (i,j) indices of guess gridpoint that bound obs location
           call get_ij(mm1, grdlat, grdlon, tgastail(ibin)%head%ij(1),         &
                       tgastail(ibin)%head%wij(1))

           if (ldebug) then
              print *, 'ij  = ', tgastail(ibin)%head%ij
              print *, 'wij = ', tgastail(ibin)%head%wij
           end if

           do k = 1,npro
              do j = 1,nchanl
                 tgastail(ibin)%head%avgker(j,k) = avgker(j,k)
              end do 
           end do

           do l = 1,nsig
              do k = 1,npro
                 tgastail(ibin)%head%avgwgt(k,l) = avgwgt(k,l)
              end do 
           end do
        end if ! (.not. last .and. ikeep == 1)
end if ! (in_curbin)

!       Link obs to diagnostics structure
        do k = 1,nchanl
           if (.not. lobsdiag_allocated) then
              if (.not. associated(obsdiags(i_tgas_ob_type,ibin)%head)) then
                 allocate(obsdiags(i_tgas_ob_type,ibin)%head, stat=istatus)
                 if (istatus /= 0) then
                    write(6,*) trim(myname), ': failure to allocate ' //      &
                               'obsdiags ', istatus
                    call stop2(260)
                 end if
                 obsdiags(i_tgas_ob_type,ibin)%tail                           &
                    => obsdiags(i_tgas_ob_type,ibin)%head
              else
                 allocate(obsdiags(i_tgas_ob_type,ibin)%tail%next,            &
                          stat=istatus)
                 if (istatus /= 0) then
                    write(6,*) trim(myname), ': failure to allocate ' //      &
                               'obsdiags ', istatus
                    call stop2(261)
                 end if
                 obsdiags(i_tgas_ob_type,ibin)%tail                           &
                    => obsdiags(i_tgas_ob_type,ibin)%tail%next
              end if

              allocate(obsdiags(i_tgas_ob_type,ibin)%tail%muse(miter+1))
              allocate(obsdiags(i_tgas_ob_type,ibin)%tail%nldepart(miter+1))
              allocate(obsdiags(i_tgas_ob_type,ibin)%tail%tldepart(miter))
              allocate(obsdiags(i_tgas_ob_type,ibin)%tail%obssen(miter))

              obsdiags(i_tgas_ob_type,ibin)%tail%indxglb     = i
              obsdiags(i_tgas_ob_type,ibin)%tail%nchnperobs  = -99999
              obsdiags(i_tgas_ob_type,ibin)%tail%luse        = .false.
              obsdiags(i_tgas_ob_type,ibin)%tail%wgtjo       = -huge(zero)

              obsdiags(i_tgas_ob_type,ibin)%tail%muse(:)     = .false.
              obsdiags(i_tgas_ob_type,ibin)%tail%nldepart(:) = -huge(zero)
              obsdiags(i_tgas_ob_type,ibin)%tail%tldepart(:) = zero
              obsdiags(i_tgas_ob_type,ibin)%tail%obssen(:)   = zero

              n_alloc(ibin) =  n_alloc(ibin) + 1
              my_diag       => obsdiags(i_tgas_ob_type,ibin)%tail
              my_diag%idv   =  is
              my_diag%iob   =  i
              my_diag%ich   =  k
           else
              if (.not. associated(obsdiags(i_tgas_ob_type,ibin)%tail)) then
                 obsdiags(i_tgas_ob_type,ibin)%tail                            &
                    => obsdiags(i_tgas_ob_type,ibin)%head
              else
                 obsdiags(i_tgas_ob_type,ibin)%tail                            &
                    => obsdiags(i_tgas_ob_type,ibin)%tail%next
              end if
              if (obsdiags(i_tgas_ob_type,ibin)%tail%indxglb /= i) then
                 write(6,*) trim(myname), ': index error'
                 call stop2(262)
              end if
           end if

if (in_curbin) then
           obsdiags(i_tgas_ob_type,ibin)%tail%luse            = luse(i)
           obsdiags(i_tgas_ob_type,ibin)%tail%muse(jiter)     = (ikeep==1)
           obsdiags(i_tgas_ob_type,ibin)%tail%nldepart(jiter) = obsinv(k)
           obsdiags(i_tgas_ob_type,ibin)%tail%wgtjo                            &
              = varinv3(k)*ratio_errors(k)**2
 
           if (.not. last .and. ikeep == 1) then
              tgastail(ibin)%head%ipos(k)      = ipos(k)
              tgastail(ibin)%head%res(k)       = obsinv(k)
              tgastail(ibin)%head%err2(k)      = varinv3(k)
              tgastail(ibin)%head%raterr2(k)   = ratio_errors(k)**2
              tgastail(ibin)%head%diags(k)%ptr                                 &
                 => obsdiags(i_tgas_ob_type,ibin)%tail

              my_head => tgastail(ibin)%head
              my_diag => tgastail(ibin)%head%diags(k)%ptr
              if (my_head%idv /= my_diag%idv .or. my_head%iob /= my_diag%iob   &
                       .or. k /= my_diag%ich) then
                call perr(myname, 'mismatching %[head,diags]%' //              &
                                  '(idv,iob,ich,ibin) =',                      &
                                  (/is,i,k,ibin/))
                call perr(myname, 'my_head%(idv,iob,ich) =',                   &
                          (/my_head%idv,my_head%iob,k/))
                call perr(myname, 'my_diag%(idv,iob,ich) =',                   &
                          (/my_diag%idv,my_diag%iob,my_diag%ich/))
                call die(myname)
              end if
           end if

           if (ltgasdiagsave .and. lobsdiagsave .and. luse(i)) then
              jdiag = irmin
              do jj = 1,miter
                 jdiag = jdiag + 1
                 if (obsdiags(i_tgas_ob_type,ibin)%tail%muse(jj)) then
                    rdiagbuf(jdiag,k,id) =  one
                 else
                    rdiagbuf(jdiag,k,id) = -one
                 end if
              end do
              do jj = 1,miter+1
                 jdiag = jdiag + 1
                 rdiagbuf(jdiag,k,id)                                          &
                    = obsdiags(i_tgas_ob_type,ibin)%tail%nldepart(jj)
              end do
              do jj = 1,miter
                 jdiag = jdiag + 1
                 rdiagbuf(jdiag,k,id)                                          &
                    = obsdiags(i_tgas_ob_type,ibin)%tail%tldepart(jj)
              end do
              do jj = 1,miter
                 jdiag = jdiag + 1
                 rdiagbuf(jdiag,k,id)                                          &
                    = obsdiags(i_tgas_ob_type,ibin)%tail%obssen(jj)
              end do
           end if
end if ! (in_curbin)
        end do ! (k = 1,nchanl)

     else ! (.not. lmaybepassive)
if (in_curbin) then
        if (ltgasdiagsave .and. lobsdiagsave .and. luse(i)) then
           rdiagbuf(irmin+1:irdim1,1:nchanl,id) = zero
        end if
end if ! (in_curbin)
 
     end if ! (l_may_be_passive)
  end do ! (i = 1,nobs)

! If requested, write to diagnostic file
  if (ltgasdiagsave .and. 0 < id) then
     filex = obstype
     write(string,100) jiter
100  format('_',i2.2)
     diag_tgas_file =  trim(dirname) // trim(filex) // '_' // trim(dplat(is))  &
                    // string
     if (init_pass) then
       open(4, file=diag_tgas_file, form='unformatted', status='unknown',      &
            position='rewind')
     else
       open(4, file=diag_tgas_file, form='unformatted', status='old',          &
            position='append')
     end if
     iextra = npro
     if (init_pass .and. mype == mype_diaghdr(is)) then
        write(6,*) trim(myname), ': write header record for ', isis, iint,     &
                   ireal, iextra,' to file ', trim(diag_tgas_file), ' ',       &
                   ianldate
        write(4) isis, dplat(is), obstype, jiter, nchanl, ianldate, iint,      &
                 ireal, iextra
        write(4) real(pchanl,r_single), real(gross,r_single),                  &
                 real(tnoise,r_single), iouse
     end if
     write(4) id
     write(4) idiagbuf(:,1:id), diagbuf(:,1:id), rdiagbuf(:,:,1:id)
     close(4)
  end if

! Jump to this line if problem with data
135 continue

  call dtime_show(trim(myname), 'diagsave:tgas', i_tgas_ob_type)

! Clean up
  deallocate(prsiobs, grdprsi, priorpro, gespro, avgker, avgwgt)

  if (ltgasdiagsave)      deallocate(rdiagbuf)
  call final_vars_

! End of routine
  return

  contains

  subroutine check_vars_(lproceed)
     logical, intent(out) :: lproceed
     integer(i_kind) :: ivar, istatus

     lproceed = .true.

!    Check to see if required guess fields are available
     if (ihave_co) then
        call gsi_chemguess_get('var::co', ivar, istatus)
        lproceed = lproceed .and. (ivar > 0)
     end if
     if (ihave_co2) then
        call gsi_chemguess_get('var::co2', ivar, istatus)
        lproceed = lproceed .and. (ivar > 0)
     end if
  end subroutine check_vars_

  subroutine init_vars_
     real(r_kind), dimension(:,:,:), pointer :: rank3 => NULL()

     integer(i_kind)  :: ifld, ier
     character(len=5) :: varname

!    If require guess vars available, extract from bundle
     if (size(gsi_chemguess_bundle) == nfldsig) then
!       Get carbon monoxide
        if (ihave_co) then
           varname = 'co'
           call gsi_bundlegetpointer(gsi_chemguess_bundle(1), trim(varname),   &
                                     rank3, ier)
           if (ier == 0) then
              if (allocated(ges_co)) then
                 write(6,*) trim(myname), ': ges_co already incorrectly ' //   &
                            'allocated'
                 call stop2(997)
              end if
              allocate(ges_co(size(rank3,1),size(rank3,2),size(rank3,3),       &
                              nfldsig))
              ges_co(:,:,:,1) = rank3
              do ifld = 2,nfldsig
                 call gsi_bundlegetpointer(gsi_chemguess_bundle(ifld),         &
                                           trim(varname), rank3, ier)
                 ges_co(:,:,:,ifld) = rank3
              end do
           else
              write(6,*) trim(myname), ': ', trim(varname), ' not found ' //   &
                         'in chem bundle, ier = ', ier
              call stop2(998)
           end if
        end if

!       Get carbon dioxide
        if (ihave_co2) then
           varname = 'co2'
           call gsi_bundlegetpointer(gsi_chemguess_bundle(1), trim(varname),   &
                                     rank3, ier)
           if (ier == 0) then
              if (allocated(ges_co2)) then
                 write(6,*) trim(myname), ': ges_co2 already incorrectly ' //  &
                            'allocated'
                 call stop2(997)
              end if
              allocate(ges_co2(size(rank3,1),size(rank3,2),size(rank3,3),      &
                               nfldsig))
              ges_co2(:,:,:,1) = rank3
              do ifld = 2,nfldsig
                 call gsi_bundlegetpointer(gsi_chemguess_bundle(ifld),         &
                                           trim(varname), rank3, ier)
                 ges_co2(:,:,:,ifld) = rank3
              end do
           else
              write(6,*) trim(myname), ': ', trim(varname), ' not found ' //   &
                         'in chem bundle, ier = ', ier
              call stop2(998)
           end if
        end if
     else
        write(6,*) trim(myname), ': inconsistent vector sizes (nfldsig, ' //   &
                   'size(gsi_chemguess_bundle) = ', nfldsig,                   &
                   size(gsi_chemguess_bundle)
        call stop2(999)
     end if
  end subroutine init_vars_

  subroutine final_vars_
     if (allocated(ges_co))   deallocate(ges_co)
     if (allocated(ges_co2))  deallocate(ges_co2)
  end subroutine final_vars_

  subroutine zavgtgas_(f, llaypro, g, w)
     real(r_kind), intent(in   ) :: f(nsig)
     logical,      intent(in   ) :: llaypro
     real(r_kind), intent(  out) :: g(npro), w(npro,nsig)

     real(r_kind)    :: avg(nlays), wgt(nlays,nsig)

     real(r_kind)    :: dz1, dz2, delp, delz, dlay
     integer(i_kind) :: iz1, iz2
     logical         :: lincrease

     avg = zero
     wgt = zero

     do j = 1,nlays
        if (grdprsi(1) < grdprsi(nlays+1)) then
          dz1 = grdprsi(j+1)
          dz2 = grdprsi(j)
        else
          dz1 = grdprsi(j)
          dz2 = grdprsi(j+1)
        end if

        iz1 = dz1
        iz2 = dz2
        if (nsig < iz1) iz1 = nsig

!       Avoid division by zero, filled by extrapolation below
!       Think this needs some work, can create zeros if two obs levels
!       are within one model level (fixme)
        if (iz1 == iz2) cycle

        dlay = zero
        do l = iz1,iz2,-1
           delz = one
           if (l == iz1) delz =         dz1 - iz1
           if (l == iz2) delz = delz - (dz2 - iz2)

           delp = prsimod(l) - prsimod(l+1)

           dlay     = dlay   +      delp*delz
           avg(j)   = avg(j) + f(l)*delp*delz
           wgt(j,l) =               delp*delz
        end do

        avg(j) = avg(j)/dlay

        do l = 1,nsig
           wgt(j,l) = wgt(j,l)/dlay
        end do

        if (ldebug) then
           print *, 'iz1  = ', iz1,  'iz2  = ', iz2
           print *, 'dlay = ', dlay, 'avg  = ', avg(j)
        end if
     end do

!    Extrapolate to any obs layers entirely ABOVE model layers
     do k = 1,nlays
        if (avg(k) /= zero) exit
     end do
     do j = 1,k-1
        avg(j)   = avg(k)
        wgt(j,:) = wgt(k,:)
     end do
     if (ldebug) print *, trim(myname), ': extrapolating to top ', k-1,        &
                          ' layers ...'

!    Extrapolate to any obs layers entirely BELOW model layers
     do k = nlays,1,-1
        if (avg(k) /= zero) exit
     end do
     do j = k+1,nlays
        avg(j)   = avg(k)
        wgt(j,:) = wgt(k,:)
     end do
     if (ldebug) print *, trim(myname), ': extrapolating to bottom ',          &
                          nlays-(k+1), ' layers ...'

!    Is it a layer profile or not (interface profile)?
     if (llaypro) then
!       Copy layer data to outputs
        g = avg
        w = wgt
     else
!       Interpolate layer data to interface data
        g(1) = avg(1)
        do k = 2,nlays
           g(k) = 0.5*(avg(k-1) + avg(k))
        end do
        g(nlays+1) = avg(nlays)

        do l = 1,nsig
           w(1,l) = wgt(1,l)
           do k = 2,nlays
              w(k,l) = 0.5*(wgt(k-1,l) + wgt(k,l))
           end do
           w(nlays+1,l) = wgt(nlays,l)
        end do
     end if

  end subroutine zavgtgas_

end subroutine setuptgas
