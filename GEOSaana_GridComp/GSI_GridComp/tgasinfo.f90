module tgasinfo
!$$$ module documentation block
!           .      .    .                                       .
!   module: tgasinfo
!  prgrmmr: weir        org: gmao                date: 2014-04-18
!
! abstract: This module contains variables and routines related
!           to the assimilation of trace gas observations
!           (presently, satellite based co and co2 observations)
!
! program history log:
!   2014-04-18  weir     - initial code
!   2014-11-10  weir     - renamed pob_tgas to pchanl_tgas
!   2016-05-16  weir     - adding o3 & nox support
!
! subroutines Included:
!   sub init_tgas     - set trace gas related variables to defaults
!   sub tgasinfo_read - read in trace gas info
!
! functions Included:
!
! variable Definitions:
!   def diag_tgas   - logical to turn on/off the diagnostic file (true = on)
!   def mype_tgas   - task id for writing out diagnostics
!
!   def jpch_tgas   - number of tgasinfo lines
!   def nusis_tgas  - sensor/intrument/satellite id
!   def ichanl_tgas - integer channel number of observation
!   def iuse_tgas   - integer flag to control usage of observation
!                     (-1 = don't use, 1 = use)
!   def pchanl_tgas - pressure of observation channel (hPa)
!   def gross_tgas  - gross error limit
!   def error_tgas  - observation error multiplicative inflation factor
!   def b_tgas      - variational qc: how many std devs is gross spread
!   def pg_tgas     - variational qw: probability of gross error [0,1)
!   def bias_tgas   - what it says
!   def vname_tgas  - observation variable name (co, co2, ch4, ...)
!   def vunit_tgas  - observation variable unit (mix, ppm, ppb, l10, ...)
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  use kinds,     only: r_kind, i_kind
  use constants, only: varlen => max_varname_length

  implicit none

! Set default to private
  private

! Set subroutines to public
  public :: init_tgas
  public :: tgasinfo_read

! Set passed variables to pubic
  public :: diag_tgas, mype_tgas, ihave_tgas, ntgas, tgnames
  public :: jpch_tgas, nusis_tgas, ichanl_tgas, iuse_tgas, pchanl_tgas
  public :: gross_tgas, error_tgas, b_tgas, pg_tgas, bias_tgas
  public :: vname_tgas, vunit_tgas

! Variables used to control reactive trace gases (no2 & so2)
  public :: nreact                    ! number of reactive gases (set below)
  public :: reactname                 ! names of reactive gases
  public :: tgas_minsigobs            ! minimum obs uncertainty (in 1e15 molec cm-2)
  public :: tgas_maxsigobs            ! maximum obs uncertainty (in 1e15 molec cm-2)
  public :: tgas_sigobsscal           ! obs uncertainty scale factor
  public :: tgas_minbgstrat           ! minimum bkg uncertainty in stratosphere 
  public :: tgas_minbgwater           ! minimum bkg uncertainty over water (sfc & free trop) 
  public :: tgas_minbglndpbl          ! minimum bkg uncertainty over land, pbl 
  public :: tgas_minbglndfree         ! minimum bkg uncertainty over land, free troposphere
  public :: tgas_bgscalstrat          ! background error scale factor for stratosphere
  public :: tgas_bgscalwater          ! background error scale factor over water 
  public :: tgas_bgscallndpbl         ! background error scale factor over land, pbl 
  public :: tgas_bgscallndfree        ! background error scale factor over land, free trop 
  public :: tgas_vadj                 ! vertical scale adjustment factor
  public :: tgas_hadjlevidx           ! horizontal scale boundary (level index)
  public :: tgas_hadjabove            ! horizontal scale adjustment above boundary level
  public :: tgas_hadjbelow            ! horizontal scale adjustment below boundary level
  public :: tgas_szamax               ! Maximum SZA
  public :: tgas_albmax               ! Maximum surface albedo
  public :: tgas_cldmax               ! Maximum cloud radiance fraction 

  integer(i_kind) :: mype_tgas, jpch_tgas, ntgas

  logical :: diag_tgas, ihave_tgas

  real(r_kind), allocatable :: pchanl_tgas(:), gross_tgas(:), error_tgas(:)
  real(r_kind), allocatable :: b_tgas(:), pg_tgas(:), bias_tgas(:)

  integer(i_kind),   allocatable :: ichanl_tgas(:), iuse_tgas(:)

  character(len=20),     allocatable :: nusis_tgas(:)
  character(len=varlen), allocatable :: vname_tgas(:), vunit_tgas(:)
  character(len=varlen), allocatable :: tgnames(:)

  ! Reactive trace gas stuff
  integer(i_kind), parameter  :: nreact = 2    ! 1=no2; 2=so2
  character(len=6), parameter :: reactname(nreact) = (/ 'no2', 'so2' /)
  real(r_kind)                :: tgas_minsigobs(nreact), tgas_maxsigobs(nreact), tgas_sigobsscal(nreact)
  real(r_kind)                :: tgas_minbgstrat(nreact), tgas_minbgwater(nreact)
  real(r_kind)                :: tgas_minbglndpbl(nreact), tgas_minbglndfree(nreact)
  real(r_kind)                :: tgas_bgscalstrat(nreact), tgas_bgscalwater(nreact)
  real(r_kind)                :: tgas_bgscallndpbl(nreact), tgas_bgscallndfree(nreact)
  real(r_kind)                :: tgas_vadj(nreact), tgas_szamax(nreact), tgas_albmax(nreact), tgas_cldmax(nreact)
  integer(i_kind)             :: tgas_hadjlevidx(nreact)
  real(r_kind)                :: tgas_hadjabove(nreact), tgas_hadjbelow(nreact)

contains
  
  subroutine init_tgas
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    init_tgas    initialize parameters for tgas data
!     prgmmr:    weir        org: gmao                date: 2014-04-18
!
! abstract:      This routine sets default values for variables used in 
!                the tgas processing routines
!
! program history log:
!   201?-??-??  weir     - initial code
!   2018-06-20  weir     - eliminating tracer name dependence
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

    use mpimod,            only: npe, mype
    use mpeu_util,         only: die
    use gsi_chemguess_mod, only: gsi_chemguess_get

    implicit none

    character(len=*), parameter :: myname = 'init_tgas'

    integer(i_kind) :: ier

    diag_tgas = .true.                              ! default is to generate tgas diagnostic file
    mype_tgas = max(0, npe-6)                       ! mpi task to write tgas summary report
    jpch_tgas = 0                                   ! number of entries read from tgasinfo 

    call gsi_chemguess_get('ghg', ntgas, ier)
    if (ier /= 0) then
!      write(0,*) myname // ': chemguess_get #1 failed, error = ', ier
!      call die(myname)

!      Report as a warning since some versions incorrectly report an error
       if (mype == 0) then
          write(6,*) myname // ': chemguess_get #1 reported error = ', ier
          write(6,*) myname // ': proceeding anyway with ntgas = ', ntgas
       end if
    end if

!   Some error checking
    ihave_tgas = (ntgas > 0)
    if (.not. ihave_tgas) return

    if (allocated(tgnames)) then
       write(0,*) myname // ': tgnames already allocated'
       call die(myname)
    end if

    allocate(tgnames(ntgas))
    call gsi_chemguess_get('ghg', tgnames, ier)

    if (ier /= 0) then
       write(0,*) myname // ': chemguess_get #2 failed, error = ', ier
       call die(myname)
    end if
  end subroutine init_tgas
  
  subroutine tgasinfo_read
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    tgasinfo_read      read tgas information file
!     prgmmr:    weir        org: gmao                date: 2014-04-18
!
! abstract:  This routine reads the tgas information file tgasinfo
!
! program history log:
!   2010-04-18  weir     - initial code
!
!   input argument list:
!     mype - mpi task id
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

    use mpimod,    only: mype
    use mpeu_util, only: die
    use obsmod,    only: iout_tgas

    implicit none

    character(len=*), parameter :: myname = 'tgasinfo_read'

    integer(i_kind)    :: lunin, j, k, istat, nlines
    character(len=1)   :: cflg
    character(len=120) :: crecord

    data lunin / 47 /

!   Determine number of entries in trace gas information file
    open(lunin, file='tgasinfo', form='formatted')
    j      = 0
    nlines = 0
    addlines:  do 
       read(lunin,100,iostat=istat,end=123) cflg, crecord
       if (istat /= 0) exit
       nlines = nlines + 1
!      Skip given any reasonable comment indicator
       if (cflg == '!' .or. cflg == '#' .or. cflg == '%') cycle
       j = j + 1
    end do addlines
123 continue
    if (0 < istat) then
       write(0,*) myname // ': *** ERROR: error reading tgasinfo ****, ' //    &
                  'error = ', istat
       close(lunin)
       call die(myname)
    end if
    jpch_tgas = j

!   Allocate arrays to hold tgas information
    allocate(nusis_tgas(jpch_tgas), ichanl_tgas(jpch_tgas),                    &
              iuse_tgas(jpch_tgas), pchanl_tgas(jpch_tgas),                    &
             gross_tgas(jpch_tgas),  error_tgas(jpch_tgas),                    &
                 b_tgas(jpch_tgas),     pg_tgas(jpch_tgas),                    &
              bias_tgas(jpch_tgas),  vname_tgas(jpch_tgas),                    &
              vunit_tgas(jpch_tgas))

!   All mpi tasks open and read tgas information file.
!   Task mype_tgas writes information to tgas runtime file
    if (mype == mype_tgas) then
       open(iout_tgas)
       write(iout_tgas,*) myname // ': jpch_tgas = ', jpch_tgas
       write(iout_tgas,120)
    end if
    rewind(lunin)
    j = 0
    do k = 1,nlines
       read(lunin,100) cflg, crecord
!      Skip given any reasonable comment indicator
       if (cflg == '!' .or. cflg == '#' .or. cflg == '%') cycle
       j = j + 1
       read(crecord,*) nusis_tgas(j), ichanl_tgas(j), iuse_tgas(j),            &
                       pchanl_tgas(j), gross_tgas(j), error_tgas(j),           &
                       b_tgas(j), pg_tgas(j), bias_tgas(j),                    &
                       vname_tgas(j), vunit_tgas(j)
       if (mype == mype_tgas) then
          write(iout_tgas,130) j, ichanl_tgas(j), nusis_tgas(j), iuse_tgas(j), &
                               pchanl_tgas(j), gross_tgas(j), error_tgas(j),   &
                               b_tgas(j), pg_tgas(j), bias_tgas(j),            &
                               vname_tgas(j), vunit_tgas(j)
       end if
    end do
    close(lunin)
    if (mype == mype_tgas) close(iout_tgas)

100 format(a1, a120)

120 format(t5, 'jc', t10, 'kc', t13, 'isis', t30, 'iouse', t41, 'pchanl',      &
           t52, 'gross', t60, 'inflate', t72, 'qc: b', t81, 'qc: pg',          &
           t93, 'bias', t98, 'chem', t105, 'units')
130 format(2x, i4, 1x, i4, 1x, a20, i2, 1x, f11.2, 1x, f9.3, 1x, f9.5, 1x,     &
           f9.3, 1x, 2(f9.5,1x), a6, 1x, a6)

!   Successful read, return to calling routine
    return
  end subroutine tgasinfo_read
  
end module tgasinfo
