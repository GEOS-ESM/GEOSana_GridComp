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
!
! subroutines Included:
!   sub init_tgas     - set trace gas related variables to defaults
!   sub tgasinfo_read - read in trace gas info
!
! functions Included:
!
! variable Definitions:
!   def diag_tgas   - logical to turn off or on the diagnostic tgas file
!                     (true = on)
!   def jpch_tgas   - was number of (levels+1) * number of satellites, but
!                     probably just 1 now
!   def mype_tgas   - task id for writing out radiance diagnostics
!   def pchanl_tgas - pressure of observation channel (hPa)
!   def gross_tgas  - gross error limit
!   def error_tgas  - tgas observation error (total column)
!   def nusis_tgas  - sensor/intrument/satellite id
!   def nupro       - integer profile number of tgas observation
!   def iuse_tgas   - integer flag to control usage of tgas data
!                     (-1 = don't use, 1 = use)
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

  use kinds, only: r_kind, i_kind

  implicit none

! Set default to private
  private

! Set subroutines to public
  public :: init_tgas
  public :: tgasinfo_read

! Set passed variables to pubic
  public :: jpch_tgas, diag_tgas, nusis_tgas, iuse_tgas, b_tgas, pg_tgas
  public :: pchanl_tgas, gross_tgas, error_tgas, mype_tgas, nupro, ihave_tgas
  public :: ihave_co, ihave_co2, ihave_ch4

  integer(i_kind) :: mype_tgas, jpch_tgas

  logical :: diag_tgas, ihave_tgas, ihave_co, ihave_co2, ihave_ch4

  real(r_kind),      allocatable :: pchanl_tgas(:), gross_tgas(:), error_tgas(:)
  real(r_kind),      allocatable :: pg_tgas(:), b_tgas(:)

  integer(i_kind),   allocatable :: nupro(:), iuse_tgas(:)
  character(len=20), allocatable :: nusis_tgas(:), tchem_tgas(:)

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
!   2010-04-18  weir     - initial code
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

    use mpimod,            only: npe                 ! the number of mpi tasks
    use gsi_chemguess_mod, only: gsi_chemguess_get

    implicit none

    integer(i_kind) :: itgas, ier

    jpch_tgas  = 0                                   ! number of enteries read from tgasinfo 
    diag_tgas  = .true.                              ! default is to generate tgas diagnostic file
    mype_tgas  = max(0, npe-6)                       ! mpi task to write tgas summary report

    ihave_tgas = .false.
    call gsi_chemguess_get('var::co',  itgas, ier)
    ihave_co   = (itgas > 0)
    ihave_tgas = ihave_tgas .or. (itgas > 0)
    call gsi_chemguess_get('var::co2', itgas, ier)
    ihave_co2  = (itgas > 0)
    ihave_tgas = ihave_tgas .or. (itgas > 0)
    call gsi_chemguess_get('var::ch4', itgas, ier)
    ihave_ch4  = (itgas > 0)
    ihave_tgas = ihave_tgas .or. (itgas > 0)

  end subroutine init_tgas
  
  subroutine tgasinfo_read
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    tgasinfo_read      read tgas information file
!     prgmmr:    weir        org: gmao                date: 2014-04-18
!
! abstract:  This routine reads the tgas information file, global_tgasinfo.txt
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

    use mpimod, only: mype
    use obsmod, only: iout_tgas

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
       if (cflg == '!') cycle
       j = j + 1
    end do addlines
123 continue
    if (istat > 0) then
       write(6,*) trim(myname), ': ***ERROR*** error reading tgasinfo, ' //    &
                  'istat = ', istat
       close(lunin)
       write(6,*) trim(myname), ': stop program execution'
       call stop2(79)
    end if
    jpch_tgas = j

!   Allocate arrays to hold tgas information
    allocate(nusis_tgas(jpch_tgas), nupro(jpch_tgas), iuse_tgas(jpch_tgas),    &
             pchanl_tgas(jpch_tgas), gross_tgas(jpch_tgas),                    &
             error_tgas(jpch_tgas), b_tgas(jpch_tgas), pg_tgas(jpch_tgas),     &
             tchem_tgas(jpch_tgas))

!   All mpi tasks open and read tgas information file.
!   Task mype_tgas writes information to tgas runtime file
    if (mype == mype_tgas) then
       open(iout_tgas)
       write(iout_tgas,*) trim(myname), ': jpch_tgas = ', jpch_tgas
    end if
    rewind(lunin)
    j = 0
    do k = 1,nlines
       read(lunin,100) cflg, crecord
       if (cflg == '!') cycle
       j = j + 1
       read(crecord,*) nusis_tgas(j), nupro(j), iuse_tgas(j), pchanl_tgas(j),  &
                       gross_tgas(j), error_tgas(j), b_tgas(j), pg_tgas(j),    &
                       tchem_tgas(j)
       if (mype == mype_tgas) then
          write(iout_tgas,130) j, nusis_tgas(j), tchem_tgas(j), nupro(j),      &
                               iuse_tgas(j), pchanl_tgas(j), gross_tgas(j),    &
                               error_tgas(j), b_tgas(j), pg_tgas(j)
       end if
    end do
    close(lunin)
    if (mype == mype_tgas) close(iout_tgas)

100 format(a1, a120)
130 format(i3, 1x, a9, 1x, a9, 1x, ' pro = ', i4, ' use = ', i2, ' pch = ',    &
           f9.3, ' gross = ', f7.3, ' error = ', f7.3, ' b_tgas = ', f7.3,     &
           ' pg_tgas = ', f7.3)

!   Successful read, return to calling routine
    return
  end subroutine tgasinfo_read
  
end module tgasinfo
