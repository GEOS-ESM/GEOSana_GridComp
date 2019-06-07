!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS        !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: qpert: Main program for tracer TLM integration

! !INTERFACE:
 
program qpert

! !USES:

#ifdef SPMD
  use mod_comm,only : mp_init
  use mod_comm,only : mp_exit
#else
  use m_mpif90,only : MP_init
#endif
  use m_mpif90,only : MP_finalize
  use m_mpif90,only : MP_comm_rank
  use m_mpif90,only : MP_comm_world

  use timingmodule, only : timing_on
  use timingmodule, only : timing_off
  use timingmodule, only : timing_prt

  use m_gQtlm, only : gQtlm_init
  use m_gQtlm, only : gQtlm_set
  use m_gQtlm, only : gQtlm_run
  use m_gQtlm, only : gQtlm_dotp
  use m_gQtlm, only : gQtlm_wrap
  use m_gQtlm, only : gQtlm_clean

  use m_ioutil, only : luavail

  use m_die,  only : MP_die
  use m_die,  only : die

  use m_zeit

  implicit none

! !DESCRIPTION: Calculates singular vectors of FVGCM.
!
! !REVISION HISTORY:
!
!  25Sep2007  Todling   Initial code.
!
!EOP
!-----------------------------------------------------------------------

! Define parameters
! -----------------
  character(len=*), parameter :: myname = 'QPert'
  integer, parameter :: ROOT = 0

  character(len=255) :: oprefix

  integer n    ! number of independent (initial state) variables
  integer m    ! number of   dependent (final   state) variables
  integer ierr
  integer myID,seed
  integer lu
  integer vectype,freq
  logical dotp,oneshot

! MPI Initialization
! ------------------
#ifdef SPMD
  call mp_init
#else
  call MP_init(ierr)
	if(ierr/=0) call MP_die(myname,'MP_init()',ierr)
#endif
  call MP_comm_rank(MP_comm_world,myID,ierr)
	if(ierr/=0) call MP_die(myname,'MP_comm_rank()',ierr)

  call init_ ( dotp, seed, vectype, oneshot, freq, oprefix )

  call gQtlm_init  ( ierr, verbose=.true., oneshot=oneshot, freq=freq )
  call gQtlm_set   ( ierr, vectype=vectype )
	if(ierr/=0) call MP_die(myname,'set()',ierr)
  if ( dotp ) then
       call gQtlm_dotp  ( seed, ierr )
  else
       call gQtlm_run   ( ierr )
  endif
  if(trim(oprefix)=='x') then
     call gQtlm_wrap  ( ierr )
  else
     call gQtlm_wrap  ( ierr, fname=oprefix )
  endif
  call gQtlm_clean ( ierr )

! Summarize timing
! ----------------
  if(myid==ROOT) call zeit_flush(6)

! Echo run status and close MPI
! -----------------------------
  if(myid==ROOT) print *, myname, ': Normal Execution.'
  lu = luavail()
  open(lu,file='status.log')
  write(lu,*) 'Normal Execution.'
  close(lu)
#if defined ( SPMD )
  call mp_exit
#else
  call MP_finalize(ierr)
	if(ierr/=0) call MP_die(myname,'MP_finalized()',ierr)
#endif


contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: init: initialize odsshuffle
!
! !DESCRIPTION:
!
! !INTERFACE:
!
      subroutine init_ ( dotp, seed, vectype, oneshot, freq, oprefix )

! !INPUT/OUTPUT PARAMETERS:
!
      implicit NONE
      logical, intent(out)  :: dotp    ! when .t., performs test of ADM/TLM
      integer, intent(out)  :: seed    ! random number seed
      integer, intent(out)  :: vectype ! type of input perturbation (g4/g5)
      logical, intent(out)  :: oneshot ! defines where to run TLM in single shot or not
      integer, intent(out)  :: freq    ! freq of output perturbation
      character(len=*), intent(out) :: oprefix ! name of output file
!
!
! !REVISION HISTORY:
!
!	25Sep2007 Todling  Initial code, based on fvpert.F90
!
!EOP

      character*4, parameter :: myname = 'init_'

      integer iret, i, lv, iarg, argc, iargc
      character*255 argv,str
      character*10  SS

!     Set defaults
!     ------------
      dotp = .false.
      seed = 0
      vectype = 4
      oneshot = .true.
      oprefix =  'x'
      freq  = -1

!     Parse command line
!     ------------------
      argc =  iargc()
      if ( argc .lt. 1 ) return
      iarg = 0
      do i = 1, 32767
         iarg = iarg + 1
         if ( iarg .gt. argc ) go to 111
         call GetArg ( iArg, argv )
         if (index(argv,'-h' ) .gt. 0 ) then
             call usage_()
         endif
         if (index(argv,'-dotp' ) .gt. 0 ) then
            dotp = .true.
         endif
         if (index(argv,'-g5' ) .gt. 0 ) then
            vectype = 5
         endif
         if (index(argv,'-inc' ) .gt. 0 ) then
            oneshot = .false.
         endif
         if (index(argv,'-freq') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage_()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            read(SS,*) freq
         endif
         if (index(argv,'-o' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage_()
            iarg = iarg + 1
            call GetArg ( iArg, oprefix )
         endif
         if (index(argv,'-seed' ) .gt. 0 ) then
            iarg = iarg + 1
            call GetArg ( iarg, str )
            read(str,*) seed
         endif
      end do
 111  continue

      return

      end subroutine init_

      subroutine usage_()
      integer ierr
      if(myid==ROOT) then
      print *
      print *, 'Usage:  fvpert [-h] [-dotp] [-g5] [-inc]'
      print *
      print *, 'where'
      print *
      print *,   '-h         echoes this Usage statement'
      print *,   '-dotp      performs dop-product check for ADM/TLM validity'
      print *,   '-seed #    seed to random # generator; seed=-1, use input perturbation'
      print *,   '-g5        to be specified when input perturbation is GEOS-5-like'
      print *,   '-inc       runs TLM incrementally (step-by-step)'
      print *,   '-freq HHMMSS frequency for output perturbation'
      print *,   '-o PREFIXFN  prefix of filename for propagated perturbation'
      print *
      endif
      call MP_finalize(ierr)
	if(ierr/=0) call MP_die(myname,'MP_finalized()',ierr)
      stop
      end subroutine usage_

end program qpert
