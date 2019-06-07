!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS        !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: qsens: Main program for tracer sensitivity studies

! !INTERFACE:
 
program qsens

! !USES:

#if defined ( SPMD )
  use mod_comm,only : mp_init
  use mod_comm,only : mp_exit
#else
  use m_mpif90,only : MP_init
#endif
  use m_mpif90,only : MP_finalize
  use m_mpif90,only : MP_comm_rank
  use m_mpif90,only : MP_comm_world

  use m_ioutil, only : luavail

  use m_zeit

  use m_die,  only : MP_die
  use m_die,  only : die

  use m_gQadm, only : gQadm_init
  use m_gQadm, only : gQadm_set
  use m_gQadm, only : gQadm_run
  use m_gQadm, only : gQadm_wrap
  use m_gQadm, only : gQadm_clean

  implicit none

! !DESCRIPTION: Calculates sensitivity vectors
!
! !REVISION HISTORY:
!
!  27Sep2007  Todling   Initial code, based on fvsens.
!
!EOP
!-----------------------------------------------------------------------

! Define parameters
! -----------------
  character(len=*), parameter :: myname = 'QSens'
  integer, parameter :: ROOT = 0

  character(len=255) :: oprefix
  integer n    ! number of independent (initial state) variables
  integer m    ! number of   dependent (final   state) variables
  integer ierr
  integer myID
  integer lu
  integer vectype, freq
  logical oneshot

! MPI Initialization
! ------------------
#if defined ( SPMD )
  call mp_init
#else
  call MP_init(ierr)
	if(ierr/=0) call MP_die(myname,'MP_init()',ierr)
#endif
  call MP_comm_rank(MP_comm_world,myID,ierr)
	if(ierr/=0) call MP_die(myname,'MP_comm_rank()',ierr)

  call init_ ( vectype, oneshot, freq, oprefix )

    call zeit_ci  ('qsens')
  call gQadm_init  ( ierr, verbose=.true., oneshot=oneshot, freq=freq )
  call gQadm_set   ( ierr, vectype=vectype )
      if(ierr/=0) call MP_die(myname,'set()',ierr)
  call gQadm_run   ( ierr )
  if(trim(oprefix)=='x') then
     call gQadm_wrap  ( ierr )
  else
     call gQadm_wrap  ( ierr, fname=oprefix )
  endif
  call gQadm_clean ( ierr )
    call zeit_co  ('qsens')

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
      subroutine init_ ( vectype, oneshot, freq, oprefix )

! !INPUT/OUTPUT PARAMETERS:
!
      implicit NONE
      integer, intent(out)  :: vectype ! type of input perturbation (g4/g5)
      logical, intent(out)  :: oneshot ! defines where to run TLM in single shot or not
      integer, intent(out)  :: freq    ! freq of output perturbation
      character(len=*), intent(out) :: oprefix ! name of output file
!
!
! !REVISION HISTORY:
!
!	27Sep2007 Todling  As in fvsens
!
!EOP

      character*4, parameter :: myname = 'init_'

      integer iret, i, lv, iarg, argc, iargc
      character*255 argv,str
      character*10  SS

!     Set defaults
!     ------------
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
      end do
 111  continue

      return

      end subroutine init_

      subroutine usage_()
      integer ierr
      if(myid==ROOT) then
      print *
      print *, 'Usage:  fvpert [-h] [-nlm] [-dotp] [-g5] -[inc]'
      print *
      print *, 'where'
      print *
      print *,   '-h         echoes this Usage statement'
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

end program qsens

