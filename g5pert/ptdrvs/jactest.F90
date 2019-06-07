!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS        !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: jactest: Main program for testing tlm/adm jacobian

! !INTERFACE:
 
program jactest

! !USES:

#ifdef SPMD
  use mod_comm,only : mp_init
  use mod_comm,only : gid
  use mod_comm,only : mp_exit
#define CPP_PRT_PREFIX  if(gid == 0)
#else
  use m_mpif90,only : MP_init
  use m_mpif90,only : MP_finalize
#define CPP_PRT_PREFIX
#endif
  use m_mpif90,only : MP_comm_rank
  use m_mpif90,only : MP_comm_world

  use m_ioutil, only : luavail

  use m_die,  only : MP_die
  use m_die,  only : die
  use timingModule, only: timing_init,timing_prt,timing_on,timing_off

  implicit none

! !DESCRIPTION: Calculates singular vectors of FVGCM.
!
! !REVISION HISTORY:
!
!  24Mar2003  Errico   Initial code, based on FastOpt "prgadmtlm".
!  11Jan2007  Todling  Merged w/ Oloso SPMD mods: didn't take all.
!
!EOP
!-----------------------------------------------------------------------

! Define parameters
! -----------------
  character(len=*), parameter :: myname = 'jactest'
  integer, parameter :: ROOT = 0

  integer n    ! number of independent (initial state) variables
  integer m    ! number of   dependent (final   state) variables
  integer ierr
  integer myID
  integer lu

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

  call timing_init
  call timing_on('All')

! Get number of the independent and dependent variables
! -----------------------------------------------------
  call setfunc( n, m )
	if(n/=m) call die(myname,'invalid dimensions(n,m)')

! Call jacobian test
! ------------------
  call jacobian_test ( n, m )

! Echo run status and close MPI
! -----------------------------
  if(myid==ROOT) print *, myname, ': Normal Execution.'
  lu = luavail()
  open(lu,file='status.log')
  write(lu,*) 'Normal Execution.'
  close(lu)

  call timing_off('All')

#if defined ( SPMD )
  call mp_exit
#else
  call MP_finalize(ierr)
	if(ierr/=0) call MP_die(myname,'MP_finalized()',ierr)
#endif

end program jactest

