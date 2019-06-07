!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS        !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: fvsvec: Main program for singular vector calculation

! !INTERFACE:
 
program fvsvec

! !USES:

#if defined ( SPMD )
  use mod_comm,only : mp_init
  use mod_comm,only : mp_exit
#else
  use m_mpif90,only : MP_init
#endif
  use m_mpif90,only : MP_comm_rank
  use m_mpif90,only : MP_comm_world
  use m_mpif90,only : MP_finalize

  use m_zeit, only : zeit_ci
  use m_zeit, only : zeit_co
  use m_zeit, only : zeit_flush

  use m_ioutil, only : luavail

  use m_die,  only : MP_die
  use m_die,  only : die

  implicit none

! !DESCRIPTION: Calculates singular vectors of FVGCM.
!
! !REVISION HISTORY:
!
!  28Sep2002  Todling   Initial code, based on FastOpt "prgadmtlm".
!  14May2007  Todling   Updated interface to svecdrv
!
!EOP
!-----------------------------------------------------------------------

! Define parameters
! -----------------
  character(len=*), parameter :: myname = 'FVSvec'
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

! Call main svec driver
! ---------------------
    call zeit_ci  ('SVDrv')
  call svdrv( ROOT, MP_comm_world )
    call zeit_co ('SVDrv')

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

end program fvsvec

