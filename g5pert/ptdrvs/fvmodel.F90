!-------------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS        !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: fvmodel: Main program to model alone

! !INTERFACE:
 
program fvmodel

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

  use m_zeit

  use m_die,  only : MP_die
  use m_die,  only : die

  implicit none

! !DESCRIPTION: Calculates singular vectors of FVGCM.
!
! !REVISION HISTORY:
!
!  28Sep2002  Todling   Initial code, based on FastOpt "prgadmtlm".
!  28Sep2007  Todling   Upgraded to lasted interface.
!
!EOP
!-----------------------------------------------------------------------

! Define parameters
! -----------------
  character(len=*), parameter :: myname = 'main'
  integer, parameter :: ROOT = 0

  integer n    ! number of independent (initial state) variables
  integer m    ! number of   dependent (final   state) variables
  integer ierr
  integer myID

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

! Integrate model forward ...
! ---------------------------
    call zeit_ci ('nlm')
  call nlmodel ()
    call zeit_co ('nlm')

! Summarize timing and Close MPI
! ------------------------------
  if(myid==ROOT) call zeit_flush(6)

#if defined ( SPMD )
  call mp_exit
#else
  call MP_finalize(ierr)
        if(ierr/=0) call MP_die(myname,'MP_finalized()',ierr)
#endif


end program fvmodel

