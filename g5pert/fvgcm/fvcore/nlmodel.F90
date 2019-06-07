!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: nlmodel:  Integrates nonlinear model
!
! !INTERFACE:

subroutine nlmodel

! !USES:

  use precision

  use stepon, only : job
  use stepon, only : nymd
  use stepon, only : nhms
  use stepon, only : fvpsasdt
  use stepon, only : nstep
  use stepon, only : ptop
  use stepon, only : ks
  use stepon, only : ak
  use stepon, only : bk
  use stepon, only : ts
  use stepon, only : oro
  use stepon, only : stepon_set
  use stepon, only : stepon_finalize

  use prognostics, only : dyn_prog

  use control,only : mod2cont
  use dependent,only : model2dependent

  use m_iostate, only : PutState

  use m_die, only : MP_die

  implicit none

! !INPUT PARAMETERS:

! integer, intent(in) :: n  ! dimension of state vector
! integer, intent(in) :: m  ! dummy

! !DESCRIPTION: Integrates nonlinear model.
!
! !REVISION HISTORY:
!
!  28Sep2002  Todling   Initial code, based on FastOpt "prgadmtlm".
!  26Jan2007  Todling   Added hour to output state filename 
!  14May2007  Todling   Introduced dyn_prog; global change.
!  28Sep2007  Todling   Updated to latest inferfaces.
!
!EOP
!-----------------------------------------------------------------------

! Define parameters
! -----------------
  character(len=*), parameter :: myname = 'nlmodel'
  integer,          parameter :: prec = 0  ! 32-bit file

  integer n,m,ierr
  character(len=255)  fname

  type(dyn_prog) :: prog

  call setfunc ( n, m, prog )			! Read  in  initial state
  call stepon_set ( prog )			! Read  in  initial state
  write(fname,'(2a,i8.8,a,i2.2,a)') trim(job), '.prog.eta.', nymd, '_', nhms/10000, 'z.hdf'

  if(fvpsasdt/=0) &
  call PutState ( fname, &
                  job, nymd, nhms, prog, fvpsasdt, nstep, &
                  ak, bk, Ts, oro, prec=prec, stat=ierr )
	if(ierr/=0) call MP_die(myname,'error from PutState(1)',ierr)

! Integrate model forward ...
! ---------------------------
  call func ( prog )

  write(fname,'(2a,i8.8,a,i2.2,a)') trim(job), '.prog.eta.', nymd, '_', nhms/10000, 'z.hdf'

  if(fvpsasdt/=0) &
  call PutState ( fname, &
                  job, nymd, nhms, prog, fvpsasdt, nstep, &
                  ak, bk, Ts, oro, prec=prec, stat=ierr )
 	if(ierr/=0) call MP_die(myname,'error from PutState(2)',ierr)

! Wrap it all up
! --------------
  call stepon_finalize ( prog )

end subroutine nlmodel

