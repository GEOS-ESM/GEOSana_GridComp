!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_gQadm --- GEOS Tracer Adjoint Model Module
!
! !INTERFACE:
!
  module m_gQadm

! !USES:

  use precision
  use mod_comm, only : gid
  use stepon,  only : stepon_set
  use stepon,  only : stepon_finalize
  use stepon,  only : job
  use stepon,  only : nymd
  use stepon,  only : nhms
  use stepon,  only : nymde
  use stepon,  only : nhmse
  use stepon,  only : fvpsasdt
  use stepon,  only : pdt
  use stepon,  only : nstep
  use stepon,  only : ninner
  use stepon,  only : nouter
  use stepon,  only : ptop
  use stepon,  only : imr,jfirst,jlast
  use stepon,  only : ks
  use stepon,  only : ak
  use stepon,  only : bk
  use stepon,  only : ts
  use stepon,  only : oro

  use prognostics_q, only : prognostics_q_initial
  use prognostics_q, only : prognostics_q_final
  use prognostics_q, only : prognostics_q_dup
  use prognostics_q, only : dynq

  use prognostics,   only : dyn_prog

  use m_model_tad, only : initial_tad
  use m_model_tad, only : qmodel_tad
  use m_model_tad, only : final_tad

  use m_iostate, only : PutState
  use m_iostate, only : getstate_init
  use m_iostate, only : getstate
  use m_iostate, only : getstate_clean

  use m_trjphys, only : physdrv1_get_init
  use m_trjphys, only : physdrv1_get_all
  use m_trjphys, only : physdrv1_get_clean

  use m_trajmng, only : GetPert
  use m_trajmng, only : PutPert

  use m_poles        

  use m_mpout,only : mpout_log

  implicit none

  PRIVATE

  PUBLIC gQadm_init
  PUBLIC gQadm_set
  PUBLIC gQadm_run
  PUBLIC gQadm_wrap
  PUBLIC gQadm_clean

! !REMARKS:
!
!  Because of the timing routines in the GCM this backward integration
!  for the sensitivity operator is done as if the model were going
!  forward in the time over the time period of interest. Inside the
!  adjoint, the clock is ticked properly in order to execute the
!  operations backwards as they should be executed. The consequece
!  of this is that the "initial" and "final" dates/times below are
!  in reality the reverse.
!
! !REVISION HISTORY:
!
!  27Sep2007 Todling - Initial code; based full ADM
!  06Mar2009 Todling - .hdf to .nc4
!
!EOP
!-------------------------------------------------------------------------

! Allocatable arrays
! ------------------
  type(dyn_prog)              :: prog   ! individual prognostic   variables
  type(dynq)                  :: xpert  ! individual perturbation variables
  type(dynq)                  :: ypert  ! individual perturbation variables

  integer :: ndim                           ! dimension of state-vector and perturbation

  integer :: nymdi                          ! initial date of integration
  integer :: nhmsi                          ! initial time of integration
  integer :: nymdf                          ! final   date of integration
  integer :: nhmsf                          ! final   time of integration
 
  character(len=*), parameter :: myname = 'm_gQadm'

  character(len=*), parameter :: gradfname = 'Jqgradf.eta.nc4'
  character(len=*), parameter :: fsensname = 'qsens.eta'
  character(len=*), parameter :: RCfile    = 'tracer.rc'

  integer, save      :: nvec_ad    = 2
  logical, save      :: memtraj    = .true.   ! load trajectory to memory
  logical, save      ::  g5pert    = .false.
  logical, save      ::  oneshot_  = .false.
  integer, save      :: freq_      = -1       ! frequency of TLM output (HHMMSS)
                                                                                                                           
  interface gQadm_init  ; module procedure init_  ; end interface
  interface gQadm_set   ; module procedure set1_  ; end interface
  interface gQadm_run   ; module procedure run_   ; end interface
  interface gQadm_wrap  ; module procedure wrap1_ ; end interface
  interface gQadm_clean ; module procedure clean_ ; end interface

  CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_ --- initialize atmospheric model adjoint
!
! !INTERFACE:

  subroutine init_ ( stat, &
                     n, verbose, oneshot, freq  )   ! optional

! !USES:

  implicit none

! !INPUT PARAMETERS:

   logical, optional, intent(in) :: verbose
   logical, optional, intent(in)  :: oneshot ! run step-by-step or not
   integer, optional, intent(in)  :: freq    ! freq of output during TLM (HHMMSS)
                                                                                                                           
! !OUTPUT PARAMETERS:

   integer,           intent(out) :: stat
   integer, optional, intent(out) :: n

! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  27Sep2007 Todling - Initial code
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  character(len=*), parameter :: myname_ = '::init_'

  character(len=255) words
  integer            i, n1, n2, ierr

   stat = 0

! Get number of the independent and dependent variables
! -----------------------------------------------------
  call setfunc( n1, n2, prog )
        if(n1/=n2)then
           stat = 90
           call mpout_log(myname_,' invalid dimensions(n,m)')
           return
         endif
  ndim = n1

! Initialize options, model, and set the first guess of the control variables
! ---------------------------------------------------------------------------
  call stepon_set ( prog )			! read in (initial nonlinear) state

! Store timings for output pert file tagging
! ------------------------------------------
  nymdi = nymd ;   nhmsi = nhms
  nymdf = nymde;   nhmsf = nhmse
    write(words,'(i8.8,1x,i6.6)') nymdi, nhmsi
    call mpout_log(myname_,': integration starts at: '//words)
    write(words,'(i8.8,1x,i6.6)') nymdf, nhmsf
    call mpout_log(myname_,': integration   ends at: '//words)

! Initialize dynamics trajectory handle
! -------------------------------------
  call getstate_init ( nymd, nhms, memtrj=memtraj, verbose=verbose )
  call getstate      ( nymd, nhms )

! Initialize physics trajectory handle
! -------------------------------------
  call physdrv1_get_init ( nymd, nhms, memphys=memtraj, verbose=verbose )
  call physdrv1_get_all  ( nymd, nhms )

! Allocate individual perturbation variables
! ------------------------------------------
  call prognostics_q_initial ( xpert )
  call prognostics_q_initial ( ypert )

! Reset inner and outer iterations
! --------------------------------
  ninner = nouter
  nouter = 1
  if(gid==0) print *,' No Checkpointing: Swap inner/outer iter'
  if(gid==0) print *,' nouter = ', nouter
  if(gid==0) print *,' ninner = ', ninner

! Return dim if requested
! -----------------------
  if ( present(n) ) then
       n = ndim
  endif

  if ( present(oneshot) ) then
       oneshot_ = oneshot
  endif
  if(oneshot_) freq_ = fvpsasdt
  if ( present(freq) ) then
       freq_ = freq
  endif

  end subroutine init_ 
!EOC

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: set1_ --- set atmospheric model adjoint (case 1: pert in file)
!
! !INTERFACE:

  subroutine set1_ ( stat, &
                     fname, vectype, yp ) ! optionals

! !USES:

  implicit none

! !INPUT PARAMETERS:

  character(len=*), optional, intent(in) :: fname   ! filename containing input perturbation
  integer,          optional, intent(in) :: vectype ! specifies g4/g5 perturbation
                                                                                                                           
! !OUTPUT PARAMETERS:

  type(dynq),     optional, intent (out) :: yp      ! perturbation structure

  integer,                  intent (out) :: stat    ! return error code

! !DESCRIPTION: 
!
! !TO DO:
!   1) unwire perturbation filename (get from rc file)
!
! !REVISION HISTORY:
!
!  27Sep2007 Todling - Initial code
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  character(len=*), parameter :: myname_ = 'set1_'

  character(len=255)  fnpert
  integer    nymdp        ! date of perturbation (not used)
  integer    nhmsp        ! time of perturbation (not used)
  integer    n, ierr

  stat = 0

! Get perturbation field
! ----------------------
!_RT  nymdp = 0; nhmsp = 0
  nymdp = nymdf; nhmsp = nhmsf
  if ( present(fname) ) then
       fnpert = trim(fname)
  else
       fnpert = gradfname
  endif
  if ( present(vectype) ) then
       if (vectype==5) then
          g5pert = .true.
       else
          g5pert = .false.
       endif
  endif

  call GetPert ( fnpert, nymdp, nhmsp, xpert, stat=ierr, RCfile=RCfile )
      if(ierr/=0)then
          stat = 90
          call mpout_log(myname_,'Error retrieving perturbation')
         return
      endif

! If so, return perturbation in user array
! ----------------------------------------
  if ( present(yp) ) then
       call prognostics_q_dup ( xpert, yp )
  endif

  end subroutine set1_ 
!EOC

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: run_ --- integrate atmospheric model adjoint
!
! !INTERFACE:

  subroutine run_ ( stat,    &
                    xp, yp ) ! optionals

! !USES:

  implicit none

! !INPUT PARAMETERS:

  type(dynq), optional, intent(in)  :: yp     ! initial perturbation for  ADM integration

! !OUTPUT PARAMETERS:

  integer,              intent(out) :: stat   ! return error code
  type(dynq), optional, intent(out) :: xp     ! final   perturbation from ADM integration

! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  27Sep2007 Todling - Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

  character(len=*), parameter :: myname_ = 'run_'

  integer, parameter :: ndt = 2   ! wired to evolve two dt at a time (for now)
  integer  ierr
  integer  nymdloc, nhmsloc, i
  integer  mynymd, mynhms
  logical  dostep

  stat = 0

! If so, return perturbation in user array
! ----------------------------------------
  if ( present(yp) ) then
       call prognostics_q_dup ( yp, ypert )
  endif

! Apply the ADM operator M^T to ypert adjoint to create xpert adjoint
! -------------------------------------------------------------------
  if ( oneshot_ ) then
       call qmodel_tad ( ypert, xpert )
  else
       call initial_tad
       call prognostics_q_dup ( ypert, xpert )
       mynymd=nymdf; mynhms=nhmsf
       dostep = ( mynymd==nymdi .and. mynhms==nhmsi )
       do while ( .not. dostep )
          call qmodel_tad ( xpert, nymdi=mynymd, nhmsi=mynhms, ntsteps=ndt )
          call tick ( mynymd, mynhms, -ndt*pdt )
          if ( freq_ > 0 ) then
             if ( mod(mynhms,freq_) == 0 ) then
                  call wrap1_ ( ierr, xp=xpert, nymdw=mynymd, nhmsw=mynhms, ofreq=freq_ )
             endif
          endif
          dostep = ( mynymd==nymdi .and. mynhms==nhmsi )
       enddo
       call final_tad
  endif

! If so, return perturbation in user array
! ----------------------------------------
  if ( present(xp) ) then
       call prognostics_q_dup ( xpert, xp )
  endif


  end subroutine run_ 
!EOC

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: wrap1_ --- collects results from ADM integration
!
! !INTERFACE:

  subroutine wrap1_ ( stat, &
                      fname, xp, nymdw, nhmsw, ofreq )

! !USES:

  implicit none

! !INPUT PARAMETERS:

  type(dynq),       optional, intent(in) :: xp     ! perturbation to write out
  character(len=*), optional, intent(in) :: fname  ! output file name (w/o date/time)
  integer,          optional, intent(in) :: nymdw  ! current (wrap) date
  integer,          optional, intent(in) :: nhmsw  ! current (wrap) time
  integer,          optional, intent(in) :: ofreq  ! ouput frequency (HHMMSS)

! !OUTPUT PARAMETERS:

  integer, intent (out) :: stat                  ! return error code

! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  27Sep2007 Todling - Initial code
!
!EOP
!-------------------------------------------------------------------------
!BOC

  character(len=*), parameter :: myname_ = 'wrap1_'

  integer    n, ierr
  integer    nymdw_,nhmsw_,ofreq_
  logical    forceflip
  character(len=255) :: fnpert
  type(dynq)         :: pert_out

  stat = 0

! Convert contiguous array xpert onto model-like state and write out evolved vector
! ---------------------------------------------------------------------------------
  call prognostics_q_initial ( pert_out )      ! allocate perturbation arrays
  if ( present(xp) ) then
       call prognostics_q_dup ( xp, pert_out )
  else
       call prognostics_q_dup ( xpert, pert_out )
  endif
  if ( present(nymdw) .and. present(nhmsw) ) then
       nymdw_ = nymdw
       nhmsw_ = nhmsw
  else
       nymdw_ = nymdi
       nhmsw_ = nhmsi
  endif
  forceflip = .false.
  if (g5pert)  forceflip=.true.
  if (present(ofreq)) then
      ofreq_ = ofreq
  else
      ofreq_ = fvpsasdt
  endif

!_RT call SetPoles( pert_out, ierr )       ! impose proper conditions on both poles

  if(nymd/=nymdi .and. nhms/=nhmsi)then
     stat = 91
     call mpout_log(myname_,' incorrect final date/time')
     return
  endif

  if (ofreq_>0) then
     if ( present(fname) ) then
          fnpert = trim(fname)
          call PutPert ( job, nymdw_, nhmsw_, pert_out, ofreq_, nstep, fnpert,  &
                         forceflip=forceflip, RCfile=RCfile )
     else
          write(fnpert,'(3a)') trim(job), '.', trim(fsensname)
          call PutPert ( job, nymdw_, nhmsw_, pert_out, ofreq_, nstep, fnpert, &
                         nymd_b=nymdf, nhms_b=nhmsf, nymd_e=nymdw_, nhms_e=nhmsw_, &
                         forceflip=forceflip, RCfile=RCfile )
     endif
  endif

  call prognostics_q_final ( pert_out )       ! deallocate perturbation arrays

  end subroutine wrap1_ 
!EOC

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: clean_ --- clean up
!
! !INTERFACE:

  subroutine clean_ ( stat )

! !USES:

  implicit none

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:

  integer, intent (out) :: stat

! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  27Sep2007 Todling - Initial code
!
!EOP
!-------------------------------------------------------------------------
!BOC

  character(len=*), parameter :: myname_ = 'clean_'

  integer    i, ierr

  stat = 0

! Deallocate individual perturbation variables
! --------------------------------------------
  call prognostics_q_final ( ypert )
  call prognostics_q_final ( xpert )

! Release trajectory
! ------------------
  call getstate_clean
  call physdrv1_get_clean

  call stepon_finalize ( prog )

  end subroutine clean_ 
!EOC

  end module m_gQadm
