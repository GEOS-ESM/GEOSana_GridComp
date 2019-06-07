!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_gQtlm --- GEOS Tracer Tangent Linear Model Module
!
! !INTERFACE:
!
  module m_gQtlm

! !USES:

  use precision

  use mod_comm, only : gid,mpi_bcst_n_real_ad

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
  use stepon,  only : ks
  use stepon,  only : ak
  use stepon,  only : bk
  use stepon,  only : ts
  use stepon,  only : oro

  use prognostics, only : dyn_prog

  use prognostics_q, only : prognostics_q_initial
  use prognostics_q, only : prognostics_q_final
  use prognostics_q, only : prognostics_q_dup
  use prognostics_q, only : prognostics_q_dotp
  use prognostics_q, only : prognostics_q_rand
  use prognostics_q, only : dynq

  use m_model_tad, only : qmodel_tad
  use m_model_tad, only : initial_tad
  use m_model_tad, only : final_tad

  use m_model_ttl, only : qmodel_ttl

  use m_cnop, only : cnop_init
  use m_cnop, only : cnop_evalf
  use m_cnop, only : cnop_clean

! use m_poles, only : SetPoles    ! To impose proper pole conditions
                                  ! uncomment line: use m_poles

  use m_iostate, only : PutState
  use m_iostate, only : getstate_init
  use m_iostate, only : getstate
  use m_iostate, only : getstate_clean
                                                                                                                             
  use m_trjphys, only : physdrv1_get_init
  use m_trjphys, only : physdrv1_get_all
  use m_trjphys, only : physdrv1_get_clean

  use m_trajmng, only : GetPert
  use m_trajmng, only : PutPert

  use m_mpout,only : mpout_log
  use m_zeit

  implicit none

  PRIVATE

  PUBLIC gQtlm_init
  PUBLIC gQtlm_set
  PUBLIC gQtlm_run
  PUBLIC gQtlm_wrap
  PUBLIC gQtlm_clean
  PUBLIC gQtlm_dotp
 
! !DESCRIPTION: Evolve perturbation with the FV tangent linear model.
!
! !REVISION HISTORY:
!
!  24Sep2005 Todling - Initial code
!  06Mar2009 Todling - .hdf to .nc4
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname    = 'm_gQtlm'
  character(len=*), parameter :: fnpertout = 'qpert.eta'
  character(len=*), parameter :: fnpertin  = 'qpert.eta.nc4'
  character(len=*), parameter :: RCfile    = 'tracer.rc'

  integer, parameter :: ROOT = 0

! Allocatable arrays
! ------------------
  type(dyn_prog)              :: prog   ! individual prognositic  variables
  type(dynq)                  :: xpert  ! individual perturbation variables
  type(dynq)                  :: ypert  ! individual perturbation variables
 
  integer ::   ndim                     ! dimension of perturbation/state vector

  integer ::   nymdi                    ! initial date of integration
  integer ::   nhmsi                    ! initial time of integration
  integer ::   nymdf                    ! final   date of integration
  integer ::   nhmsf                    ! final   time of integration

  logical, save      :: memtraj    = .true.   ! load trajectory to memory
  logical, save      ::  g5pert    = .false.
  logical, save      ::  oneshot_  = .false.
  integer, save      :: freq_      = -1       ! frequency of TLM output (HHMMSS)

  interface gQtlm_init  ; module procedure init_  ; end interface
  interface gQtlm_set   ; module procedure set1_  ; end interface
  interface gQtlm_run   ; module procedure run_   ; end interface
  interface gQtlm_dotp  ; module procedure dotp_  ; end interface
  interface gQtlm_wrap  ; module procedure wrap1_ ; end interface
  interface gQtlm_clean ; module procedure clean_ ; end interface

  CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_ --- initialized atmospheric tangent linear model 
!
! !INTERFACE:

  subroutine init_ ( stat, &
                     n, verbose, oneshot, freq )   ! optionals

! !USES:

  implicit none

! !INPUT PARAMETERS:

   logical, optional, intent(in)  :: verbose
   logical, optional, intent(in)  :: oneshot ! run step-by-step or not
   integer, optional, intent(in)  :: freq    ! freq of output during TLM (HHMMSS)

! !OUTPUT PARAMETERS:
 
   integer,           intent(out) :: stat
   integer, optional, intent(out) :: n

! !DESCRIPTION: Evolve perturbation with the FV tangent linear model.
!
! !REVISION HISTORY:
!
!  24Sep2007 Todling - Initial code
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  character(len=*), parameter :: myname_ = myname//'::init_'

  character(len=255) words
  integer    i, n1, n2, ierr

  ierr = 0; stat =0

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
  call stepon_set ( prog )

! Store timings for output pert file tagging
! ------------------------------------------
  nymdi = nymd ;   nhmsi = nhms
  nymdf = nymde;   nhmsf = nhmse
    write(words,'(i8.8,1x,i6.6)') nymdi, nhmsi
    call mpout_log(myname,': integration starts at: '//words)
    write(words,'(i8.8,1x,i6.6)') nymdf, nhmsf
    call mpout_log(myname,': integration   ends at: '//words)

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
! !ROUTINE: set1_ --- get initial perturbation from file
!
! !INTERFACE:

  subroutine set1_ ( stat, &
                     fname, vectype, xp )         ! optionals

! !USES:

  implicit none

! !INPUT PARAMETERS:

  character(len=*), optional, intent(in) :: fname   ! filename containing input perturbation
  integer,          optional, intent(in) :: vectype ! specifies g4/g5 perturbation

! !OUTPUT PARAMETERS:

  integer, intent(out)                  :: stat   ! return error code
  type(dynq),     optional, intent(out) :: xp     ! perturbation type

! !DESCRIPTION: Get initial perturbation from file.
!
! !REVISION HISTORY:
!
!  24Sep2007 Todling - Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  character(len=*), parameter :: myname_ = myname//'::set1_'

  character(len=255)  fnpert
  integer    nymdp        ! date of perturbation (not used)
  integer    nhmsp        ! time of perturbation (not used)
  integer    n, ierr

  ierr = 0; stat = 0

! Read in perturbation
! --------------------
!_RT  nymdp = 0; nhmsp = 0
  nymdp = nymdi; nhmsp = nhmsi
  if ( present(fname) ) then
       fnpert = trim(fname)
  else
       fnpert = fnpertin
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
  
! If so, return perturbation in user arrays
! -----------------------------------------
  if ( present(xp) ) then
       call prognostics_q_dup ( xpert, xp )  
  endif

  end subroutine set1_
!EOC

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: run_ --- run atmospheric tangent linear model
!
! !INTERFACE:

  subroutine run_( stat, & 
                   xp, yp ) ! optionals

! !USES:

  implicit none

! !INPUT PARAMETERS:

  type(dynq),     optional, intent(in) :: xp         ! initial perturbation

! !OUTPUT PARAMETERS:

  integer,                  intent(out) :: stat      ! return error code
  type(dynq),     optional, intent(out) :: yp        ! final perturbation

! !DESCRIPTION: Evolve perturbation with the FV tangent linear model.
!
! !REVISION HISTORY:
!
!  24Sep2007 Todling - Initial code
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  character(len=*), parameter :: myname_ = myname//'::run_'

  integer, parameter :: ndt = 2   ! wired to evolve two dt at a time (for now)
  integer  ierr
  integer  mynymd, mynhms
  logical  dostep

  ierr = 0; stat = 0

! If so, get input perturbation from user
! ---------------------------------------
  if ( present(xp) ) then
       call prognostics_q_dup ( xp, xpert )
  endif

! Evolve perturbation with TLM
! ----------------------------
  if ( oneshot_ ) then
       call qmodel_ttl ( xpert, ypert )
  else
       call prognostics_q_dup ( xpert, ypert )  
       mynymd=nymdi; mynhms=nhmsi
       dostep = ( mynymd==nymdf .and. mynhms==nhmsf )
       do while ( .not. dostep )
          call qmodel_ttl ( ypert, nymdi=mynymd, nhmsi=mynhms, ntsteps=ndt )
          call tick ( mynymd, mynhms, ndt*pdt )
          if ( freq_ > 0 ) then
             if ( mod(mynhms,freq_) == 0 ) then
                  call wrap1_ ( ierr, yp=ypert, nymdw=mynymd, nhmsw=mynhms, ofreq=freq_ )
             endif
          endif
          dostep = ( mynymd==nymdf .and. mynhms==nhmsf )
       enddo
  endif

! If so, return evolved perturbation in user array
! ------------------------------------------------
  if ( present(yp) ) then
       call prognostics_q_dup ( ypert, yp )
  endif

  end subroutine run_
!EOC

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: dotp_ --- test TLM/ADM: dot product
!
! !INTERFACE:

  subroutine dotp_ ( seed, stat, & 
                     xp, yp ) ! optionals

! !USES:

  implicit none

! !INPUT PARAMETERS:

  integer,                  intent(in) :: seed       ! random number seed
  type(dynq),     optional, intent(in) :: xp         ! initial perturbation

! !OUTPUT PARAMETERS:

  integer,                  intent(out) :: stat      ! return error code
  type(dynq),     optional, intent(out) :: yp        ! final perturbation

! !DESCRIPTION: Evolve perturbation with the FV tangent linear model.
!
! !REVISION HISTORY:
!
!  24Sep2007 Todling - Initial code
!  28Sep2007 Todling - Fixed issue w/ calling ADM more than once
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  character(len=*), parameter :: myname_ = myname//'::dotp_'

  integer, parameter :: many = 2
  integer    i, n, m, imany, ierr
  real       d0, d1, d2, d11
  logical    dummy
  type(dynq) zpert

  ierr = 0; stat = 0; m = ndim

  call zeit_ci ('dotp')

! Initialize auxiliar structure
! -----------------------------
  call prognostics_q_initial ( zpert )

! If so, get input perturbation from user
! ---------------------------------------
  if ( present(xp) ) then
       call prognostics_q_dup ( xp, xpert )  
  endif

  call initial_tad
  do imany = 1, many

!   Set random perturbation
!   -----------------------
    if ( imany==1 ) then
        if(seed>0) call prognostics_q_rand ( xpert, seed=seed )  
        call prognostics_q_dup ( xpert, ypert )
    else
        call prognostics_q_dup ( ypert, xpert )
    endif

!   Position prog in beginning of time window
!   -----------------------------------------
    call getstate ( nymdi, nhmsi, prog, dummy )

!   Run TLM
!   -------
       d0 = prognostics_q_dotp ( ypert, ypert )
       call zeit_ci ('ttl')
!   call qmodel_ttl ( ypert, nymdi=nymdi, nhmsi=nhmsi, ntsteps=2 )
    call qmodel_ttl ( ypert )
       call zeit_co ('ttl')
       d1 = prognostics_q_dotp ( ypert, ypert )

!   Run ADM
!   -------
!   call initial_tad
       call zeit_ci ('tad')
!   call qmodel_tad ( ypert, nymdi=nymdi, nhmsi=nhmsi, ntsteps=2 )
    call qmodel_tad ( zpert, ypert, setup=.false.  )
       call zeit_co ('tad')
!   call final_tad
       d2 = prognostics_q_dotp ( xpert, zpert )
 
!   Report results of dot-product test
!   ----------------------------------
       if (gid==ROOT) then
                     print*
                     print*, '                d0 = ', d0
                     print*, '                d1 = ', d1
                     print*, '                d2 = ', d2
                     print*, '        abs(d1-d2) = ', abs(d1-d2)
           if(d0>0.) print*, 'abs(d1-d2)/abs(d0) = ', abs(d1-d2)/abs(d0)
                     print*
       endif

  enddo ! < many >
  call final_tad

! If so, return evolved perturbation in user array
! ------------------------------------------------
  if ( present(yp) ) then
       call prognostics_q_dup ( ypert, yp )  
  endif

! Initialize auxiliar structure
! -----------------------------
  call prognostics_q_final ( zpert )

  call zeit_co ('dotp')

  end subroutine dotp_
!EOC

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: wrap1_ --- Output evolved perturbation to a file
!
! !INTERFACE:

  subroutine wrap1_ ( stat, &
                      fname, yp, nymdw, nhmsw, ofreq ) ! optionals

  use stepon,  only : imr
  use stepon,  only : jnp
  use stepon,  only : nl
  use stepon,  only : jfirst
  use stepon,  only : jlast

! !USES:

  implicit none

! !INPUT PARAMETERS:

  type(dynq),       optional, intent(in) :: yp     ! perturbation vector to write out
  character(len=*), optional, intent(in) :: fname  ! output file name (w/o date/time)
  integer,          optional, intent(in) :: nymdw  ! current (wrap) date
  integer,          optional, intent(in) :: nhmsw  ! current (wrap) time
  integer,          optional, intent(in) :: ofreq  ! ouput frequency (HHMMSS)

! !OUTPUT PARAMETERS:

  integer,  intent(out) :: stat            ! return error code

! !DESCRIPTION: Output evolved perturbation to a file
!
! !REVISION HISTORY:
!
!  24Sep2007 Todling - Initial code
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  character(len=*), parameter :: myname_ = myname//'::wrap1_'

  integer    i,j,k,n, m, ierr
  integer    nymdw_,nhmsw_,ofreq_
  logical    forceflip
  character(len=255) :: fnpert
  type(dynq) :: pert_out

  stat = 0; m = ndim

  call prognostics_q_initial ( pert_out )

! Write out evolved vector(s)
! ---------------------------
  if ( present(yp) ) then
       call prognostics_q_dup ( yp,    pert_out )  
  else
       call prognostics_q_dup ( ypert, pert_out )  
  endif
  if ( present(nymdw) .and. present(nhmsw) ) then
       nymdw_ = nymdw
       nhmsw_ = nhmsw
  else
       nymdw_ = nymdf
       nhmsw_ = nhmsf
  endif
  forceflip = .false.
  if (g5pert) forceflip=.true.
  if (present(ofreq)) then
      ofreq_ = ofreq
  else
      ofreq_ = fvpsasdt
  endif 

! Diagnose ps from perturbation delp
! ----------------------------------
  if (ofreq_>0) then
     if ( present(fname) ) then
          fnpert = trim(fname)
          call PutPert ( job, nymdw_, nhmsw_, pert_out, ofreq_, nstep, fnpert,  &
                         forceflip=forceflip, RCfile=RCfile )
     else
          write(fnpert,'(3a)')     trim(job), '.', trim(fnpertout)
          call PutPert ( job, nymdw_, nhmsw_, pert_out, ofreq_, nstep, fnpert, &
                         nymd_b=nymdi, nhms_b=nhmsi, nymd_e=nymdw_, nhms_e=nhmsw_, &
                         forceflip=forceflip, RCfile=RCfile )
     endif
  endif

  call prognostics_q_final ( pert_out )

  end subroutine wrap1_
!EOC

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: clean_ --- Clean up workspace
!
! !INTERFACE:

  subroutine clean_ ( stat )

! !USES:

  implicit none

! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:

   integer, intent(out) :: stat 

! !DESCRIPTION: Clean up workspace
!
! !REVISION HISTORY:
!
!  24Sep2007 Todling - Initial code
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  character(len=*), parameter :: myname_ = myname//'::clean_'
  integer i, ierr

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

end module m_gQtlm
