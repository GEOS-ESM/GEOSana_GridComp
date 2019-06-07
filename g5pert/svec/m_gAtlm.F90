!#define _SINGLE_
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_gAtlm --- GEOS Atmospheric Tangent Linear Model Module
!
! !INTERFACE:
!
  module m_gAtlm

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
  use stepon,  only : ptop
  use stepon,  only : ks
  use stepon,  only : ak
  use stepon,  only : bk
  use stepon,  only : ts
  use stepon,  only : oro

  use prognostics, only : prognostics_initial
  use prognostics, only : prognostics_final
  use prognostics, only : prognostics_dup
  use prognostics, only : prognostics_dotp
  use prognostics, only : prognostics_rand
  use prognostics, only : prognostics_axpy
  use prognostics, only : prognostics_scal
  use prognostics, only : dyn_prog

  use m_model_ad, only : amodel_ad
  use m_model_ad, only : initial_ad
  use m_model_ad, only : final_ad

  use m_model_tl, only : amodel_tl
  use m_model_tl, only : initial_tl
  use m_model_tl, only : final_tl

  use m_fdiff,    only : fdiff_init
  use m_fdiff,    only : fdiff_run
  use m_fdiff,    only : fdiff_clean

  use m_cnop, only : cnop_init
  use m_cnop, only : cnop_evalf
  use m_cnop, only : cnop_clean

  use m_poles, only : SetPoles    ! To impose proper pole conditions
                                  ! uncomment line: use m_poles

  use m_delp2ps, only : delp2ps

! use cnvpert_tl, only : g4tog5_tl
! use cnvpert_tl, only : g5tog4_tl

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
  use m_stdio, only : stdout
  use m_zeit

  implicit none

  PRIVATE

  PUBLIC gAtlm_init
  PUBLIC gAtlm_set
  PUBLIC gAtlm_run
  PUBLIC gAtlm_wrap
  PUBLIC gAtlm_clean
  PUBLIC gAtlm_dotp
 
  SAVE

! !DESCRIPTION: Evolve perturbation with the FV tangent linear model.
!
! !REVISION HISTORY:
!
!  07Apr2005 Todling - Initial code
!  20Sep2005 Elena N.- Added a call to impose proper conditions on poles
!                      Commented out now.
!  09May2007 Todling - Modularized
!  11Dec2007 Todling - Default is to now NOT load trajectory to memory
!  06Mar2009 Todling - .hdf to .nc4
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname    = 'm_gAtlm'
  character(len=*), parameter :: fnpertout = 'fvpert.eta'
  character(len=*), parameter :: fnpertin  = 'fvpert.eta.nc4'

  integer, save      :: nvec_tl = 2 
  integer, parameter :: ROOT = 0

! Allocatable arrays
! ------------------
  type(dyn_prog)              :: prog   ! individual prognositic  variables
  type(dyn_prog)              :: xpert  ! individual perturbation variables
  type(dyn_prog)              :: ypert  ! individual perturbation variables
 
  integer ::   ndim                     ! dimension of perturbation/state vector

  integer ::   nymdi                    ! initial date of integration
  integer ::   nhmsi                    ! initial time of integration
  integer ::   nymdf                    ! final   date of integration
  integer ::   nhmsf                    ! final   time of integration

  logical, parameter :: nonlinear_ = .false.  ! get this from rc file/command line
  logical, save      :: memtraj_   = .false.  ! load trajectory to memory
  logical, save      ::  g5pert    = .false.
  logical, save      ::  oneshot_  = .false.
  logical, save      :: differencing_ = .false.
  integer, save      :: freq_      = -1       ! frequency of TLM output (HHMMSS)

  interface gAtlm_init  ; module procedure init_  ; end interface
  interface gAtlm_set   ; module procedure set1_  ; end interface
  interface gAtlm_run   ; module procedure run_   ; end interface
  interface gAtlm_dotp  ; module procedure dotp_  ; end interface
  interface gAtlm_wrap  ; module procedure wrap1_ ; end interface
  interface gAtlm_clean ; module procedure clean_ ; end interface

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
                     n, verbose, oneshot, freq, memtraj, differencing )   ! optionals

! !USES:

  implicit none

! !INPUT PARAMETERS:

   logical, optional, intent(in)  :: verbose
   logical, optional, intent(in)  :: oneshot ! run step-by-step or not
   integer, optional, intent(in)  :: freq    ! freq of output during TLM (HHMMSS)
   logical, optional, intent(in)  :: memtraj ! when .t., loads trajectory to memory
   logical, optional, intent(in)  :: differencing ! when .t., TLM replaced by 2-NLMs

! !OUTPUT PARAMETERS:
 
   integer,           intent(out) :: stat
   integer, optional, intent(out) :: n

! !DESCRIPTION: Evolve perturbation with the FV tangent linear model.
!
! !REVISION HISTORY:
!
!  09May2007 Todling - Initial code
!  28Aug2007 Todling - Added oneshot option
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

  if ( present(memtraj) ) memtraj_ = memtraj

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
  call getstate_init ( nymd, nhms, memtrj=memtraj_, verbose=verbose )
  call getstate      ( nymd, nhms )

! Initialize physics trajectory handle
! -------------------------------------
  call physdrv1_get_init ( nymd, nhms, memphys=memtraj_, verbose=verbose )
  call physdrv1_get_all  ( nymd, nhms )

! Allocate individual perturbation variables
! ------------------------------------------
  call prognostics_initial ( xpert )
  call prognostics_initial ( ypert )

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
  if ( present(differencing) ) then
       differencing_ = differencing
       if (differencing_) then
           call fdiff_init ( )
       endif
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
  type(dyn_prog), optional, intent(out) :: xp     ! perturbation type

! !DESCRIPTION: Get initial perturbation from file.
!
! !REVISION HISTORY:
!
!  09May2007 Todling - Initial code
!  28Aug2007 Todling - Add ability to read in g5 perturbation
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
  nymdp = 0; nhmsp = 0
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
  if ( g5pert ) then
       call GetPert ( fnpert, nymdp, nhmsp, xpert, pick=.false., vectype=vectype, forceflip=.true., stat=ierr )
  else
       call GetPert ( fnpert, nymdp, nhmsp, xpert, pick=.false., stat=ierr )
  endif
    if(ierr/=0)then
           stat = 90
           call mpout_log(myname_,'Error retrieving perturbation') 
           return
    endif
  
! If so, return perturbation in user arrays
! -----------------------------------------
  if ( present(xp) ) then
       call prognostics_dup ( xpert, xp )  
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
                   nonlinear, xp, yp ) ! optionals

! !USES:

  implicit none

! !INPUT PARAMETERS:

  logical,        optional, intent(in) :: nonlinear  ! determines how to evolve pert
  type(dyn_prog), optional, intent(in) :: xp         ! initial perturbation

! !OUTPUT PARAMETERS:

  integer,                  intent(out) :: stat      ! return error code
  type(dyn_prog), optional, intent(out) :: yp        ! final perturbation

! !DESCRIPTION: Evolve perturbation with the FV tangent linear model.
!
! !REVISION HISTORY:
!
!  09May2007 Todling - Initial code
!  28Aug2007 Todling - Add ability to propagate g5 perturbation
!  12Dec2009 Todling - Allow incremental integration to reproduce single-shot
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  character(len=*), parameter :: myname_ = myname//'::run_'

  integer, parameter :: ndt = 2   ! wired to evolve two dt at a time (for now)
  integer    n, m, ierr
  integer    mynymd, mynhms
  integer    bfnymd,bfnhms  ! before final date/time: final time - ndt*pdf
  logical    dostep,which,tlfirst,tllast

  ierr = 0; stat = 0; m = ndim

! If so, get input perturbation from user
! ---------------------------------------
  if ( present(xp) ) then
       call prognostics_dup ( xp, xpert )  
  endif

! Determines whether to evolve pert w/ TLM or two nonlinear integrations
! ----------------------------------------------------------------------
  which = nonlinear_
  if ( present(nonlinear) ) then
       which = nonlinear
  endif
  
  if ( which ) then

!   Evolve perturbation with two non-linear integrations
!   ----------------------------------------------------
!_RTcall cnop_init  ( m, ierr )
!_RTcall cnop_evalf ( m, xpert, ierr, fnlpert=ypert )
!_RTcall cnop_clean ( )

  else

!   Evolve perturbation with TLM
!   ----------------------------
    if ( oneshot_ ) then
            call zeit_ci ('tlm')
         call amodel_tl ( xpert, ypert, g5pert=g5pert )
            call zeit_co ('tlm')
    else ! the follow incremental integration is set to reproduce the oneshot integration
         if(gid==0) print *, 'incremental integration ...'
         call initial_tl
         call prognostics_dup ( xpert, ypert )  
         mynymd=nymdi; mynhms=nhmsi
         dostep = ( mynymd==nymdf .and. mynhms==nhmsf )
         bfnymd = nymdf
         bfnhms = nhmsf
         call tick ( bfnymd, bfnhms, -ndt*pdt )
         tlfirst = .true.
         if ( freq_ > 0 ) then
            tllast  = .true.
         else
            tllast  = .false.
         endif
         do while ( .not. dostep )
            if ( mynymd==bfnymd .and. mynhms==bfnhms ) tllast = .true.
               call zeit_ci ('tlm')
            write(6,'(a,i8.8,1x,i6.6)')'model_tl: tlm nymd,nhms=',mynymd,mynhms
            call amodel_tl ( ypert, nymdi=mynymd, nhmsi=mynhms, ntsteps=ndt, g5pert=g5pert, tlfirst=tlfirst, tllast=tllast )
               call zeit_co ('tlm')
            call tick ( mynymd, mynhms, ndt*pdt )
            if ( freq_ > 0 ) then
               if ( mod(mynhms,freq_) == 0 ) then
                    call wrap1_ ( ierr, yp=ypert, nymdw=mynymd, nhmsw=mynhms, ofreq=freq_ )
               endif
            else
               tlfirst=.false.
            endif
            dostep = ( mynymd==nymdf .and. mynhms==nhmsf )
         enddo
         write(6,'(a,i8.8,1x,i6.6)')'model_tl: tlm nymd,nhms=',mynymd,mynhms
         call final_tl
    endif

  endif

! If so, return evolved perturbation in user array
! ------------------------------------------------
  if ( present(yp) ) then
       call prognostics_dup ( ypert, yp )  
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
                     xp, yp, checkpoint_tl, checkpoint_ad ) ! optionals

! !USES:

  implicit none

! !INPUT PARAMETERS:

  integer,                  intent(in) :: seed          ! random number seed
  type(dyn_prog), optional, intent(in) :: xp            ! initial perturbation
  integer,        optional, intent(in) :: checkpoint_tl ! determines checkpoint level for TLM
  integer,        optional, intent(in) :: checkpoint_ad ! determines checkpoint level for ADM

! !OUTPUT PARAMETERS:

  integer,                  intent(out) :: stat      ! return error code
  type(dyn_prog), optional, intent(out) :: yp        ! final perturbation

! !DESCRIPTION: Evolve perturbation with the FV tangent linear model.
!
! !REVISION HISTORY:
!
!  07Jun2007 Todling - Initial code
!  18Jul2008 Todling - Add checkpoint; add differencing integration
!  12Dec2009 Todling - Allow incremental integration to reproduce single-shot
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  character(len=*), parameter :: myname_ = myname//'::dotp_'

  integer, parameter :: ndt = 2   ! wired to evolve two dt at a time (for now)
  integer    i, n, m, ierr
  integer    mynymd, mynhms, ntk, ntick
  integer    ainymd, ainhms
  integer    bfnymd, bfnhms
  integer    inymd, inhms
  real       d0, d1, d2, d11
  logical    dummy,dostep
  logical    tlfirst,tllast
  logical    adfirst,adlast

  ierr = 0; stat = 0; m = ndim

   call zeit_ci ('dotp')

! Position prog in beginning of time window
! -----------------------------------------
  call getstate ( nymdi, nhmsi, prog, dummy )

! If so, get input perturbation from user
! ---------------------------------------
  if ( present(xp) ) then
       call prognostics_dup ( xp, xpert )  
  endif

! Set random perturbation
! -----------------------
  if ( seed .ge. 0 ) then
      call prognostics_rand ( xpert, seed=seed )  
  endif

! If so, convert perturbations to GEOS-5 type
! -------------------------------------------
!  if ( g5pert ) then
!       call fake_g5pert_ ( prog, xpert )
!  endif


!   Evolve perturbation with TLM
!   ----------------------------
    call prognostics_dup ( xpert, ypert )

!   Run TLM
!   -------
       d0 = prognostics_dotp ( ypert, ypert )
    if ( oneshot_ ) then
             call zeit_ci ('tlm')
       if ( differencing_ ) then
          call fdiff_run ( prog, ypert, checkpoint=checkpoint_tl, g5pert=g5pert )
       else
          call amodel_tl ( ypert, g5pert=g5pert )
       endif
             call zeit_co ('tlm')
    else
       call initial_tl
#ifdef _SINGLE_
          call zeit_ci ('tlm')
       call amodel_tl ( ypert, nymdi=nymdi, nhmsi=nhmsi, ntsteps=ndt, g5pert=g5pert )
          call zeit_co ('tlm')
#else
         if(gid==0) print *, 'incremental integration ...'
         mynymd=nymdi; mynhms=nhmsi
         dostep = ( mynymd==nymdf .and. mynhms==nhmsf )
         bfnymd = nymdf
         bfnhms = nhmsf
         call tick ( bfnymd, bfnhms, -ndt*pdt )
         tlfirst = .true.
         tllast  = .false.
         do while ( .not. dostep )
            if ( mynymd==bfnymd .and. mynhms==bfnhms ) tllast = .true.
            write(6,'(a,i8.8,1x,i6.6)')'model_tl: tlm nymd,nhms=',mynymd,mynhms
               call zeit_ci ('tlm')
            call amodel_tl ( ypert, nymdi=mynymd, nhmsi=mynhms, ntsteps=ndt, g5pert=g5pert, tlfirst=tlfirst, tllast=tllast )
               call zeit_co ('tlm')
            call tick ( mynymd, mynhms, ndt*pdt )
            tlfirst = .false.
            dostep = ( mynymd==nymdf .and. mynhms==nhmsf )
         enddo
         write(6,'(a,i8.8,1x,i6.6)')'model_tl: tlm nymd,nhms=',mynymd,mynhms
#endif
         call final_tl
    endif
       d1 = prognostics_dotp ( ypert, ypert )

!   Run ADM
!   -------
    if ( oneshot_ ) then
          call zeit_ci ('adm')
       call amodel_ad ( ypert, g5pert=g5pert, checkpoint=checkpoint_ad )
          call zeit_co ('adm')
    else
       call initial_ad
#ifdef _SINGLE_
          call zeit_ci ('adm')
       call amodel_ad ( ypert, nymdi=nymdi, nhmsi=nhmsi, ntsteps=ndt, g5pert=g5pert, checkpoint=checkpoint_ad )
          call zeit_co ('adm')
#else
       mynymd=nymdf; mynhms=nhmsf
       call tick ( mynymd, mynhms, -ndt*pdt )
       write(6,'(a,i8.8,1x,i6.6)')'model_ad: tlm nymd,nhms=',mynymd,mynhms
       inymd=nymdi;inhms=nhmsi
       ainymd=inymd;ainhms=inhms
       call tick ( inymd, inhms, -ndt*pdt )
       adlast  = .true.
       adfirst = .false.
       dostep  = ( mynymd==inymd .and. mynhms==inhms )
       do while ( .not. dostep )
          if ( mynymd==ainymd .and. mynhms==ainhms ) adfirst = .true.
             call zeit_ci ('adm')
          call amodel_ad ( ypert, nymdi=mynymd, nhmsi=mynhms, ntsteps=ndt, g5pert=g5pert, checkpoint=checkpoint_ad, adfirst=adfirst, adlast=adlast )
             call zeit_co ('adm')
          call tick ( mynymd, mynhms, -ndt*pdt )
          write(6,'(a,i8.8,1x,i6.6)')'model_ad: tlm nymd,nhms=',mynymd,mynhms
          adlast=.false.
          dostep = ( mynymd==inymd .and. mynhms==inhms )
       enddo
#endif
       call final_ad
    endif
       d2 = prognostics_dotp ( xpert, ypert )
 
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

! If so, return evolved perturbation in user array
! ------------------------------------------------
  if ( present(yp) ) then
       call prognostics_dup ( ypert, yp )  
  endif

  call zeit_co ('dotp')

  end subroutine dotp_
!EOC

  subroutine fake_g5pert_ ( prog, pert )
  use stepon,  only : imr
  use stepon,  only : jnp
  use stepon,  only : nl
  use stepon,  only : jfirst
  use stepon,  only : jlast
  use stepon,  only : kfirst
  use stepon,  only : klast
  use stepon,  only : ng_d
  use stepon,  only : ng_s
  use stepon,  only : coslon
  use stepon,  only : sinlon
  implicit none
  type(dyn_prog) prog
  type(dyn_prog) pert
  integer ims, ime
  integer jms, jme
  integer kms, kme
  real(r8), allocatable :: ua(:,:,:)
  real(r8), allocatable :: va(:,:,:)
! real(r8)              :: ua(imr,jfirst-1:jlast,nl)
! real(r8)              :: va(imr,jfirst:jlast,nl)

! Pretend input Theta is actually Temperature
! -------------------------------------------
  call th2t_tl ( prog%delp, pert%delp, prog%pt, pert%pt )

! Pretend winds are on the A-grid
! -------------------------------
   allocate ( ua(imr,jfirst-1:jlast,nl), va(imr,jfirst:jlast,nl) )
   ua = 0._r8; va = 0._r8
   call d2a3d ( pert%u(:,jfirst:jlast+1,:), pert%v, ua(:,jfirst:jlast,:), va,  &
                imr, jnp, nl, jfirst, jlast, ng_d, ng_s, coslon, sinlon )
   pert%u = 0._r8
   pert%v = 0._r8
   pert%u(:,jfirst-1:jlast,:) = ua(:,jfirst-1:jlast,:)
   pert%v(:,jfirst:jlast,:)   = va(:,jfirst  :jlast,:)
!  call a2d3d  ( ua, va, pert%u, pert%v, &
!                imr, jnp, nl, jfirst, jlast, ng_d, ng_s, coslon, sinlon )
   deallocate(ua,va)

  end subroutine fake_g5pert_

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

  type(dyn_prog),   optional, intent(in) :: yp     ! perturbation vector to write out
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
!  09May2007 Todling - Initial code
!  28Aug2007 Todling - Add ability to write out g5 perturbation
!  19Jan2010 Todling - Flip fields on way out; appy pole correction(!)
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  character(len=*), parameter :: myname_ = myname//'::wrap1_'

  integer    i,j,k,n, m, vectype, ierr
  integer    nymdw_,nhmsw_,ofreq_
  logical    forceflip_
  character(len=255) :: fnpert
  type(dyn_prog) :: pert_out
  real(r8), allocatable :: ps(:,:)

  stat = 0; m = ndim

  call prognostics_initial ( pert_out )

! Write out evolved vector(s)
! ---------------------------
  if ( present(yp) ) then
       call prognostics_dup ( yp,    pert_out )  
  else
       call prognostics_dup ( ypert, pert_out )  
  endif
  if ( present(nymdw) .and. present(nhmsw) ) then
       nymdw_ = nymdw
       nhmsw_ = nhmsw
  else
       nymdw_ = nymdf
       nhmsw_ = nhmsf
  endif
  vectype = 4
  forceflip_=.false.
  if (g5pert)  then
      vectype=5
      forceflip_=.true.
  endif
  if (present(ofreq)) then
      ofreq_ = ofreq
  else
      ofreq_ = fvpsasdt
  endif 

! Diagnose ps from perturbation delp
! ----------------------------------
  if (ofreq_>0) then
     allocate ( ps(imr,jfirst:jlast) )
     call delp2ps ( pert_out%delp, ps )
     call SetPoles( pert_out, ierr )     !_RT: I don't like this but will leave it here for now
     if ( present(fname) ) then
          fnpert = trim(fname)
          call PutPert ( job, nymdw_, nhmsw_, pert_out, ofreq_, nstep, &
                         ak, bk, Ts, oro, ps, fnpert, vectype=vectype, forceflip=forceflip_ )
     else
          write(fnpert,'(3a)')     trim(job), '.', trim(fnpertout)
          call PutPert ( job, nymdw_, nhmsw_, pert_out, ofreq_, nstep, &
                         ak, bk, Ts, oro, ps, fnpert, vectype=vectype, forceflip=forceflip_,&
                         nymd_b=nymdi, nhms_b=nhmsi, nymd_e=nymdw_, nhms_e=nhmsw_ )
     endif
     deallocate(ps)
  endif

  call prognostics_final ( pert_out )

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
!  09May2007 Todling - Initial code
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  character(len=*), parameter :: myname_ = myname//'::clean_'
  integer i, ierr

  stat = 0

  call fdiff_clean ()

! Deallocate individual perturbation variables
! --------------------------------------------
  call prognostics_final ( ypert )
  call prognostics_final ( xpert )

! Release trajectory
! ------------------
  call getstate_clean
  call physdrv1_get_clean

  call stepon_finalize ( prog )

  end subroutine clean_
!EOC

end module m_gAtlm
