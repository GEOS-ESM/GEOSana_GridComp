!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_fdiff --- GEOS Atmospheric Finite Difference Pert Model
!
! !INTERFACE:
!
  module m_fdiff

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

  use cnvpert_tl, only : g4tog5_tl
  use cnvpert_tl, only : g5tog4_tl

  use m_mpout,only : mpout_log
  use m_stdio, only : stdout
  use m_zeit

  implicit none

  PRIVATE

  PUBLIC fdiff_init
  PUBLIC fdiff_run
  PUBLIC fdiff_clean
 
  SAVE

! !DESCRIPTION: Evolve perturbation with the FV tangent linear model.
!
! !REVISION HISTORY:
!
!  19Jul2008 Todling - Initial code
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname    = 'm_fdiff'

  integer, parameter :: ROOT = 0

! Allocatable arrays
! ------------------
  type(dyn_prog) :: uprog  ! store result of 1st NLM run for checkpointing
 
  integer, save  :: checkpoint_  = 0        ! default: don't checkpoint
  logical, save  :: checkpointed = .false.  ! status of checkpoint
  logical, save  :: initialized  = .false.  ! status of internal state
  logical, save  :: g5pert_ = .true.        ! expects to get GEOS-5 perturbation

  interface fdiff_init ; module procedure init_ ; end interface
  interface fdiff_run  ; module procedure run_  ; end interface
  interface fdiff_clean; module procedure clean_; end interface


CONTAINS


   subroutine init_ ( )
   implicit none
   end subroutine init_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: run_ --- Finite difference TLM
!
! !INTERFACE:

   subroutine run_ ( prog, xp, checkpoint, g5pert )
 
! !USES:

   implicit none

! !INPUT PARAMETERS:

   type(dyn_prog),    intent(in)    :: prog
   type(dyn_prog),    intent(inout) :: xp
   integer, optional, intent(in)    :: checkpoint
   logical, optional, intent(in)    :: g5pert

! !DESCRIPTION: Evolve perturbation using finite-differencing of two
!               consecutive non-linear runs.
!
! !REMARKS: This routines represents a linear approximation to the TLM
!           operator - it is not a tangent linear integration (it's only
!           tangent linear when running for a single a time step).
!           A tangent linear correspondent would have to save every
!           step of the prog evolution - it can be done, but it's 
!           memory intense.
!
! !REVISION HISTORY:
!
!  19Jul2008 Todling - Initial code
!
!EOP
!-------------------------------------------------------------------------

!  Local variables
!  ---------------
   character(len=*), parameter :: myname_ = myname//"::evalf_"
   logical, parameter :: verbose_ = .true.
   type(dyn_prog) :: pprog
 
   if (present(checkpoint)) then
       checkpoint_ = checkpoint
   endif
   if (present(g5pert)) then
       g5pert_ = g5pert
   endif

!  Initialize unperturbed state vector
!  -----------------------------------
     if(.not.initialized) then
        call prognostics_initial ( uprog )
        initialized=.true.
     endif
     call prognostics_initial ( pprog )

!  If needed, convert perturbation
!  -------------------------------
     if (g5pert_) then
         call g5tog4_tl ( prog, xp )
     endif

!  If originally integrated state not in memory ...
!  ------------------------------------------------
     if ( .not. checkpointed ) then

!    Initialize new state vector
!    ---------------------------
       call stepon_set ( uprog )

!    --------------------------------------
!    Integrate model from initial condition
!    --------------------------------------

       call zeit_ci ( 'model' )
       call func ( uprog )
       call zeit_co ( 'model' )
           if(verbose_.and.gid==ROOT) write(stdout,'(2a,i8.8,1x,i6.6)') trim(myname_), ': complete unperturbed integration: ', nymd, nhms

     endif

!  Initialize new state vector
!  ---------------------------
     call stepon_set ( pprog )

!  Initialize perturbed state vector
!  ---------------------------------
     call prognostics_axpy ( 1.0, xp, pprog )

!  ------------------------------------------------
!  Integrate model from perturbed initial condition
!  ------------------------------------------------

     call zeit_ci ( 'model' )
     call func ( pprog )
     call zeit_co ( 'model' )
         if(verbose_.and.gid==ROOT) write(stdout,'(2a,i8.8,1x,i6.6)') trim(myname_), ': complete   perturbed integration: ', nymd, nhms
         call flush(6)

     call prognostics_dup ( pprog, xp )
     call prognostics_axpy ( -1.0, uprog, xp )

!  If needed, convert perturbation back
!  ------------------------------------
     if (g5pert_) then
         call g4tog5_tl ( prog, xp )
     endif

!  Clean up
!  --------
     call prognostics_final ( pprog )
     if (checkpoint_==0) then
       call prognostics_final ( uprog )
       initialized=.false.
     else
       checkpointed = .true.
     endif

     end subroutine run_

   subroutine clean_ ()
   implicit none
   if(checkpoint_>0.and.initialized) then
      call prognostics_final ( uprog )
      initialized=.false.
   endif
   end subroutine clean_

end module m_fdiff
