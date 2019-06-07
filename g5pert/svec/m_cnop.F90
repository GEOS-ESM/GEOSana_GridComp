!-------------------------------------------------------------------------
!    NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_cnop --- Various operators related to CNOPs calculation
!
! !INTERFACE:
                                                                                                                         
  module m_cnop
                                                                                                                         
! !USES:
                                                                                                                         

  use precision

  use stepon,  only : stepon_set
  use control, only : mod2cont

  use m_sVecDef, only : svecnormI  ! Initial state norm
  use m_sVecDef, only : svecnormF  ! Final state norm
  use m_sVecDef, only : cnop_sigma ! Magnitude of initial perturbation
  use m_sVecDef, only : spg_verb   ! controls verbose
  use m_sVecDef, only : cnop_tol   ! tolerance for convergence criterium

  use control, only : ibeg_delp, iend_delp
  use control, only : ibeg_pt, iend_pt

  use m_svnorms, only : norm_svecBC
  use m_svnorms, only : norm_svecD
  use m_svnorms, only : norm_zvecsize
  use m_svprojs, only : proj_svec
                                                                                                                         
  use m_mpif90, only : comm => MP_comm_world

  use m_die, only : die
  use m_stdio, only : stdout
  use m_stdio, only : stderr

  use m_random

  use timingModule
                                                                                                                         
  implicit NONE

#include "lapack.H"
                                                                                                                         
! !PUBLIC MEMBER FUNCTIONS:
                                                                                                                         
  PRIVATE
  
  PUBLIC cnop_init
  PUBLIC cnop_set
  PUBLIC cnop_evalf
  PUBLIC cnop_evalg
  PUBLIC cnop_proj
  PUBLIC cnop_clean

  PUBLIC cnop_fcnt
                   
  interface cnop_init   ; module procedure init_; end interface	
  interface cnop_set    ; module procedure set_; end interface
  interface cnop_evalf  ; module procedure evalf_; end interface
  interface cnop_evalg  ; module procedure evalg_; end interface
  interface cnop_proj   ; module procedure proj_; end interface
  interface cnop_clean  ; module procedure clean_; end interface

  logical,  save                            :: done_init = .false. ! sets initialz. step                                                                                                                         
  logical,  save                            :: done_intg = .false. ! controls model integration
  
  integer,  save                            :: nn    ! dim of full state vector
  real,     save                            :: ynrm0 ! norm of initial perturbation
  
  real(r8), save, allocatable, dimension(:) :: xx    ! initial unperturbed state
  real(r8), save, allocatable, dimension(:) :: yfu   ! final state of unperturbed integration
  real(r8), save, allocatable, dimension(:) :: yfp   ! final state of   perturbed integrations 
  real(r8), save, allocatable, dimension(:) :: ydf   ! diff of unpert and pert integrations 

                                                                                                                         
!
! !DESCRIPTION: Operations to calculate Conditional Nonlinear Optimal
!               Perturbations (CNOPs).
!
! !REVISION HISTORY:
!
!  14Aug2006  Todling   Initial code.
!
!EOP
!-------------------------------------------------------------------------
                                                                                                                         
  character(len=*), parameter :: myname = 'm_cnop'

  logical, parameter :: cnop_debug = .true.
  logical, parameter :: proj_debug = .true.
  integer :: seed = 0

  integer, save :: cnop_fcnt

  integer :: ROOT = 0  ! this is here only for now

  real, parameter :: mytol  = 1.e-10
  real, parameter :: mytiny = 1.e-30
  
  CONTAINS

!-------------------------------------------------------------------------
!    NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_ --- Allocate arrays memory
!
! !INTERFACE:
                                                                                                                         
subroutine init_ ( n, rc )
  implicit none
                                                                                                                         
! !INPUT PARAMETERS:

   integer n, rc                                                                                                                    
! !DESCRIPTION:  Allocate state vector
!
! !REVISION HISTORY:
!
!   14Aug2006  Todling    Initial code.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname_ = myname//'::init_'
  
  integer ierr

  nn = n
  rc = 0
  if (.not.done_init) then
      allocate ( xx(n), yfu(n), yfp(n), ydf(n), stat=ierr )
        if(ierr/=0) then
          rc = 1
          return
        endif
      done_init = .true.
      if(cnop_debug) write(stdout,'(2a)') trim(myname_), ': Allocation done'
  endif
  if(spg_verb) write(stdout,'(2a)') trim(myname_), ': complete.'
    
  end subroutine init_

!-------------------------------------------------------------------------
!    NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: set_ --- Initialize perturbation for CNOP calculation
!
! !INTERFACE:
                    
  subroutine set_( comm, ROOT, n, x, rc, xpert )

! !USES:

  implicit none
  
! !INPUT PARAMETERS:  

  integer, intent(in) :: comm
  integer, intent(in) :: ROOT
  integer, intent(in) :: n
  
! !INPUT/OUTPUT PARAMETERS: 
  
  real(r8), optional :: xpert(n)
  real(r8) x(n)
  
! !OUTPUT PARAMETERS:

  integer, intent(out) :: rc 
  
! !DESCRIPTION:  Initialize perturbation for CNOP calculation.
!
! !REVISION HISTORY:
!
!   14Aug2006  Todling    Initial code.
!
!EOP
!-------------------------------------------------------------------------  
  character(len=*), parameter :: myname_ = myname//'::set_'

  integer  ierr
  real(r8) znorm
  real(r8) _SDOT

  rc = 0

  if ( n/= nn ) then
     rc = 1
     write(stderr,'(2a)') trim(myname_), ': Error dim(x)' 
     return
  endif

! Set fields and parameter for evalf
! ----------------------------------
  cnop_fcnt = 0
  xx = x

! If so, set initial perturbation to random field
! -----------------------------------------------
  if (present(xpert) ) then

      call zufalli ( seed )       ! initialize random number generator
      call normalen(n,xpert)
!_RT  call norm_svecBC ( comm, ROOT, svecnormI, n, xx, xpert, ynrm0, 'SQ', ierr )  ! whether L2 or not
      if ( cnop_debug ) then
           znorm = _SDOT ( n, xpert, 1, xpert, 1 )
           write(stdout,'(1x,2a,1p,e13.6)') myname_, ': norm2(direct) = ', ynrm0
           write(stdout,'(1x,2a,1p,e13.6)') myname_, ': norm2(sdot)   = ', znorm
           if( abs(znorm-ynrm0) > mytol ) call die(myname_,'unmatched norm',ierr)
      endif

  endif
  if(spg_verb) write(stdout,'(2a)') trim(myname_), ': complete.'
  
  end subroutine set_
    
!-------------------------------------------------------------------------
!    NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: evalf_ --- Evaluate functional for CNOP problem
!
! !INTERFACE:

  subroutine evalf_ ( n, xpert, rc, f, fnlpert )

! !USES:

  implicit none
  
! !INPUT PARAMETERS:

  integer, intent(in) ::  n          ! dimension of control vector
  real(r8) :: xpert(n)               ! input (control) perturbation vector

! !OUTPUT PARAMETERS:

  integer, intent(out) :: rc          ! return error code
  real(r8), optional   :: f           ! function is evaluated when present;
                                      !  otherwise only calculates kernel
  real(r8), optional   :: fnlpert(:)  ! final non-linearly evolved perturbation
 
! !DESCRIPTION:
!
!     This subroutine computes the objective function.
!
!     On Entry:
!
!     n     integer,
!           size of the problem,
!
!     xpert real perturbation xpert(n),
!           point at which the function will be evaluated.
!
!     On Return
!
!     f     real,
!           scalar function value at xpert,
!
!     rc    integer,
!           return code (termination parameter):
!           0= the function was succesfully evaluated,
!           1= some error occurs in the function evaluation.
!
!     Functionality:
!       call cnop_init
!       call cnop_evalf
!       call cnop_clean
!
! !REMARKS:
!
!        This function allows only for an approximate calculation 
!        of CNOPs. In actuality this function should be invoking 
!        the full GCM (the same one that creates the trajectory
!        for the adjoint). As implemented here, this functional
!        is simply a stub for a more general implementation to be
!        carried on in the future.
!
! !REVISION HISTORY:
!
!   14Aug2006  Todling    Initial code.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname_ = myname//'::evalf_'

!  local variables
!  ---------------
   integer i, ierr
   real(r8)  _SDOT
   real(r8)  znorm, xnorm, myfunc
   real(r8), allocatable :: x(:)
   real(r8)              :: ynorm   ! proper norm of diff between pert
                                    ! and unpert integrations (final time)
 
   rc = 0

   if ( n/=nn ) call die(myname_,' input dim(n) inconsistent w/ nn ')

   allocate ( x(n), stat=ierr )
      if(ierr/=0) call die(myname_,'Alloc()',ierr)

    if ( cnop_debug ) then
         znorm = _SDOT ( n, xpert, 1, xpert, 1 )
         write(stdout,'(1x,2a,1p,e13.6)') myname_, ': norm2(xpert) = ', znorm
         if(znorm<mytiny) call die('evalf_: input xpert is negligible')
    endif

!  --------------------------------------
!  Integrate model from initial condition
!  -------------------------------------- 
     
     if (.not.done_intg) then
         call timing_on ( 'Model' )    
         call func ( n, xx, n, yfu )
         call timing_off ( 'Model' )
         done_intg = .true.
         if(spg_verb) write(stdout,'(2a)') trim(myname_), ': complete unperturbed integration.'
     endif
     
!  ------------------------------------------------   
!  Integrate model from perturbed initial condition
!  ------------------------------------------------      
     
          x = xx + xpert       ! add perturbation to initial condition

     yfp = 0.0 
     call timing_on ( 'Model' )  
     call func ( n, x, n, yfp )
     call timing_off ( 'Model' )      
      
          ydf = yfp - yfu      ! diff between perturbed and unperturbed integrations     

     if(present(fnlpert)) fnlpert = ydf  ! returned non-linearly evolved perturbation
 
!    Calculate norm of integrated diff but do not scale different vector
!    -------------------------------------------------------------------
!_RT call proj_svec ( comm, ROOT, n, ydf, ierr )   ! Apply local projection operator
        if(ierr/=0) call die(myname_,'Error from proj_svec',ierr)

!_RT call norm_svecBC ( comm, ROOT, svecnormF, n, xx, ydf, ynorm, 'NA', ierr )
!_RT    if(ierr/=0) call die(myname_,'Error from norm_svecBC',ierr)  

!_RT call proj_svec ( comm, ROOT, n, ydf, ierr )   ! Apply local projection operator
        if(ierr/=0) call die(myname_,'Error from proj_svec',ierr)

        if ( cnop_debug ) then
             xnorm = _SDOT ( n, ydf, 1, ydf, 1 )
             write(stdout,'(1x,2a,1p,e13.6)') myname_, ': norm2(ydf) = ', xnorm
             if(xnorm<mytiny) call die('evalf_: ydf is negligible')
        endif

     if(spg_verb) write(stdout,'(2a)') trim(myname_), ': complete perturbed integration.'
     
     myfunc = 1.0 / ynorm
     if (present(f)) then
         f = myfunc
         cnop_fcnt = cnop_fcnt + 1 
     endif

     if (cnop_debug) then
         write(stdout,'(1x,2a,1p,e13.6)') trim(myname_), ': J = ', myfunc
         if(ynorm<mytiny) call die('evalf_: y is negligible')
     endif

     deallocate ( x, stat=ierr )
      
     end subroutine evalf_

!-------------------------------------------------------------------------
!    NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: evalg_ --- Calculate gradient of cost function
!
! !INTERFACE:

     subroutine evalg_ ( n, xpert, g, rc )

! !USES:

     implicit none

! !INPUT PARAMETERS:

     integer, intent(in) :: n           ! dimension of control vector
     real(r8)            :: xpert(n)    ! input (control) perturbation vector

! !OUTPUT PARAMETERS:

     integer, intent(out) :: rc         ! return error code
     real(r8)             :: g(n)       ! gradient of cost function

! !DESCRIPTION:
!      
!     This subroutine computes the gradient of the 
!     objective function.
!
!     On Entry:
!
!     n     integer,
!           size of the problem,
!
!     xpert real xpert(n),
!           point at which the gradient will be evaluated.
!
!     On Return
!
!     g     real g(n),
!           gradient vector at xpert,
!
!     inform integer,
!           termination parameter:
!           0= the gradient was succesfully evaluated,
!           1= some error occurs in the gradient evaluation.
!
! !REVISION HISTORY:
!
!   14Aug2006  Todling    Initial code.
!
!EOP
!-------------------------------------------------------------------------

     character(len=*), parameter :: myname_ = myname//'::evalg_'

!    Local variables
!    ---------------
     integer i, ierr
     real(r8), allocatable :: x(:)
     real(r8) f, ynorminv, dumn1, dumn2
     real(r8) _SDOT

     rc = 0

     if ( n/=nn ) call die(myname_,' input dim(n) inconsistent w/ nn ')

!    Calculate diff between perturbed and unperturbed integrations
!      NOTE: This call calculates yfp and ydf needed for the 
!            following ADM call; these are "saved" vectors.
!    -------------------------------------------------------------
     call evalf_ ( n, xpert, ierr, f=ynorminv )
       if(ierr/=0)  call die(myname_,'Failed in calculating functional',ierr)
       if(spg_verb) write(stdout,'(2a)') trim(myname_), ': complete eval(func).'
     
!    Apply the operator S to calculate g
!    -----------------------------------
!_RT call proj_svec ( comm, ROOT, n, ydf, ierr )   ! Apply local projection operator
        if(ierr/=0) call die(myname_,'Error from proj_svec',ierr)

!_RT call norm_svecBC ( comm, ROOT, svecnormF, n, xx, ydf, dumn1, 'BC', ierr )
        if(ierr/=0) call die(myname_,'Error from norm_svecBC',ierr)

!_RT call proj_svec ( comm, ROOT, n, ydf, ierr )   ! Apply local projection operator
        if(ierr/=0) call die(myname_,'Error from proj_svec',ierr)

        if ( cnop_debug ) then
             if(dumn1<mytiny) call die('evalf_: ydf is negligible')
        endif


     allocate ( x(n), stat=ierr )
        if(ierr/=0) call die(myname_,'Alloc()',ierr)

     g=0.0
     call timing_on  ( 'func_AD' )
!_RT call func_ad ( n, x, g, n, yfp, ydf )
     call timing_off ( 'func_AD' )
     if(spg_verb) write(stdout,'(2a)')  trim(myname_), ': complete ADM integration.'

     deallocate ( x, stat=ierr )
        if(ierr/=0) call die(myname_,'Dealloc()',ierr)

!    Properly scale gradient
!    -----------------------
     g =  - 2.0 * g * ynorminv**2

        if ( cnop_debug ) then
             dumn1 = _SDOT ( n, g, 1, xpert, 1 )
             dumn2 = _SDOT ( n, g, 1,     g, 1 )
             write(stdout,'(1x,2a,1p,e13.6)') myname_, ': final pert norm2(y) = ', 1.0/ynorminv
             write(stdout,'(1x,2a,1p,e13.6)') myname_, ':            norm2(g) = ', dumn2
             write(stdout,'(1x,2a,1p,e13.6)') myname_, ':       norm(g,xpert) = ', dumn1
        endif

     if(spg_verb) write(stdout,'(2a)') trim(myname_), ': complete.'

     end subroutine evalg_

!-------------------------------------------------------------------------
!    NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: proj_ --- Project arbitrary point on to feasible set
!
! !INTERFACE:

     subroutine proj_ ( n, xpert, rc )

! !USES:

     implicit none
     
! !INPUT/OUTPUT PARAMETERS:

     integer, intent(in) :: n
     real(r8)            :: xpert(n)
 
! !OUTPUT PARAMETERS:

     integer, intent(out) :: rc

! !DESCRIPTION:
!
!     This subroutine computes the projection of an arbitrary
!     point onto the feasible set. 
!
!     On Entry:
!
!     n     integer,
!           size of the problem,
!
!     xpert real xpert(n),
!           point that will be projected.
!
!     On Return
!
!     xpert real xpert(n),
!           projected point,
!
!     rc    integer,
!           termination parameter:
!           0= the projection was succesfuly done,
!           1= some error occurs in the projection.
!
! !REVISION HISTORY:
!
!   14Aug2006  Todling    Initial code.
!
!EOP
!-------------------------------------------------------------------------
 
     character(len=*), parameter :: myname_ = myname//'::proj_'

!    local variables
!    ---------------
     integer i, ierr
     real(r8)  sig, znorm, xnorm
     real(r8)  _SDOT
     real(r8), allocatable :: xaux (:)
     character(len=10) :: scaled
  
     if(cnop_debug) write(stdout,'(1x,2a)') trim(myname_), ': starting projection ...'

     rc = 0
     
     if ( n/=nn ) then
          rc = 1
          return
     endif

     allocate ( xaux(n), stat=ierr )
        if (ierr/=0) then
	   write(stderr,'(2a,i7)') trim(myname_), ' Alloc(), Error = ', ierr
	   rc = 2
	   return
	endif   

     xaux  = xpert

     if ( proj_debug ) then
          znorm = _SDOT ( n, xpert, 1, xpert, 1 )
          write(stdout,'(a,1p,2(a,e13.6))') myname_, ': on entry norm2 = ', znorm, & 
                                                      ' sigma = ', cnop_sigma
     endif
      
!_RT call proj_svec   ( comm, ROOT, n, xaux, ierr )   ! Apply local projection operator
!_RT call norm_svecBC ( comm, ROOT, svecnormF, n, xx, xaux, znorm, 'SQ', ierr )
!_RT call proj_svec   ( comm, ROOT, n, xaux, ierr )   ! Apply local projection operator

     if ( proj_debug ) then
          write(stdout,'(a,1p,2(a,e13.6))') myname_, ':   after E-norm = ', znorm, & 
                                                     ' sigma = ', cnop_sigma
     endif
      
!    Rescale vector if its norm not within bounds; ||z||_S <= sigma
!    --------------------------------------------------------------
     scaled = 'Unscaled'
     if ( znorm > cnop_sigma ) then

         sig = sqrt(cnop_sigma) / sqrt(znorm)
         call _SCAL ( n, sig, xaux, 1 )
         xpert = xaux

         scaled = 'Scaled'
     elseif ( znorm < mytiny ) then
         rc = 3
         deallocate ( xaux, stat=ierr )
         return
     end if
      
     if ( proj_debug ) then
          znorm = _SDOT ( n, xpert, 1, xpert, 1 )
          write(stdout,'(a,1p,2(a,e13.6),1x,a)') myname_, ': on exit  norm2 = ', znorm, & 
                                                          ' sigma = ', cnop_sigma, scaled

     endif

     deallocate ( xaux, stat=ierr )
        if (ierr/=0) then
	   write(stderr,'(2a,i7)') trim(myname_), ' Dealloc(), Error = ', ierr
	   rc = 2
	   return
	endif 
	      
     if(spg_verb) write(stdout,'(2a)') trim(myname_), ': complete.'

     end subroutine proj_

!-------------------------------------------------------------------------
!    NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: clean_ --- Deallocate arrays memory
!
! !INTERFACE:
                                                                                                                         
subroutine clean_ ( )

  implicit none
                                                                                                                         
! !INPUT PARAMETERS:
                                                                                                                         
! !DESCRIPTION:  Deallocate state vector
!
! !REVISION HISTORY:
!
!   14Aug2006  Todling    Initial code.
!
!EOP
!-------------------------------------------------------------------------
                                                                                                                         
  character(len=*), parameter :: myname_ = myname//'::clean_'
                                                                                                                         
  integer  ier
                                                                                                                         
  deallocate ( xx, yfu, yfp, ydf, stat=ier )
    if(ier/=0) call die (myname_,'Dealloc()')
                                                                                                                         
end subroutine clean_

end module m_cnop

   

