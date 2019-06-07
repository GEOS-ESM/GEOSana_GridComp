!  
! !REVISION HISTORY:
!
!  14May2007 Todling Introduced dyn_prog; global change.
!  17May2007 Todling Largely revampped from original code:
!                    - turned into model
!                    - interfaces 1 and 2
!  31May2007 Todling Add hook to handle g5-perturbation
!  
!
module m_model_ttl

!==============================================
! referencing used modules
!==============================================
use precision
use prognostics_q, only : prognostics_q_zero
use prognostics_q, only : prognostics_q_dup
use prognostics_q, only : prognostics_q_initial
use prognostics_q, only : prognostics_q_final
use prognostics_q, only : dynq
use prognostics, only : prognostics_initial
use prognostics, only : prognostics_final
use prognostics, only : dyn_prog
use stepon, only : stepon_set
use stepon, only : nymd
use stepon, only : nhms
use stepon, only : ninner
use stepon, only : nouter
use stepon, only : nstep
use stepon_ttl, only : stepon_do_ttl
use m_zeit, only : zeit_ci
use m_zeit, only : zeit_co

implicit none

PRIVATE
PUBLIC qmodel_ttl

interface qmodel_ttl; module procedure &
          model_tl1_, &
          model_tl2_
end interface qmodel_ttl

contains
subroutine model_tl1_( pert_tl, nymdi, nhmsi, ntsteps )

implicit none

!==============================================
! declare arguments
!==============================================
type(dynq) :: pert_tl
integer,optional,intent(in)::nymdi   ! initial date
integer,optional,intent(in)::nhmsi   ! initial time
integer,optional,intent(in)::ntsteps ! number of time steps

type(dyn_prog) :: prog
type(dynq)     :: xpert
type(dynq)     :: ypert

logical, save :: setup = .true.
integer i,nymds,nhmss,nouters,nstepsv
logical reset

!----------------------------------------------
! RESET TIME IN TLM
!----------------------------------------------
reset = present(nymdi) .and. present(nhmsi) .and. present(ntsteps)
if ( reset ) then
     nstepsv= ninner+(nouter-1)*ninner
     nymds  = nymd ; nhmss=nhms ; nouters=nouter
     nymd   = nymdi; nhms =nhmsi; nouter =ntsteps
     setup  = .false.
endif

!----------------------------------------------
! RESET GLOBAL TANGENT VARIABLES
!----------------------------------------------
call prognostics_q_initial ( xpert )
call prognostics_q_initial ( ypert )

!----------------------------------------------
! TANGENT LINEAR AND FUNCTION STATEMENTS
!----------------------------------------------

call prognostics_q_dup ( pert_tl, xpert )

call model_tl2_ ( xpert, ypert, setup=setup )

call prognostics_q_zero ( pert_tl )
call prognostics_q_dup  ( ypert, pert_tl )

!----------------------------------------------
! FINALIZE GLOBAL TANGENT VARIABLES
!----------------------------------------------
call prognostics_q_final ( ypert )
call prognostics_q_final ( xpert )

!----------------------------------------------
! RESET TIME BACK TO ORIGINAL IN TLM
!----------------------------------------------
if ( reset ) then
     nstep  = nstepsv
     nymd   = nymds; nhms =nhmss; nouter =nouters
     setup  = .true.
endif

end subroutine model_tl1_

subroutine model_tl2_( xpert, ypert, setup )

implicit none

!==============================================
! declare arguments
!==============================================
type(dynq) :: xpert
type(dynq) :: ypert
logical,optional,intent(in) :: setup

type(dyn_prog) :: prog
logical, save :: setup_  = .true.
integer i

if (present(setup)) then
    setup_ = setup
endif
!----------------------------------------------
! RESET GLOBAL TANGENT VARIABLES
!----------------------------------------------
call prognostics_initial ( prog )

!----------------------------------------------
! TANGENT LINEAR AND FUNCTION STATEMENTS
!----------------------------------------------
if (setup_) then
    call stepon_set ( prog )
endif

   call zeit_ci('stpqtlm')
call stepon_do_ttl  ( prog, xpert, ypert )
   call zeit_co('stpqtlm')

!----------------------------------------------
! FINALIZE GLOBAL TANGENT VARIABLES
!----------------------------------------------
call prognostics_final ( prog )

end subroutine model_tl2_

end module m_model_ttl

