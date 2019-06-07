!  14May2007  Todling   Introduced dyn_prog; global change.
module dependentq
  use precision
  use prognostics
  use control, only : ibeg_q   , iend_q
#ifdef SPMD
  use stepon, only: ng_d, ng_s
  use mod_comm, only: mp_gather4d
#endif

 real(r8), pointer :: q(:,:,:,:)

contains

subroutine dependent_number( m )

  implicit none
  integer :: m

  m = 0
  m = m + imr * jnp * nl

end subroutine dependent_number

subroutine model2dependent( m, y, prog )
! Project the final model state onto a dependent vector.
! For singular vector detection the final norm should be defined
! here be a (matrix free) matrix vector multiplication


  implicit none
  integer  :: m
  real(r8) :: y(m)
  type(dyn_prog), TARGET :: prog

  q => prog%q
  !----------------------------------------------------------------------- 

#ifdef SPMD
  call mp_gather4d( q(:,:,:,1), y(ibeg_q), &
                    imr, jnp, nl, 1, jfirst, jlast, &
                    kfirst, klast, ng_d, ng_d, 0 )
#else
  y(ibeg_q :iend_q )   = reshape(q (:,jfirst:jlast,:,1),(/1+iend_q -ibeg_q /))
#endif

end subroutine model2dependent

end module dependentq
