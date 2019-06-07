!  14May2007  Todling   Introduced dyn_prog; global change.
!  13Aug2007  Todling   Extended to handle nc>1
!  05Sep2007  Todling   Added handling the general nc number of q
module dependent

  use precision
  use prognostics, only : dyn_prog
  use prognostics, only : imr,jnp,nl,nc
  use prognostics, only : jfirst,jlast
  use prognostics, only : kfirst,klast

  use control, only : ibeg_u   , iend_u
  use control, only : ibeg_v   , iend_v
  use control, only : ibeg_delp, iend_delp
  use control, only : ibeg_pt  , iend_pt
  use control, only : ibeg_q   , iend_q
  use control, only : ibeg_tr  , iend_tr

#ifdef SPMD
  use stepon, only: ng_d, ng_s
  use mod_comm, only: mp_gather4d, mp_bcst_n_real
#endif

PRIVATE

  PUBLIC dependent_number
  PUBLIC model2dependent 

  real(r8), pointer ::    u(:,:,:)   ! zonal wind on D-grid
  real(r8), pointer ::    v(:,:,:)   ! meridional wind
  real(r8), pointer ::   pt(:,:,:)   ! virtual potential temperature
  real(r8), pointer :: delp(:,:,:)   ! pressure thickness (pascal)
  real(r8), pointer ::    q(:,:,:,:) ! specific humidity & tracer mixing ratios
  real(r8), pointer ::   ps(:,:)     ! surface pressure (pascal)

contains

subroutine dependent_number( m )

  implicit none
  integer :: m, m_nc

  m = 0
  m = m + imr * jnp * nl
  m = m + imr * jnp * nl
  m = m + imr * jnp * nl
  m = m + imr * jnp * nl
  m = m + imr * jnp * nl
  do m_nc = 2,nc
     m = m + imr * jnp * nl
  enddo

end subroutine dependent_number

subroutine model2dependent( m, y, prog )
! Project the final model state onto a dependent vector.
! For singular vector detection the final norm should be defined
! here be a (matrix free) matrix vector multiplication


  implicit none
  integer  :: m, m_nc
  real(r8) :: y(m)
  type(dyn_prog), TARGET :: prog

! Set pointers
! ------------
       u   => prog%u
       v   => prog%v
      pt   => prog%pt
      delp => prog%delp
      q    => prog%q

  !----------------------------------------------------------------------- 

#ifdef SPMD
  call mp_gather4d( u   , y(ibeg_u)   , imr, jnp, nl, 1, jfirst, jlast, &
                    kfirst, klast, ng_d, ng_s, 0 )
  call mp_gather4d( v   , y(ibeg_v)   , imr, jnp, nl, 1, jfirst, jlast, &
                    kfirst, klast, ng_s, ng_d, 0 )
  call mp_gather4d( pt  , y(ibeg_pt)  , imr, jnp, nl, 1, jfirst, jlast, &
                    kfirst, klast, ng_d, ng_d, 0 )
  call mp_gather4d( q(:,:,:,1), y(ibeg_q), &
                    imr, jnp, nl, 1, jfirst, jlast, &
                    kfirst, klast, ng_d, ng_d, 0 )
  call mp_gather4d( delp, y(ibeg_delp), imr, jnp, nl, 1, jfirst, jlast, &
                    kfirst, klast, 0   , 0   , 0 )
  if (nc .gt. 1) then
      do m_nc = 2, nc
         call mp_gather4d( q(:,:,:,m_nc), y(ibeg_tr(m_nc-1)), &
                           imr, jnp, nl, 1, jfirst, jlast, &
                           kfirst, klast, ng_d, ng_d, 0 )
      enddo
  endif

! broadcast dependent vector
! --------------------------
  call mp_bcst_n_real( y, m )

#else
  y(ibeg_delp:iend_delp)=reshape(delp(:,jfirst:jlast,:),(/1+iend_delp-ibeg_delp/))
  y(ibeg_q :iend_q )   = reshape(q (:,jfirst:jlast,:,1),(/1+iend_q -ibeg_q /))
  y(ibeg_pt:iend_pt)   = reshape(pt(:,jfirst:jlast,:)  ,(/1+iend_pt-ibeg_pt/))
  y(ibeg_v :iend_v )   = reshape(v (:,jfirst:jlast,:)  ,(/1+iend_v -ibeg_v /))
  y(ibeg_u :iend_u )   = reshape(u (:,jfirst:jlast,:)  ,(/1+iend_u -ibeg_u /))
#endif

end subroutine model2dependent

end module dependent
