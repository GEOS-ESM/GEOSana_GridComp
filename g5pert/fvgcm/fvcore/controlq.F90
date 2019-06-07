!  14May2007  Todling   Introduced dyn_prog; global change.
module controlq

  use precision
  use prognostics, only: dyn_prog
  use prognostics, only : imr,jnp,nl
  use prognostics, only : jfirst,jlast
  use prognostics, only : kfirst,klast
  use stepon, only: ng_d, ng_s
#ifdef SPMD
  use stepon,   only: ng_d, ng_s
  use mod_comm, only: mp_scatter4d, gid
  use mod_comm, only: mp_gather4d
#endif

  implicit none

  PRIVATE

! public mapvars
  public cont2mod
  public control_number 
  public mod2cont

  public ibeg_q   , iend_q

  interface mapvars ; module procedure mapvars_    ; end interface

  integer, save :: ibeg_q   , iend_q
  integer, save :: ndim

  logical, save :: iamset = .false.

  real(r8), pointer :: q(:,:,:,:)

contains
subroutine mapvars_ ( imr, jnp, nl )

  integer, intent(in) :: imr
  integer, intent(in) :: jnp
  integer, intent(in) :: nl 

  ibeg_q  = 1
  iend_q  = ibeg_q  + imr * jnp * nl - 1

end subroutine mapvars_

subroutine cont2mod( n, x, prog )

  implicit none
  integer  :: n
  real(r8) :: x(n)
  type(dyn_prog), TARGET :: prog
  real(r8) :: qhelp(imr,jfirst     :jlast     ,nl)

  q => prog%q

#ifdef SPMD

  call mp_scatter4d( x(ibeg_q)   , qhelp, imr, jnp, nl, 1, jfirst, jlast, &
                              kfirst, klast, 0   , 0   , 0 )
  q   (:,jfirst:jlast,:,1) = qhelp(:,jfirst:jlast,:)

#else

  q(:,jfirst:jlast,:,1) = reshape(x(ibeg_q :iend_q ),(/ imr,jnp,nl /))

#endif

end subroutine cont2mod

subroutine mod2cont( n, x, prog )

  implicit none
  integer  :: n
  real(r8) :: x(n)
  type(dyn_prog), TARGET :: prog

  q => prog%q

  !-----------------------------------------------------------------------

#ifdef SPMD
  call mp_gather4d( q(:,:,:,1), x(ibeg_q), &
                    imr, jnp, nl, 1, jfirst, jlast, &
                    kfirst, klast, ng_d, ng_d, 0 )
#else
  x(ibeg_q :iend_q )   = reshape(q (:,jfirst:jlast,:,1),(/1+iend_q -ibeg_q /))
#endif

end subroutine mod2cont

subroutine control_number( n )
  implicit none
  integer :: n

  if (iamset) then
      n = ndim
      return
  endif

  n = 0
  n = n + imr * jnp * nl

  ndim = n
  print *,' #control = ', n

  call mapvars_ ( imr, jnp, nl )

  iamset = .true.

end subroutine control_number

end module controlq
