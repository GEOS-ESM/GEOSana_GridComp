!  14May2007  Todling   Introduced dyn_prog; global change.

subroutine setfunc( n, m, prog )

  use control, only: control_number
  use stepon, only: stepon_initialize
  use dependent, only: dependent_number
  use prognostics, only : dyn_prog

  implicit none
  integer   :: n, m
  type(dyn_prog) :: prog

  call stepon_initialize ( prog )
  call control_number( n )
  call dependent_number( m )

end subroutine setfunc
!-----------------------------------------------------------------------
