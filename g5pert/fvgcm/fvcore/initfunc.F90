!  14May2007  Todling   Introduced dyn_prog; global change.
subroutine initfunc( n, x, prog )
  use precision
  use control, only: mod2cont
  use timingModule
  use prognostics, only: dyn_prog

  implicit none

  integer  :: n
  real(r8) :: x(n)
  type(dyn_prog) :: prog

  call mod2cont( n, x, prog )

  call timing_init

end subroutine initfunc
