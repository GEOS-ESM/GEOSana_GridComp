!  14May2007  Todling   Introduced dyn_prog; global change.
subroutine postfunc( prog )
!********************************************************************
  use precision
  use stepon, only: stepon_finalize
  use prognostics, only: dyn_prog
  use timingModule

  implicit none

  integer  :: n, m
  type(dyn_prog) :: prog

  call timing_prt

  call stepon_finalize ( prog )

end subroutine postfunc
