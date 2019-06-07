! 15Aug2006 - Todling - f77 Interface to SPG calls
subroutine dspg_evalf ( n, x, f, rc )
use m_cnop, only : cnop_evalf
implicit none
integer n, rc
real x(n), f
call cnop_evalf ( n, x, rc, f=f )
end

subroutine dspg_evalg ( n, x, g, rc )
use m_cnop, only : cnop_evalg
implicit none
integer n, rc
real x(n), g(n)
call cnop_evalg ( n, x, g, rc )
end

subroutine dspg_proj ( n, x, rc )
use m_cnop, only : cnop_proj
implicit none
integer n, rc
real x(n)
call cnop_proj ( n, x, rc )
end
