!--------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!--------------------------------------------------------------------------
!BOP
!
! !ROUTINE: OP --- Matrix Vector multiplication routine for Parallel NAG 
!
! !INTERFACE:
!
subroutine op(n,p,q,r,comm)
! !USES:
  use precision
  use m_admtlm
  implicit none

! !INPUT PARAMETERS:
  integer, intent(in) :: n,comm
  real (r8) , intent(in) :: p(n),q(n)
  real (r8) , intent(inout):: r(n)

!
! !DESCRIPTION: Matrix Vector multiplication routine for Parallel NAG
!
! !REVISION HISTORY:
!  23Jul2007 Kim    - Initial implementation
!  23Aug2007 Todling- Declared implicit none; rid of atx93
!
!EOP
!---------------------------------------------------------------------------
  integer i

  do i= 1, n
     r(i)= q(i) 
  enddo

! Call for the TLM and ADM routines
  call admtlm(r,n)

return
end subroutine op

!--------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1, GEOS/DAS !
!--------------------------------------------------------------------------
!BOP
!
! !ROUTINE: OPM --- Strandard matrix operation routine for Parallel NAG 
!
! !INTERFACE:
subroutine opm(n,a,b,comm)
! !USES:
  use precision
  implicit none

! !INPUT PARAMETERS:
  integer, intent(in) :: n,comm
  real(r8) , intent(in) :: a(n)
  real(r8) , intent(out):: b(n)

!
! !DESCRIPTION: Standard Matrix operation option of planso solver
!
! !REVISION HISTORY:
!  23Jul2007 Kim    - Initial implementation
!  23Aug2007 Todling- Declared implicit none
!
!EOP
!---------------------------------------------------------------------------
  integer i

  do i= 1, n
     b(i)= a(i)
  enddo

return
end subroutine opm


