!-----------------------------------------------------------------------
!BOP
! !ROUTINE: avgc --- Interpolate onto C grid
!
! !INTERFACE:

   subroutine avgc(aa, ac, im, jm, jfirst, jlast, wk)

! !USES:

   use precision

   implicit none

! !INPUT PARAMETERS:
   integer im, jm                      ! Dimensions
   integer jfirst, jlast               ! Latitude strip
   real(r8) aa(im,jfirst-1:jlast+1)    ! Field on A grid (ghosted N*1 S*1)

! !INPUT/OUTPUT PARAMETERS:
! Work array, passed for better cache performance
   real(r8) wk(im,jfirst-1:jlast+1)

! !OUTPUT PARAMETERS:
   real(r8), intent(out):: ac(im,jfirst:jlast+1)  ! Field on C grid (ghosted N for cd_core)

! !DESCRIPTION:
!     Perform averaging of AA at cell centers to AC at cell corners
!     AC may or may not be AA
!
! !CALLED FROM:
!     cd_core
!
! !REVISION HISTORY:
!     SJL 99.04.13:  Delivery
!     WS  99.09.09:  Documentation; indentation, cleaning, jfirst:jlast
!     WS  00.05.14:  Renamed ghost variables as per Kevin's definitions
!
!EOP
!-----------------------------------------------------------------------
!BOC
      integer i, j, js1g1, js2g0, jn1g1

      js1g1 = max(jfirst-1,1)
      js2g0 = max(jfirst,2)
      jn1g1 = min(jlast+1,jm)

     do j=js1g1,jn1g1
           wk(1,j) = aa(1,j) + aa(im,j)
        do i=2,im
           wk(i,j) = aa(i,j) + aa(i-1,j)
        enddo
     enddo

     do j=js2g0,jn1g1
        do i=1,im
           ac(i,j) = wk(i,j) + wk(i,j-1)
        enddo
     enddo

   return
!EOC
   end
!-----------------------------------------------------------------------

