 subroutine fillz(im, km, nq, q, dp)

! !USES:
 use precision
 implicit none

! !INPUT PARAMETERS:
   integer, intent(in) :: im                ! No. of longitudes
   integer, intent(in) :: km                ! No. of levels
   integer, intent(in) :: nq                ! Total number of tracers
   real(r8), intent(in) ::  dp(im,km)       ! pressure thickness

! !INPUT/OUTPUT PARAMETERS:
   real(r8), intent(inout) :: q(im,km,nq)   ! tracer mixing ratio

! !DESCRIPTION:
!   Check for "bad" data and fill from east and west neighbors
!
! !BUGS:
!   Currently this routine only performs the east-west fill algorithm.
!   This is because the N-S fill is very hard to do in a reproducible
!   fashion when the problem is decomposed by latitudes.
!
! !REVISION HISTORY:
!   00.04.01   Lin        Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
   integer i, k, ic
   real(r8) qup, qly, dup

   do ic=1,nq
! Top layer
      do i=1,im
         if( q(i,1,ic) < 0.) then
             q(i,2,ic) = q(i,2,ic) + q(i,1,ic)*dp(i,1)/dp(i,2)
             q(i,1,ic) = 0.
          endif
      enddo

! Interior
      do k=2,km-1
         do i=1,im
         if( q(i,k,ic) < 0. ) then
! Borrow from above
             qup =  q(i,k-1,ic)*dp(i,k-1)
             qly = -q(i,k  ,ic)*dp(i,k  )
             dup =  min( 0.5*qly, qup )        !borrow no more than 50%
             q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1) 
! Borrow from below: q(i,k,ic) is still negative at this stage
             q(i,k+1,ic) = q(i,k+1,ic) + (dup-qly)/dp(i,k+1) 
             q(i,k  ,ic) = 0.
          endif
          enddo
      enddo
 
! Bottom layer
      k = km
      do i=1,im
         if( q(i,k,ic) < 0.) then
! Borrow from above
             qup =  q(i,k-1,ic)*dp(i,k-1)
             qly = -q(i,k  ,ic)*dp(i,k  )
             dup =  min( qly, qup )
             q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1) 
             q(i,k,ic) = 0.
          endif
      enddo
   enddo
!EOC
end subroutine fillz
