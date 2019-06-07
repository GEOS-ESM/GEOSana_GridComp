!-----------------------------------------------------------------------
!BOP
! !ROUTINE: avgp2 --- Average 2D field at the poles
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine avgp2(p,im,jm,jfirst,jlast)
!****6***0*********0*********0*********0*********0*********0**********72
! !USES:
      use precision
      implicit none
! !INPUT PARAMETERS:
      integer(i4) im, jm                   ! Horizontal dimensions
      integer(i4) jfirst, jlast            ! Strip size
! !INPUT/OUTPUT PARAMETERS:
      real(r8) p(im,jfirst:jlast)          ! field
! !DESCRIPTION:
!    Author: Shian-Jiann Lin NASA/GSFC
!
! !REVISION HISTORY:
!   SJL 99.04.13:  Delivery
!   WS  99.11.03:  Documentation; indentation; jfirst:jlast
!
!EOP
!-----------------------------------------------------------------------
!BOC
      integer(i4) i, j
      real(r8) sine(3)
      real(r8) sum1, sum2
      real(r8) sum3, sum4
      real(r8) rim
      real(r8) pi, ph5     ! was DOUBLE PRECISION

      pi  = 4.d0 * datan(1.d0)

      do j=2,3
        ph5  = -0.5d0*pi + (dble(j-1)-0.5d0)*(pi/dble(jm-1))
        sine(j) = sin(ph5)
      enddo

      rim = 1./ float(im)

! WS 99.11.03 : replaced sump2 with 2 calls to F90 SUM

!!!      call sump2(p(1,1),p(1,jm),im,sum1,sum2)
!!!      sum1 = sum1*(1.+sine(2))
!!!      sum2 = sum2*(1.+sine(2))
!!!      call sump2(p(1,2),p(1,jm-1),im,sum3,sum4)
!!!      sum1 = rim * ( sum1 + sum3*(sine(3)-sine(2)) ) / (1.+sine(3))
!!!      sum2 = rim * ( sum2 + sum4*(sine(3)-sine(2)) ) / (1.+sine(3))


      if ( jfirst .eq. 1 ) then
        sum1 = sum(p(1:im,1))*(1.+sine(2))
        sum3 = sum(p(1:im,2))
        sum1 = rim * ( sum1 + sum3*(sine(3)-sine(2)) ) / (1.+sine(3))
        do i=1,im
          p(i,  1)  = sum1
          p(i,  2)  = sum1
        enddo
      endif

      if ( jlast .eq. jm ) then
        sum2 = sum(p(1:im,jm))*(1.+sine(2))
        sum4 = sum(p(1:im,jm-1))
        sum2 = rim * ( sum2 + sum4*(sine(3)-sine(2)) ) / (1.+sine(3))
        do i=1,im
          p(i,jm-1) = sum2
          p(i,jm)   = sum2
        enddo
      endif

      return
!EOC
      end subroutine avgp2
!-----------------------------------------------------------------------
