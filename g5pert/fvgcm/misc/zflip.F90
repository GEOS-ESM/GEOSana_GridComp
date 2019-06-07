!****6***0*********0*********0*********0*********0*********0**********72
      subroutine zflip(q,im,km,nc)
!****6***0*********0*********0*********0*********0*********0**********72
! This routine flip the array q (in the vertical).
      use precision
      implicit none
      integer(i4) im, km, nc
      integer(i4) i, k, ic
      real(r8) q(im,km,nc)
!
      real(r8) qtmp
!
      do ic = 1, nc
        do i = 1, im
          do k = 1, (km+1)/2
            qtmp = q(i,k,ic)
            q(i,k,ic) = q(i,km+1-k,ic)
            q(i,km+1-k,ic) = qtmp
          end do
        end do
      end do

      return
      end
