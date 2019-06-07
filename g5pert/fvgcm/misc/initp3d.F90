      subroutine initp3d(im, jfirst, jlast, kfirst, klast, q, def)

      use precision
      implicit none

      integer im, jfirst, jlast, kfirst, klast
      real(r8)    q(im,kfirst:klast,jfirst:jlast)

      real(r8)    def
      integer i, j, k

!$omp  parallel do         &
!$omp  default(shared)     &
!$omp  private(i,j,k)

      do j=jfirst,jlast
        do k=kfirst,klast
          do i=1,im
            q(i,k,j) = def
          enddo
        enddo
      enddo

      return
      end
