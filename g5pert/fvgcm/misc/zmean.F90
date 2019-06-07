      subroutine zmean(u,im,jm,jfirst,jlast,ux)
      use precision
      implicit none
      integer im, jm, jfirst, jlast
      integer i, j
      real(r8) u(im,jfirst:jlast),ux(jfirst:jlast)
      real(r8) rim

      rim = 1./ float(im)

      do j=jfirst,jlast
         ux(j) = u(im,j)
         do i=1,im-1
            ux(j) = ux(j) + u(i,j)
         enddo
         ux(j) = ux(j) * rim
      enddo

      return
      end
