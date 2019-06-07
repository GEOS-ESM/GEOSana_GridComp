      subroutine zsmean(u,im,jfirst,jlast,us,ux)
      use precision
      implicit none
      integer im, jfirst, jlast
      integer i, j
      real(r8) u(im,jfirst:jlast), us(im,jfirst:jlast)
      real(r8) ux(jfirst:jlast)
      real(r8)    rim

      rim = 1. / float(im)

      do j=jfirst,jlast
        ux(j) = u(im,j)
        do i=1,im-1
          ux(j) = ux(j) + u(i,j)
        enddo
        ux(j) = ux(j) * rim

        do i=1,im
          us(i,j) = u(i,j) -  ux(j)
        enddo

      enddo

      return
      end
