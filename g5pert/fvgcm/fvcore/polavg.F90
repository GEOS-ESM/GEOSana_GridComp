      subroutine polavg(p, im, jm, jfirst, jlast)

      use precision
      implicit none

      integer im, jm, jfirst, jlast
      real(r8) p(im,jfirst:jlast)
      real(r8) sum1

      integer i

      if ( jfirst == 1 ) then 
          sum1 = 0.
        do i=1,im
          sum1 = sum1 + p(i,1)
        enddo
          sum1 = sum1/im

        do i=1,im
          p(i,1) = sum1
        enddo
      endif

      if ( jlast == jm ) then
          sum1 = 0.
        do i=1,im
          sum1 = sum1 + p(i,jm)
        enddo
          sum1 = sum1/im

        do i=1,im
          p(i,jm) = sum1
        enddo
      endif

      return
      end
