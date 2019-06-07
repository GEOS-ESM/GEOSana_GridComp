      subroutine xpavg(p,im)

      use precision

      implicit none
      integer im
      integer i
      real(r8) p(im)
      real(r8) sum1

         sum1 = 0.
      do i=1,im
         sum1 = sum1 + p(i)
      enddo

      sum1 = sum1/im

      do i=1,im
         p(i) = sum1
      enddo

      return
      end
