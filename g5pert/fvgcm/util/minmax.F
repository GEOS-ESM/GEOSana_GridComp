      subroutine minmax (name, m, n, l, f)
      
      implicit         none

      character*(*)      name

      integer          m, n, l
      integer          i, j, k

      real             f(m,n,l)
      real             fmax
      real             fmin

      fmax = f(1,1,1)
      fmin = f(1,1,1)

      do k = 1, l
        do j = 1, n
          do i = 1, m
            fmax = max(fmax,f(i,j,k))
            fmin = min(fmin,f(i,j,k))
          end do
        end do
      end do

      write(6,*) name, ' max = ', fmax, ' min = ', fmin
c     print '(a8, 2(1x,a7,e16.8))',
c    .      name, ' max = ', fmax, ' min = ', fmin

      return
      end
