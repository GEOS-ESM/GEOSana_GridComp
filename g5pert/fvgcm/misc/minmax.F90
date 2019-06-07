      subroutine minmax (name, f, m, n, l)

      use precision
      implicit         none

      character*(*)        name

      integer(i4)          m, n, l
      integer(i4)          i, j, k

      real(r8)             f(m,n,l)
      real(r8)             fmax
      real(r8)             fmin
      real(r8)             mean
      real(r8)             big
      parameter (big = 1.e20)
      integer(i4)          count

      fmax = f(1,1,1)
      fmin = fmax
      mean = 0.
      count = 0

      do k = 1, l
        do j = 1, n
          do i = 1, m
            fmax = max(fmax,f(i,j,k))
            fmin = min(fmin,f(i,j,k))
            if( abs(f(i,j,k)) .lt. big ) then
                mean = mean + f(i,j,k)
                count = count + 1
            endif
          end do
        end do
      end do

      if( count .ne. 0 ) mean = mean / count

      write(*,*) name, ' max = ', fmax, ' min = ', fmin, ' mean=',mean

      return
      end

