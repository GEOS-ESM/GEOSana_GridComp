      subroutine init2dz(jm, km, jfirst, jlast, q, def)

      use precision
      implicit      none

      integer jm, km, jfirst, jlast
      real(r8)    q(jfirst:jlast, km)
      real(r8)    def
      integer j,  k

!$omp  parallel do        &
!$omp  default(shared)    &
!$omp  private(j,k)

      do k=1,km
        do j=jfirst,jlast
          q(j,k) = def
        enddo
      enddo

      return
      end

