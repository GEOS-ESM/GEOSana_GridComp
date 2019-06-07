      subroutine par_vecsum(jm, jfirst, jlast, InVector, te0)

! !ROUTINE: par_vecsum --- Calculate vector sum bit-wise consistently

#if defined ( SPMD )
      use mod_comm, only : mp_sum1d
#endif
      use precision
      implicit none
! Input:
      integer jm                   ! global latitudes
      integer jfirst               ! first latitude on this PE
      integer jlast                ! last latitude on this PE
      real InVector(jfirst:jlast)  ! input vector to be summed

! OUTPUT:
      real(r8) te0                 ! sum of all vector entries
! Local
      integer j

#if defined ( SPMD ) 
      call mp_sum1d(jm, jfirst, jlast, InVector, te0)
#else
      te0 = 0.0
      do j=1,jm
        te0 = te0 + InVector(j)
      enddo
#endif

      return
      end
