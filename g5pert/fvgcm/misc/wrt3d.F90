      subroutine wrt3d(iout, nrec, im, jm, km, a3, jfirst, jlast)

      use precision

#if defined( SPMD )
      use mod_comm, only: mp_wrt3d
#endif

      implicit none
      integer iout, im, jm, km, nrec
      real r_zero
      parameter ( r_zero = 1.e-25 )
      integer i, j, k
      integer jfirst, jlast

      real(r8) a3(im, jfirst:jlast, km)    ! local slice
      real*4 a2(im,jm)

#if defined( SPMD  )
      call mp_wrt3d(iout, nrec, r_zero, a3, im, jm, km, &
                    jfirst, jlast, 1, km, 0)
#else

!-------------------------------------
! This is a 1st attemp in parallel I/O
!-------------------------------------
!**** SLOWER !!!  !$omp parallel do private(i, j, k, a2)
      do 50 k=1,km
         do j=1,jm
            do i=1,im
              if( abs(a3(i,j,k)) < r_zero ) then
                 a2(i,j) = 0.
              else
                 a2(i,j) = a3(i,j,k)
              endif
            enddo
         enddo
         write(iout, rec=nrec+k) a2
50    continue
#endif
      nrec = nrec + km
      return
      end
