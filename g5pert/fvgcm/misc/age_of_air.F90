      subroutine age_of_air(im, jm, km, jfirst, jlast, ng, time, pe, q)

      use precision
      implicit none
      integer im
      integer jm
      integer km
      integer jfirst
      integer jlast  
      integer ng

! q is the age tracer
! Need to be converted to mixing ratio (mass of tracer / dry_air-mass)
! Ignore this inconsistency for now.

      real(r8) q(im,jfirst-ng:jlast+ng,km)
      real(r8) pe(im,km+1,jfirst:jlast)
      real(r8) time                        ! accumulated time since init

! Local
      integer i, j, k
      real(r8) p_source
      real(r8) a
      real(r8) tiny
      parameter ( tiny = 1.e-6 )
      parameter ( p_source = 75000. )
      parameter ( a = 5.e-6 / 60. )
      real(r8) pm

!$omp parallel do private(i, j, k, pm)

      do k=1,km
        do j=jfirst, jlast
          do i=1,im
            pm = 0.5 * ( pe(i,k,j) + pe(i,k+1,j))
            if( time < tiny ) then
                q(i,j,k) = 0.
            elseif( pm >= p_source ) then
                q(i,j,k) = a * time
            endif
          enddo
        enddo          ! j-loop
      enddo             ! k-loop

      return
      end
