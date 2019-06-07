      subroutine aoa_wrt(im, jm, km, jfirst, jlast, ng,  time,   &
                         q, iout)
      use precision
      implicit none

      integer im
      integer jm
      integer km
      integer jfirst
      integer jlast  
      integer ng
! q is "specific humidity"
! Need to be converted to mixing ratio (mass of tracer / dry_air-mass)
! Ignore this inconsistency for now.
      real(r8) q(im,jfirst-ng:jlast+ng,km)
      real(r8) time                        ! accumulated time since init
      integer iout                     ! unit to put data

! Local
      integer i, j, k
      real(r8)  qx(jfirst:jlast,km)
      real(r8) age(jfirst:jlast,km)
      real(r8) p_source
      real(r8) a
      real(r8) ra
      real(r8) tiny
      real(r8) ryear
      parameter ( tiny = 1.e-6 )
      parameter ( p_source = 75000. )
      parameter ( a = 5.e-6 / 60., ra = 1./a )
      parameter (ryear = 1./(365.*24.*3600) )
      real(r8) qsum
      real(r8) rim
      real(r8) qmax, qmin

      rim = 1./float(im)

!$omp  parallel do             &
!$omp  default(shared)         &
!$omp  private(i, j, k, qsum)

      do k=1,km
        do j=jfirst, jlast
          qsum = q(1,j,k)
          do i=2,im
            qsum = qsum + q(i,j,k) 
          enddo
          qx(j,k) =   qsum * rim
          age(j,k) =  (time - ra * qx(j,k)) * ryear
        enddo          ! j-loop
      enddo             ! k-loop

! Zonal mean tracer & age (years)

      call wrt3dr(iout, 1, jm, km,  qx, jfirst, jlast)
      call wrt3dr(iout, 1, jm, km, age, jfirst, jlast)
      call pmaxmin('Age-of-Air', age, qmin, qmax, jlast-jfirst+1, km, 1.)

      return
      end
