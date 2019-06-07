! Parallelized utility routine for computing/printing
! max/min of an input array
!
      subroutine pmaxmin( qname, a, pmin, pmax, im, jt, fac )

      use precision
#if defined( SPMD )
#define CPP_PRT_PREFIX  if(gid.eq.0)
      use mod_comm, only : gid, mp_reduce_max
#else
#define CPP_PRT_PREFIX
#endif

      implicit none

      character*(*)  qname
      integer im, jt
      integer i, j
      real(r8) a(im,jt)

      real(r8) qmin(jt), qmax(jt)
      real(r8) pmax, pmin
      real(r8) fac                     ! multiplication factor
      real(r8) pm1(2)

!$omp parallel do private(i, j, pmax, pmin)

      do j=1,jt
         pmax = a(1,j)
         pmin = a(1,j)
         do i=2,im
            pmax = max(pmax, a(i,j))
            pmin = min(pmin, a(i,j))
         enddo
         qmax(j) = pmax
         qmin(j) = pmin
      enddo
!
! Now find max/min of amax/amin
!
      pmax = qmax(1)
      pmin = qmin(1)
      do j=2,jt
         pmax = max(pmax, qmax(j))
         pmin = min(pmin, qmin(j))
      enddo


#if defined( SPMD )
      pm1(1) = pmax
      call mp_reduce_max(1, pm1)
      pmax=pm1(1)
      pm1(1) = -pmin
      call mp_reduce_max(1, pm1)
      pmin=-pm1(1)
#endif

      CPP_PRT_PREFIX write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

      return
      end
