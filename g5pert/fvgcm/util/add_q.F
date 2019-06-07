      program main
!
! Util to add tracer to d_rst file (eg, age-of-air tracer)
!

!Usage:
! on Linux with FFC f90: f90 -CcdRR8 add_q.F

      implicit none

      integer im, jm, km                 ! Dimensions
      integer nq                         ! Number of tracers

!     parameter (im= 72, jm=46, km=18)
      parameter (im=144, jm=91, km=55)
!     parameter (im=288, jm=181, km=55)
!     parameter (im=576, jm=361, km=32)

      integer nq_old
      integer nq_new
      parameter ( nq_old = 1)
      parameter ( nq_new = 2)

      integer iuic                       ! Unit number
      integer iout

      integer nymd, nhms, nstep          ! Time stamps
      real    u(im,jm,km)
      real    v(im,jm,km)
      real    pt(im,jm,km)
      real    delp(im,jm,km)
      real    q(im,jm,km,nq_old)
      real    qq(im,jm,km,nq_new)
      real    ps(im,jm)

      integer i,j,k, iq
      real tmp, qmax, qmin

      iuic = 81
      iout = 66

      open (unit=iuic,file='d_rst',form='unformatted',status='old')
      open (unit=iout,file='d_rst_q',form='unformatted',
     &      status='unknown')

!--------------------------------------
      read(iuic)  nstep, nymd, nhms
      write(*,*)  nstep, nymd, nhms
      nstep = 2 
      write(iout) nstep, nymd, nhms

      read(iuic)  ps, delp, u, v, pt
      write(iout) ps, delp, u, v, pt
      call pmaxmin('PS ', ps, qmin, qmax, im*jm,1, 0.01)

      if(nq_old > 0) read(iuic) q
!--------------------------------------

! Change date/time
!     nymd = 
      

! Add tracer
      do iq=1,nq_old
      do k=1,km
         do j=1,jm
            do i=1,im
               qq(i,j,k,iq) = q(i,j,k,iq)
            enddo
         enddo
      enddo
      enddo

      do iq=nq_old+1, nq_new
      do k=1,km
         do j=1,jm
            do i=1,im
               qq(i,j,k,iq) = 0.
            enddo
         enddo
      enddo
      enddo

      write(iout) qq

      end

! Parallelized utility routine for computing/printing
! max/min of an input array
!
      subroutine pmaxmin( qname, a, pmin, pmax, im, jm, fac )

      implicit none

      character*(*)  qname
      integer im, jm
      integer i, j
      real a(im,jm)

      real qmin(jm), qmax(jm)
      real pmax, pmin
      real fac                     ! multiplication factor

      do j=1,jm
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
         do j=2,jm
            pmax = max(pmax, qmax(j))
            pmin = min(pmin, qmin(j))
         enddo

      write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

      return
      end
