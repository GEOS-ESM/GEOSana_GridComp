      program main

!Usage:
! on Linux with FFC f90: f90 -CcdRR8 d_rst.F

      implicit none

      integer im, jm, km                 ! Dimensions
      integer nq                         ! Number of tracers
c     parameter (im=144, jm=91, km=55)
      parameter (im=288, jm=181, km=55)
c     parameter (im=576, jm=361, km=32)
      parameter ( nq = 1)

      integer iuic                       ! Unit number
      integer iout
      integer nymd, nhms, nstep          ! Time stamps
      real    u(im,jm,km)
      real    v(im,jm,km)
      real    pt(im,jm,km)
      real    delp(im,jm,km)
      real    q(im,jm,km,nq)
      real    ps(im,jm)

      real    tmpu(im,jm,km)

      integer i,j,k,ks
      real pkz(im,jm,km), pk(im,jm,km+1), pe(im,km+1,jm)
      real peln(im,km+1,jm)
      real rair, cpair
      real cappa
      real zvir
      real rh2o
      real ptop
      real tmp, qmax, qmin
!     real gravit
      integer ipf
      integer ndeg

      write(*,*) 'Perform polar filter? (0: No; 1: Yes)'
      read(*,*) ipf
 
!     ptop = 1.                  ! pascal
      write(*,*) 'PTOP=? (Pascal)'
      read(*,*) ptop
      
      cpair  = 1004.64
!     gravit = 9.80616
      rair   = 287.04
      cappa  = rair/cpair

      rh2o   = 4.61e2
      zvir   = rh2o/rair - 1.

      iuic = 81
      iout = 66

      open (unit=iuic,file='d_rst',form='unformatted',status='old')
      open (unit=iout,file='dyn.data',form='unformatted',
     &      status='unknown')

      read(iuic) nstep, nymd, nhms
      write(*,*) nstep, nymd, nhms
      read(iuic) ps, delp, u, v, pt
      if(nq .gt. 0) read(iuic) q

      if ( ipf .eq. 1 ) then
            ndeg  = 60
            call plft2d(im, jm, ps, 2, jm-1, ndeg)
            call polavg(ps, im, jm, 1, jm)
         do k=1,km
            call plft2d(im, jm, delp(1,1,k), 2, jm-1, ndeg)
            call plft2d(im, jm,   pt(1,1,k), 2, jm-1, ndeg)
            call plft2d(im, jm,    u(1,1,k), 2, jm  , ndeg)
            call plft2d(im, jm,    v(1,1,k), 2, jm-1, ndeg)
            call plft2d(im, jm,  q(1,1,k,1), 2, jm-1, ndeg)
! Do pole average
            call polavg(delp(1,1,k), im, jm, 1, jm)
            call polavg(  pt(1,1,k), im, jm, 1, jm)
            call polavg( q(1,1,k,1), im, jm, 1, jm)
         enddo
      write(*,*) 'Warning: original d_rst file will be overwritten'
      rewind iuic
      write(iuic) nstep, nymd, nhms
      write(iuic) ps, delp, u, v, pt
      if(nq .gt. 0) write(iuic) q
      endif

      call pmaxmin('PS ', ps, qmin, qmax, im*jm,1, 0.01)

! Convert pt to tmpu

      do j=1,jm
            tmp = ptop**cappa
            do i=1,im
               pe(i,1,j) = ptop
               pk(i,j,1) = tmp
            enddo

            do k=2,km+1
            do i=1,im
               pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
               pk(i,j,k) = pe(i,k,j)**cappa
            enddo
            enddo
      enddo

        call pkez(1, im, km, 1, jm,
     &            pe, pk, cappa, ks, peln, pkz, .false.)

      do k=1,km
         do j=1,jm
            do i=1,im
               tmpu(i,j,k) = pt(i,j,k) * pkz(i,j,k) /
     &                       ( 1. + zvir*q(i,j,k,1))
            enddo
         enddo

      write(*,*) 'K=',k
      call pmaxmin('Delp ', delp(1,1,k), qmin, qmax, im*jm,1, 0.01)
      call pmaxmin('U ', u(1,1,k), qmin, qmax, im*jm,1, 1.)
      call pmaxmin('V ', v(1,1,k), qmin, qmax, im*jm,1, 1.)
      call pmaxmin('PT ', pt(1,1,k), qmin, qmax, im*jm,1, 1.)
      call pmaxmin('Q ', q(1,1,k,1), qmin, qmax, im*jm,1, 1.)
      call pmaxmin('T ', tmpu(1,1,k), qmin, qmax, im*jm,1, 1.)
      
      enddo

! Output grads readable 32-bit data

      call wrt3d(iout,im,jm,1,ps)
      call wrt3d(iout,im,jm,km,delp)
      call wrt3d(iout,im,jm,km,u)
      call wrt3d(iout,im,jm,km,v)
      call wrt3d(iout,im,jm,km,pt)
      call wrt3d(iout,im,jm,km,tmpu)
      call wrt3d(iout,im,jm,km,q)
      write(*,*) 'End writting grads data'

      end


!-----------------------------------------------------------------------
!BOP
! !ROUTINE: pkez --- Calculate solution to hydrostatic equation
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine pkez(nx, im, km, jfirst, jlast,
     &                pe, pk, akap, ks, peln, pkz, eta)
!****6***0*********0*********0*********0*********0*********0**********72
!
! !USES:
      implicit none

! !INPUT PARAMETERS:
      integer nx                          ! SMP decomposition in x
      integer im, km                  ! Dimensions
      integer jfirst, jlast               ! Latitude strip
      real  pe(im, km+1, jfirst:jlast)    ! Edge pressure
      integer ks
      logical eta     ! Is on ETA coordinate?
                      ! True:  input pe;     output pk, pkz, peln
                      ! False: input pe, pk; output pkz
      real akap

! !INPUT/OUTPUT PARAMETERS:
      real   pk(im,jfirst:jlast,km+1)
      real  pkz(im,jfirst:jlast,km)

! !OUTPUT
      real peln(im, km+1, jfirst:jlast)   ! log pressure (pe) at layer edges

! !DESCRIPTION:
!
!
! !CALLED FROM:
!     te_map and fvgcm 
!
! !REVISION HISTORY:
!
!     WS  99.05.19 : Removed fvcore.h
!     WS  99.07.27 : Limited region to jfirst:jlast
!     WS  99.10.22 : Deleted cp as argument (was not used)
!     WS  99.11.05 : Documentation; pruning of arguments
!     SJL 00.01.02: SMP decomposition in i
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local
      real  pk2(im, km+1)
      real pek
      real lnp
      integer i, j, k
      integer ixj, jp, it, i1, i2

      it = im / nx
      jp = nx * ( jlast - jfirst + 1 )

#if ( defined OpenMP )
!$omp  parallel do
!$omp& default(shared)
!$omp& private(ixj, i1, i2, i, j, k, pek, lnp, pk2)
#endif

#if ( defined IRIX64 )
c$doacross local(ixj, i1, i2, i, j, k, pek, lnp, pk2),
c$&        shared(it, jp)
#endif

! WS 99.07.27 : Limited region to jfirst:jlast

!     do 1000 j=jfirst, jlast
!        i1 = 1
!        i2 = im

      do 1000 ixj=1,jp

         j  = jfirst + (ixj-1) / nx
         i1 =  1 + it * mod(ixj-1, nx)
         i2 = i1 + it - 1

        if ( eta ) then

! <<<<<<<<<<< Eta cordinate Coordinate  >>>>>>>>>>>>>>>>>>>
          pek =     pk(i1,j,1)
          lnp = log(pe(i1,1,j))

          do i=i1,i2
             pk2(i,1)   = pek
            peln(i,1,j) = lnp
          enddo

          if(ks .ne. 0) then
            do k=2, ks+1
              pek = pe(i1,k,j)**akap
              lnp = log(pe(i1,k,j))
              do i=i1,i2
                 pk2(i,k)   = pek
                peln(i,k,j) =  lnp
              enddo
            enddo

            do k=1, ks
              pek = (       pk2(i1,k+1)   - pk2(i1,k))   /
     &              (akap*(peln(i1,k+1,j) - peln(i1,k,j)) )
              do i=i1,i2
                 pkz(i,j,k) = pek
              enddo
            enddo
          endif

          do k=ks+2,km
#if !defined( VECTOR_MATH )
            do i=i1,i2
               pk2(i,k) = pe(i,k,j)**akap
            enddo
#else
            call vlog(pk2(i1,k), pe(i1,k,j), it)
            do i=i1,i2
               pk2(i,k) = akap * pk2(i,k)
            enddo
            call vexp(pk2(i1,k), pk2(i1,k), it)
#endif
          enddo

          do i=i1,i2
             pk2(i,km+1) = pk(i,j,km+1)
          enddo

          do k=ks+2,km+1
            do i=i1,i2
               peln(i,k,j) =  log(pe(i,k,j))
            enddo
          enddo

          do k=ks+1,km
            do i=i1,i2
               pkz(i,j,k) = (pk2(i,k+1) - pk2(i,k)) /
     &                      (akap*(peln(i,k+1,j) - peln(i,k,j)) )
            enddo
          enddo

          do k=2,km
            do i=i1,i2
               pk(i,j,k) = pk2(i,k)
            enddo
          enddo

        else

! <<<<<<<<<<< General Coordinate  >>>>>>>>>>>>>>>>>>>

          pek =     pk(i1,j,1)
          lnp = log(pe(i1,1,j))

          do i=i1,i2
             pk2(i,1) = pek
             peln(i,1,j) = lnp
          enddo

          do k=2,km+1
             do i=i1,i2
                peln(i,k,j) =  log(pe(i,k,j))
                pk2(i,k) =  pk(i,j,k)
             enddo
          enddo

          do k=1,km
             do i=i1,i2
                pkz(i,j,k) = (       pk2(i,k+1) - pk2(i,k) )  /
     &                       (akap*(peln(i,k+1,j) - peln(i,k,j)) )
             enddo
          enddo

        endif
1000  continue

      return
!EOC
      end
!-----------------------------------------------------------------------


! Parallelized utility routine for computing/printing
! max/min of an input array
!
      subroutine pmaxmin( qname, a, pmin, pmax, im, jm, fac )

      implicit none

#if defined( MPI_ON )
#include "mpif.h"
      real tmp 
      integer ierror
#endif

      character*(*)  qname
      integer im, jm
      integer i, j
      real a(im,jm)

      real qmin(jm), qmax(jm)
      real pmax, pmin
      real fac                     ! multiplication factor

#if ( defined OpenMP )
!$omp  parallel do
!$omp& default(shared)
!$omp& private(i,j, pmax, pmin)
#endif

#if ( defined IRIX64 )
c$doacross   local(i,j),
c$&          share(im,jm, pmax, pmin)
#endif

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

#if defined( MPI_ON )
      call mpi_allreduce( pmax, tmp, 1, MPI_DOUBLE_PRECISION,
     &                    MPI_MAX, comm, ierror )
      pmax = tmp
      call mpi_allreduce( pmin, tmp, 1, MPI_DOUBLE_PRECISION,
     &                    MPI_MIN, comm, ierror )
      pmin = tmp
#endif

      write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

      return
      end

      subroutine wrt3d(iout,im,jm,km,a3)

      implicit none
      integer iout, im, jm, km
      integer i, j, k
      real a3(im,jm,km)
      real*4 a2(im,jm)

! Flip vertically
!!!        do 50 k=1,km
        do 50 k=km,1, -1
          do j=1,jm
            do i=1,im
              if(abs(a3(i,j,k)) .lt. 1.e-25) then
                a2(i,j) = 0.
              else
                a2(i,j) = a3(i,j,k)
              endif
            enddo
          enddo
          write(iout) a2
50      continue

      return
      end

      subroutine plft2d(im, jm, p, JS, JN, ndeg)
!
! This is a weak LOCAL polar filter.
! Developer: Shian-Jiann Lin

      implicit none

      integer im
      integer jm
      integer js, jn, ndeg
      real*8 p(im,jm)

      integer i, j, n, ideg, jj, jc
      real*8 cosp(jm),cose(jm)
      real*8 a(0:im/2+1)

      real*8 sine(jm),sinp(jm)
      real*8, allocatable, save :: se(:), sc(:)

      real*8 pi, dp, dl, e0, ycrit, coszc, smax, rn, rn2, esl, tmp

      data IDEG /0/

      if(IDEG .ne. ndeg) then
	IDEG = ndeg
C  (e0 = 2.6)
      e0 = 0.5 * sqrt(27.)
      PI = 4. * ATAN(1.)

      allocate( sc(jm), se(jm))

      call setrig(im, jm, DP, DL, cosp, cose, sinp, sine)

      ycrit = IDEG*PI/180.
      coszc = cos(ycrit)

      smax = (jm-1)/2
      write(6,*) 'Critical latitude in local pft = ',ndeg

      a(0) = 1.
      do n=1,im/2+1
      rn = n
      rn2 = 2*n
      a(n) = sqrt(rn2+1.) * ((rn2+1.)/rn2)**rn
      enddo

      do j=2,jm-1
      sc(j) = coszc / cosp(j)

      IF(sc(j) .gt.1. .and. sc(j) .le. 1.5 ) THEN
         esl = 1./ sc(j)
         sc(j) =  1. +  (1.-esl) / (1.+esl)
      ELSEIF(sc(j) .gt.1.5 .and. sc(j) .le. e0 ) THEN
         esl = 1./ sc(j)
         sc(j) =  1. + 2./ (27.*esl**2 - 2.)
      ELSEIF(sc(j).gt. e0) THEN
C Search
      do jj=1,im/2
      if(sc(j).le. a(jj)) then
      jc = jj
      write(*,*) 'jc=', jc
      goto 111
      endif
      enddo

      jc = im/2 + 1
111   continue
      tmp = ((sc(j) - a(jc-1))/(a(jc) - a(jc-1)))**0.25
      sc(j) =  jc + min(1., tmp)
c     sc(j) =  min(smax,sc(j))
      ENDIF
      enddo
C ====================================================
      do j=2,jm
      se(j) = coszc / cose(j)
      IF(se(j) .gt.1. .and. se(j) .le. 1.5 ) THEN
         esl = 1./ se(j)
         se(j) =  1. + (1.-esl) / (1.+esl)
      ELSEIF(se(j) .gt.1.5 .and. se(j) .le. e0 ) THEN
         esl = 1./ se(j)
         se(j) =  1. + 2./ (27.*esl**2 - 2.)
      ELSEIF(se(j).gt. e0) THEN
C Search
      do jj=1,im/2
      if(se(j) .le. a(jj)) then
      jc = jj
      goto 222
      endif
      enddo

      jc = im/2 + 1
222   continue
      tmp = ((se(j) - a(jc-1))/(a(jc) - a(jc-1)))**0.25
      se(j) =  jc + min(1., tmp)
c     se(j) =  min(smax,se(j))
      ENDIF
      enddo

      do i=1,im
        se( 2) = sc(2)
        se(jm) = sc(jm-1)
      enddo

      do j=2,jm-1
         write(*,*) j,sc(j)
      enddo
      ENDIF

      if(JN .eq. (jm-1) ) then
C Cell-centered variables
         call lpft(im, jm, p, 2, jm-1, Sc)
      else
C Cell-edge variables
         call lpft(im, jm, p, 2, jm, Se)
      endif
      return
      end


      subroutine lpft(im, jm, p, j1, j2, s)
      implicit none

      integer im, jm, j1, j2
      real*8  p(im,jm)
      real*8  s(jm)

! Local
      integer i, j, n, nt

      real*8 ptmp(0:im+1)
      real*8 q(0:im+1)
      real*8 frac, rsc, bt

      do 2500 j=j1,j2
      if(s(j) .gt. 1.02) then

        NT  = INT(S(j))
        frac = S(j) - NT
        NT = NT-1

        rsc = 1. / (1.+frac)
        bt = 0.5 * frac

        do i=1,im
           ptmp(i) = p(i,j)
        enddo

        ptmp(0)    = p(im,j)
        ptmp(im+1) = p(1 ,j)

        if( NT .lt. 1 ) then
          do i=1,im
             p(i,j) = rsc * (ptmp(i) + bt*(ptmp(i-1)+ptmp(i+1)))
          enddo
        else
          do i=1,im
             q(i) = rsc * (ptmp(i) + bt*(ptmp(i-1)+ptmp(i+1)))
          enddo

           do 500 N=1,NT
              q(0) = q(im)
              do i=1,im
                 ptmp(i) = q(i) + q(i-1)
              enddo
                 ptmp(im+1) = ptmp(1)

             if ( n .eq. nt ) then
               do i=1,im
                  p(i,j) = 0.25*(ptmp(i) + ptmp(i+1))
               enddo
             else
               do i=1,im
                  q(i) = 0.25*(ptmp(i) + ptmp(i+1))
               enddo
             endif
500        continue
        endif
      endif
2500  continue

      return
      end

      subroutine setrig(im, jm, dp, dl, cosp, cose, sinp, sine)

      implicit none

      integer im, jm
      integer j, jm1
      real*8 sine(jm),cosp(jm),sinp(jm),cose(jm)
      real*8 dp, dl
      double precision pi, ph5

      jm1 = jm - 1
      pi  = 4.d0 * datan(1.d0)
      dl  = (pi+pi)/dble(im)
      dp  = pi/dble(jm1)

      do 10 j=2,jm
         ph5  = -0.5d0*pi + (dble(j-1)-0.5d0)*(pi/dble(jm1))
10    sine(j) = dsin(ph5)

      cosp( 1) =  0.
      cosp(jm) =  0.

      do 80 j=2,jm1
80    cosp(j) = (sine(j+1)-sine(j)) / dp

c Define cosine at edges..

      do 90 j=2,jm
90    cose(j) = 0.5 * (cosp(j-1) + cosp(j))
      cose(1) = cose(2)

      sinp( 1) = -1.
      sinp(jm) =  1.

      do 100 j=2,jm1
100   sinp(j) = 0.5 * (sine(j) + sine(j+1))

      return
      end

      subroutine polavg(p,im,jm,jfirst,jlast)

      implicit none
      integer im, jm, jfirst, jlast
      integer i
      real*8 p(im,jfirst:jlast)
      real*8 sum1, sum2

      if ( jfirst .eq. 1 ) then 
        sum1 = SUM(p(1:im,1))/im
        do i=1,im
          p(i,  1) = sum1
        enddo
      endif

      if ( jlast .eq. jm ) then
        sum2 = SUM(p(1:im,jm))/im

        do i=1,im
          p(i,jm) = sum2
        enddo
      endif
      return
      end




