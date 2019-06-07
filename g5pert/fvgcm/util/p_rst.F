      program main
      implicit none

      integer im, jm                 ! Dimensions
      integer nq                         ! Number of tracers
      integer mx
      parameter (im=144, jm=91)
      parameter ( nq = 1)
      parameter ( mx = 4)

      integer iu                       ! Unit number
      integer iout
      integer i,j,k

      real oro(im, jm)        !CCM surface flag
      real pblht(im, jm)      !pbl height (m?)
      real tpert(im, jm)      !temperature pert.
      real qpert(im, nq, jm)

      real srfrad(im,jm)
      real asdir(im, jm)      ! albedo: shortwave, direct
      real asdif(im, jm)      ! albedo: shortwave, diffuse
      real aldir(im, jm)      ! albedo: longwave, direct
      real aldif(im, jm)      ! albedo: longwave, diffuse
      real lwup(im, jm)       ! longwave up radiative flux
      real ts(im, jm)         ! sfc temp (merged w/ocean if coupled)
      real tssub(im,mx,jm)
      real sicthk(im,jm)
      real snowh(im,jm)      ! ccm snow depth (liquid water)

      real fsns(im,jm)
      real precst(im,jm)
      real precc(im,jm)
      real precl(im,jm)
      real flwds(im,jm)

      real sols(im,jm)
      real soll(im,jm)
      real solsd(im,jm)
      real solld(im,jm)

      real ts0(im,jm)

      real qmax, qmin

      iu = 81
      iout = 67

      open (unit=iu,file='p_rst',status='old')
      open (unit=iout,file='phys.data',status='unknown')

            read(iu) oro
            read(iu) pblht
            read(iu) tpert
            read(iu) qpert
            read(iu) snowh
            read(iu) lwup

            read(iu) srfrad
            read(iu) fsns
            read(iu) precst
            read(iu) flwds
            read(iu) precc
            read(iu) precl
            read(iu) sols
            read(iu) soll
            read(iu) solsd
            read(iu) solld

            read(iu) ts
            read(iu) tssub
            read(iu) sicthk
 
            read(iu) asdir
            read(iu) asdif
            read(iu) aldir
            read(iu) aldif

      call pmaxmin('ORO   ',  oro, qmin, qmax, im*jm,1, 1.)
      call pmaxmin('pblht ', pblht, qmin, qmax, im*jm,1, 1.)
      call pmaxmin('snowh ', snowh, qmin, qmax, im*jm,1, 1.)
      call pmaxmin('precc ', precc, qmin, qmax, im*jm,1, 86400000.)
      call pmaxmin('precl ', precl, qmin, qmax, im*jm,1, 86400000.)
      call pmaxmin('ts    ',    ts, qmin, qmax, im*jm,1, 1.)

! Output grads readable 32-bit data

      call wrt3d(iout,im,jm,1,oro)
      call wrt3d(iout,im,jm,1,pblht)
      call wrt3d(iout,im,jm,1,tpert)
      call wrt3d(iout,im,jm,1,qpert)
      call wrt3d(iout,im,jm,1,snowh)
      call wrt3d(iout,im,jm,1,lwup )
      call wrt3d(iout,im,jm,1,srfrad)
      call wrt3d(iout,im,jm,1,fsns )
      call wrt3d(iout,im,jm,1,precst)
      call wrt3d(iout,im,jm,1,flwds)
      call wrt3d(iout,im,jm,1,precc)
      call wrt3d(iout,im,jm,1,precl)
      call wrt3d(iout,im,jm,1,sols )
      call wrt3d(iout,im,jm,1,soll )
      call wrt3d(iout,im,jm,1,solsd)
      call wrt3d(iout,im,jm,1,solld)
      call wrt3d(iout,im,jm,1,ts   )

      do k=1,mx
         do j=1,jm
            do i=1,im
               ts0(i,j) = tssub(i,k,j)
            enddo
         enddo
         call wrt3d(iout,im,jm,1,ts0)
      enddo

      call wrt3d(iout,im,jm,1,sicthk)
      call wrt3d(iout,im,jm,1,asdir)
      call wrt3d(iout,im,jm,1,asdif)
      call wrt3d(iout,im,jm,1,aldir)
      call wrt3d(iout,im,jm,1,aldif)

      write(6,*) 'End writting grads data'

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

#define CPP_PRT_PREFIX

      implicit none

      character*(*)  qname
      integer im, jm
      integer i, j
      real a(im,jm)

      real qmin(jm), qmax(jm)
      real pmax, pmin
      real fac                     ! multiplication factor
      real pm1(2), pm2(2)

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

      CPP_PRT_PREFIX write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

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
