!-----------------------------------------------------------------------
!BOP
! !ROUTINE: rst_dyn --- Read or Write the restart file
!
! !INTERFACE:
!****6***0*********0*********0*********0*********0*********0**********72
      subroutine rst_dyn(itp, iuic, dflnm, job, pick,
     &                   im, jm, km, kfirst, klast,
     &                   jfirst, jlast, ccm_on, delp, u, v, pt,
     &                   nq, q, ps, ptop, nymd, nhms, nstep, ng_s, ng_d)
!****6***0*********0*********0*********0*********0*********0**********72
!
! !USES:
#if defined( SPMD )
#define CPP_PRT_PREFIX  if(gid.eq.0)
        use mod_comm
#else
#define CPP_PRT_PREFIX
#endif
      implicit none

! !INPUT PARAMETERS:
      integer im, jm, km                 ! Dimensions
      integer jfirst, jlast              ! latitude strip
      integer kfirst, klast              ! vertical strip
      logical ccm_on                     ! running with CCM physics?
      integer nq                         ! Number of tracers
      integer iuic                       ! Unit number
      character*(*) dflnm                ! dynamics filename
      character*(*) job                  ! job name
      logical, intent(in) :: pick        ! determines whether or not to pick particular date/time       
      integer itp                        ! Read (>=0) or Write (<0)
      integer ng_s, ng_d

! !OUTPUT PARAMETERS:
      integer nymd, nhms, nstep          ! Time stamps

      real :: ps(im,jfirst:jlast)
      real :: u(im,jfirst-ng_d:jlast+ng_s,kfirst:klast)
      real :: v(im,jfirst-ng_s:jlast+ng_d,kfirst:klast)
      real :: pt(im,jfirst-ng_d:jlast+ng_d,kfirst:klast)
      real :: delp(im,jfirst:jlast,kfirst:klast)
      real :: q(im,jfirst-ng_d:jlast+ng_d,kfirst:klast,nq)

! !DESCRIPTION:
!     Gather the data if necessary and write the restart file
!
! !REVISION HISTORY:
!   WS  99.09.07:  Now calls readrst and writerst
!   WS  99.09.08:  Documentation; added im,jm,km as arguments
!   WS  99.10.28:  Switched from bcast to scatter
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local variables

      integer i, j, k, kext
      integer ic
      real ptop
      real pmax, pmin

      kext = klast-kfirst+1

      if(itp.ge.0) then

      call readrst(iuic, dflnm, job, pick, im, jm, km, nq, jfirst, jlast,
     &             nstep, nymd, nhms, ps, delp, u, v, pt, q, ng_s, ng_d)


      if(itp .eq. 0) then
         nstep = 0
         CPP_PRT_PREFIX write(6,*) 'nstep reset to zero in rst_dyn()'
      endif

!BW: will modify this part later


      call pmaxmin1('U ', u, pmin, pmax, im, jm, km, jfirst, jlast,
     &              ng_d, ng_s, 1.)
      call pmaxmin1('V ', v, pmin, pmax, im, jm, km, jfirst, jlast,
     &              ng_s, ng_d, 1.)
      call pmaxmin1('PT', pt, pmin, pmax, im, jm, km, jfirst, jlast,
     &              ng_d, ng_d, 1.)

      if(nq .ne. 0) then
      do ic=1,nq
        call pmaxmin1('Q ', q(1,jfirst-ng_d,kfirst,ic), pmin, pmax, 
     &                im, jm, km, jfirst, jlast, ng_d, ng_d, 1.)
      enddo
      endif

	do k=kfirst,klast

        if ( ccm_on .and. nstep .eq. 0) then
          do j=jfirst,jlast
            do i=1,im
              q(i,j,k,1) = max(1.2e-12, q(i,j,k,1))
            enddo
          enddo

      	call polavg(delp(1,jfirst,k),  im, jm, jfirst, jlast)
      	call polavg(  pt(1,jfirst,k),  im, jm, jfirst, jlast)
      	call polavg(  q(1,jfirst,k,1), im, jm, jfirst, jlast)
        endif
	enddo

      endif

      return
!EOC
      end
!-----------------------------------------------------------------------

       subroutine pmaxmin1(qname, a, pmin, pmax, im, jm, km, 
     &            jfirst, jlast, ng_d, ng_s, fac)

       use precision
       implicit none
       character*(*)  qname
       integer im, jm, km, jfirst, jlast, ng_d, ng_s
       real(r8) a(im,jfirst-ng_d:jlast+ng_s,km)
       real(r8) tmp(im, jfirst:jlast,km)
       real(r8) pmax, pmin, fac
       integer i, j, k

!$omp parallel do private(i, j, k)
      do k=1, km
         do j=jfirst, jlast
            do i=1, im
                tmp(i,j,k)=a(i,j,k)
            enddo
         enddo
      enddo

      call pmaxmin(qname, tmp, pmin, pmax, im*(jlast-jfirst+1),km, 1.)

      return

      end

