!-----------------------------------------------------------------------
!BOP
! !ROUTINE: drymadj --- Total dry air mass
!
! !INTERFACE:

      subroutine drymadj( im, jm, km, jfirst, jlast, ng,          &
                          moun, ptop, ps, delp, pe, nq, q, id )

      use precision
#if defined( SPMD )
#define CPP_PRT_PREFIX  if(gid.eq.0)
      use mod_comm, only : gid
#else
#define CPP_PRT_PREFIX
#endif
      use gmean, only: gmean_of

      implicit   none

! !INPUT PARAMETERS:
      integer im, jm, km     ! Dimensions
      integer jfirst, jlast  ! Latitude strip
      integer   id           ! 0:  checking total dry mass only
                             ! 1:  checking total dry mass and adjust
      integer nq             ! Number of tracers         
      integer ng
      logical moun
      real(r8)  ptop

! !INPUT/OUTPUT PARAMETERS:
      real(r8) delp(im,jfirst:jlast,km)     !
      real(r8) pe(im,km+1,jfirst:jlast)     !
      real(r8)   q(im,jfirst-ng:jlast+ng,km,nq) 
      real(r8)  ps(im,jfirst:jlast)        ! surface pressure

! !DESCRIPTION:
!  Perform adjustment of the total dry-air-mass while preserving total
!  tracer mass
!  Developer: S.-J. Lin
!
! !REVISION HISTORY:
!   WS  99.10.26:  Revision; documentation; removed fvcore.h
!   WS  99.11.02:  Limited to jfirst:jlast
!   SJL 00.03.20:
!   WS  00.07.08:  Use precision module
!   WS  00.07.10:  Incorporate simplifications to PILGRIM
!
!EOP
!---------------------------------------------------------------------
!BOC

      real(r8) psd(im,jfirst:jlast)     ! surface pressure  due to dry air mass

      real(r8) drym                     ! dry air mass in pascals
!     parameter ( drym = 98222. )       ! setting for US NAVY 10 min data
      parameter ( drym = 98288. )       ! setting for USGS

      integer   i, j, k
      real(r8)  psmo
      real(r8)  psdry
      real(r8)  dpd
      integer ic

!$omp  parallel do          &
!$omp  default(shared)      &
!$omp  private(i,j,k)

      do 1000 j=jfirst,jlast
        do i=1,im
          psd(i,j) = ptop
        enddo

        if(nq .ne. 0) then
          do k=1, km
            do i=1,im
              psd(i,j) = psd(i,j) + delp(i,j,k)*(1.-q(i,j,k,1))
            enddo
          enddo
        else
          do k=1, km
            do i=1,im
              psd(i,j) = psd(i,j) +  delp(i,j,k)
            enddo
          enddo
        endif
1000  continue

! Check global maximum/minimum
      psmo = gmean_of(im, jm, jfirst, jlast, ps(1,jfirst) )
      CPP_PRT_PREFIX write(6,*)                                      &
              'Total (moist) surface pressure before adjustment = ', &
               0.01*psmo
      psdry = gmean_of(im, jm, jfirst, jlast, psd(1,jfirst))
      CPP_PRT_PREFIX write(6,*) 'mean (dry) surface pressure = ', 0.01*psdry
      CPP_PRT_PREFIX write(6,*) 'TPW (kg/m**2) =', (psmo-psdry)/9.80616

      if( id .eq. 0) return

      if(moun) then
        dpd = drym - psdry
      else
        dpd = 1000.*100. - psdry
      endif

      CPP_PRT_PREFIX write(6,*) 'dry mass to be added (pascals) =', dpd

!$omp  parallel do            &
!$omp  default(shared)        &
!$omp  private(i,j, ic)

      do 2000 j=jfirst,jlast

         do ic=1,nq
            do i=1,im
               q(i,j,km,ic) = q(i,j,km,ic)*delp(i,j,km) /     &
                             (delp(i,j,km)+dpd)
            enddo
         enddo

! Adjust the lowest Lagrangian layer
         do i=1,im
            delp(i,j,km) = delp(i,j,km) + dpd
            pe(i,km+1,j) = pe(i,km,j) + delp(i,j,km)
            ps(i,j) = pe(i,km+1,j)
         enddo
2000  continue

      psmo = gmean_of(im, jm, jfirst, jlast, ps(1,jfirst) )
      CPP_PRT_PREFIX write(6,*)                                     &
              'Total (moist) surface pressure after adjustment = ', &
               0.01*psmo

      return
!EOC
      end
!---------------------------------------------------------------------
