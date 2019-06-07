!-----------------------------------------------------------------------
!BOP
! !ROUTINE: trac2d --- Remap Lagrangian to fixed coordinates
!
! !INTERFACE:
      subroutine trac2d(dp1, q, nq, cx, cy, mfx, mfy, iord, jord,       &
                        ng, sine, cosp, acosp, acap, rcap, fill,        &
                        im, jm, km, jfirst, jlast, va, flx,             &
			trac2d_tape_rec                                 &
			)
 
! !USES:
      use precision
      use tp_core
      use fill_module

#if defined( SPMD )
      use mod_comm, only: mp_send4d_ns, mp_recv4d_ns,                   &
                          mp_send3d_ns, mp_recv3d_ns,                   &
                          mp_send2_s, mp_recv2_n, mp_reduce_max, mp_barrier
#endif
      implicit none

! !INPUT PARAMETERS:

      integer im, jm, km, jfirst, jlast
      integer ng                      ! Max number of ghost latitudes

      real(r8) dp1(im,jfirst:jlast,km)
      real(r8)  cx(im,jfirst-ng:jlast+ng,km)
      real(r8)  cy(im,jfirst:jlast+1,km)
      real(r8) mfx(im,jfirst:jlast,km)
      real(r8) mfy(im,jfirst:jlast+1,km)

      real(r8)  sine(jm)
      real(r8)  cosp(jm)
      real(r8) acosp(jm)

      integer nq
      integer iord,  jord
      logical fill
      real(r8) acap, rcap
      integer trac2d_tape_rec

! !INPUT/OUTPUT PARAMETERS:
      real(r8) q(im,jfirst-ng:jlast+ng,km,nq)

! Input work arrays
      real(r8)  va(im,jfirst:jlast,km)
      real(r8) flx(im,jfirst:jlast,km)

! !DESCRIPTION:
!
!  Perform large-time-step tracer transport using accumulated Courant
!  numbers (cx, cy) and the mass fluxes (mfx, mfy) within the Lagrangian
!  layers.  This routine is 100\% parallel in the vertical direction
!  (with SMP).  Merdional Courant number will be further split, if
!  necessary, to ensure stability.  Cy <= 1 away from poles; Cy $\le$
!  1/2 at the latitudes closest to the poles.
!
! !CALLED FROM:
!     fvcore
!
! !REVISION HISTORY:
!
!   SJL 99.04.13:  Delivery
!   WS  99.05.26:  Added jfirst:jlast concept; im, jm, km as parameters
!                  replaced IMR, JMR, JNP, NL with IM, JM-1, JM and KM
!   WS  99.09.27:  Documentation; indentation; jfirst:jlast 
!   WS  99.09.30:  Ghosting; loop limits; full parallelization; tested
!   SJL 99.10.15:  nsplt migrated to outermost loop to remove bug
!   SJL 99.12.19:  Local 2D arrays trimmed!
!   WS  00.05.14:  Renamed ghost indices as per Kevin's definitions
!   WS  00.07.13:  Changed PILGRIM API
!
!EOP
!---------------------------------------------------------------------
!BOC
      real(r8)  tiny
      parameter ( tiny = 1.e-10 )

! Local variables:
      real(r8) cymax(km)
      real(r8) cy_global
      real(r8) frac
      real(r8) cmax
      real(r8) sum1, sum2
      real(r8) dp2(im,jfirst:jlast,km)

! Local 2d arrays
      real(r8) a2(im,jfirst:jlast)
      real(r8) fx(im,jfirst:jlast)
      real(r8) fy(im,jfirst:jlast+1)
      logical ffsl(jm,km)
      integer i, j, k
      integer it, iq, nsplt
      integer js1gd, js2g0, js2gd, jn2g0,jn2gd,jn1g1,jn1gd

      integer, parameter :: nsplt_max = 5    ! FastOpt

#define TRAC2D_STORING
#ifdef  TRAC2D_STORING
!$taf init trac2d_tape    = static, 1
!$taf init trac2d_tape_ns = static, nsplt_max
      trac2d_tape_rec = 0
#endif

#if defined( SPMD )
      call mp_barrier
      call mp_send3d_ns(im, jm, jfirst, jlast, 1, km, ng, ng, cx, 1) 
      call mp_send2_s  (im, jm, jfirst, jlast, 1, km,  0,  1, cy, mfy)
#endif
      js2g0 = max(2,jfirst)
      jn2g0 = min(jm-1,jlast)
      jn1g1 = min(jm,jlast+1)
      js1gd = max(1,jfirst-ng)     ! NG latitudes on S (starting at 1)
      js2gd = max(2,jfirst-ng)     ! NG latitudes on S (starting at 2)
      jn2gd = min(jm-1,jlast+ng)   ! NG latitudes on S (ending at jm-1)
      jn1gd = min(jm,jlast+ng)     ! NG latitudes on N (ending at jm)

!$omp parallel do private(i,j,k,cmax)
      do k=1,km
         cymax(k) = 0.
         do j=js2g0,jlast
            cmax = 0.
            do i=1,im
               cmax = max( abs(cy(i,j,k)), cmax)
            enddo
            cymax(k) = max(cymax(k), cmax*(1. + sine(j)**16) )
         enddo
      enddo

#if defined( SPMD )
      call mp_barrier
      call mp_recv3d_ns(im, jm, jfirst, jlast, 1, km, ng, ng, cx, 1) 
      call mp_recv2_n  (im, jm, jfirst, jlast, 1, km,  0,  1, cy, mfy)
      call mp_reduce_max(km, cymax)
#endif

!$taf store cx,cy,mfy,cymax = trac2d_tape, rec=trac2d_tape_rec+1

#if defined( SPMD )
      call mp_barrier
      call mp_send4d_ns(im, jm, jfirst, jlast, 1, km, nq, ng, ng, q) 
#endif


! find global max cymax
            cy_global = cymax(1)
      if ( km /= 1 ) then                ! if NOT shallow water test case
         do k=2,km
            cy_global = max(cymax(k), cy_global)
         enddo
      endif

      nsplt = int(1. + cy_global)
      frac  = 1. / float(nsplt)

!$omp parallel do private(i,j,k)
       do 4000 k=1,km

        if( nsplt /= 1 ) then
!!!     do j=2,jm-1
        do j=js2gd,jn2gd                  
          do i=1,im
            cx(i,j,k) =  cx(i,j,k) * frac      ! cx ghosted on N*ng S*ng
          enddo
        enddo

!!!     do j=2,jm-1
        do j=js2g0,jn2g0
          do i=1,im
            mfx(i,j,k) = mfx(i,j,k) * frac
          enddo
        enddo

!!!     do j=2,jm
        do j=js2g0,jn1g1                     
          do i=1,im
             cy(i,j,k) =  cy(i,j,k) * frac    ! cy ghosted on N
            mfy(i,j,k) = mfy(i,j,k) * frac    ! mfy ghosted on N
          enddo
        enddo
        endif

!!!     do j=2,jm-1
        do j=js2g0,jn2g0
          do i=1,im
             if(cy(i,j,k)*cy(i,j+1,k) > 0.) then
                if( cy(i,j,k) > 0.) then
                   va(i,j,k) = cy(i,j,k)
                else
                   va(i,j,k) = cy(i,j+1,k)      ! cy ghosted on N
                endif
             else
              va(i,j,k) = 0.
             endif
          enddo
        enddo

! Check if FFSL extension is needed.

!!!        do 2222 j=2,jm-1
        do 2222 j=js2gd,jn2gd             ! flux needed on N*ng S*ng
          ffsl(j,k) = .false.
          do i=1,im
            if(abs(cx(i,j,k)) > 1.) then  ! cx ghosted on N*ng S*ng
              ffsl(j,k) = .true.
              go to 2222
            endif
          enddo
2222    continue

! Scale e-w mass fluxes
!!!        do j=2,jm-1
        do j=js2g0,jn2g0
          if( ffsl(j,k) ) then
            do i=1,im
              flx(i,j,k) = mfx(i,j,k)/sign(max(abs(cx(i,j,k)),tiny),cx(i,j,k))
            enddo
          else
            do i=1,im
              flx(i,j,k) = mfx(i,j,k)
            enddo
          endif
        enddo
4000  continue

!$taf store flx = trac2d_tape, rec=trac2d_tape_rec+1

      do 6000 it=1, nsplt

#if defined( SPMD )
      if ( it /= 1 ) then
      call mp_send4d_ns(im, jm, jfirst, jlast, 1, km, nq, ng, ng, q) 
      endif
#endif

!$taf store dp1 = trac2d_tape_ns, rec=it+trac2d_tape_rec*nsplt

!$omp parallel do private(i, j, k, sum1, sum2)
      do 3000 k=1,km
          do j=js2g0,jn2g0
            do i=1,im-1
              dp2(i,j,k) =  dp1(i,j,k) + mfx(i,j,k) - mfx(i+1,j,k) +  &
                           (mfy(i,j,k) - mfy(i,j+1,k)) * acosp(j)
           enddo
           dp2(im,j,k) = dp1(im,j,k) + mfx(im,j,k) - mfx(1,j,k) +     &
                        (mfy(im,j,k) - mfy(im,j+1,k)) * acosp(j)
          enddo

! Poles
          if ( jfirst == 1  ) then
!
!           sum1 = SUM( mfy(1:im, 2,k) )
!
            sum1 = 0.
            do i=1,im
               sum1 = sum1 + mfy(i,2,k)
            enddo

            sum1 = - sum1 * rcap
            do i=1,im
              dp2(i,1,k) = dp1(i, 1,k) +  sum1
            enddo
          endif

          if ( jlast  == jm ) then
!           sum2 = SUM( mfy(1:im,jm,k) )
            sum2 = 0.
            do i=1,im
               sum2 = sum2 + mfy(i,jm,k)
            enddo

            sum2 = sum2 * rcap
            do i=1,im
              dp2(i,jm,k) = dp1(i,jm,k) +  sum2
            enddo
          endif
3000  continue

#if defined( SPMD )
      call mp_barrier
      call mp_recv4d_ns(im, jm, jfirst, jlast, 1, km, nq, ng, ng, q) 
#endif

!$taf store q = trac2d_tape_ns, rec=it+trac2d_tape_rec*nsplt

!$omp parallel do private(i, j, k, iq, fx, fy, a2)
      do 5000 k=1,km
!$taf loop = parallel
         do 555 iq=1,nq
            call tp2c(a2, va(1,jfirst,k), q(1,jfirst-ng,k,iq),      &
                      cx(1,jfirst-ng,k) , cy(1,jfirst,k),           &
                      im, jm, iord, jord, ng,                       &
                      fx, fy, ffsl(1,k), rcap, acosp,               &
                      flx(1,jfirst,k), mfy(1,jfirst,k),             & 
                      cosp, 1, jfirst, jlast )

            do j=jfirst,jlast
              do i=1,im
                q(i,j,k,iq) = q(i,j,k,iq)*dp1(i,j,k) + a2(i,j)
              enddo
            enddo

         if (fill) call fillxy (q(1,jfirst,k,iq),  im, jm, jfirst, jlast,   &
                                acap, cosp, acosp)

            do j=jfirst,jlast
              do i=1,im
                q(i,j,k,iq) = q(i,j,k,iq) / dp2(i,j,k)
              enddo
            enddo
555       continue
      
          if( it /= nsplt ) then
            do j=jfirst,jlast
              do i=1,im
                dp1(i,j,k) = dp2(i,j,k)
              enddo
            enddo
          endif
5000  continue
6000  continue

      return
!EOC
      end
!-----------------------------------------------------------------------
