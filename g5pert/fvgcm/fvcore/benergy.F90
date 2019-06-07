!-----------------------------------------------------------------------
!BOP
module benergy             ! FastOpt

! !USES:
      use precision

      implicit none

! Geometric arrays
      real (r8), allocatable :: sine(:), cosp(:), sinp(:), cose(:)

      real (r8) acap

contains

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: benergy_initialize --- Initialization for benergy
!
! !INTERFACE:

   subroutine benergy_initialize( im, jm )
   implicit none

! !INPUT PARAMETERS:
      integer im, jm                  ! Dimensions

! !DESCRIPTION:
!    Allocates arrays and sets variables
!EOP
!---------------------------------------------------------------------
!BOC
! Local
      real (r8) dp, dl

   allocate( sine(jm), cosp(jm), sinp(jm), cose(jm) )
   call setrig(im, jm, dp, dl, cosp, cose, sinp, sine)

! WS 99.05.25 : changed real conversion of IMR to IM
! scaled polar cap area.
   acap = im*(1.+sine(2)) / dp

   deallocate( cose, sine, sinp )

!EOC
   end subroutine benergy_initialize
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!BOP
! !ROUTINE: benergy_finalize --- Finalize benergy
   subroutine benergy_finalize
   implicit none

! !DESCRIPTION:
!    Deallocates arrays
!EOP
!---------------------------------------------------------------------
!BOC

   deallocate( cosp )

!EOC
   end subroutine benergy_finalize
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
!BOP
! !ROUTINE: benergy --- Calculate the total energy
!
      subroutine benergy_do(im, jm, km, u, v, pt, delp, pe, pk, pkz, phis,  &
                            ng_d, ng_s, cp, te0, te, dz, jfirst, jlast)
! !USES:

      use precision

#if defined( SPMD )
      use mod_comm, only : mp_send_s, mp_recv_n, mp_barrier
#endif
      implicit none

! !INPUT PARAMETERS:
      integer im, jm, km, jfirst, jlast          ! Dimensions
      integer ng_d, ng_s
      real(r8) cp

! Ghosted prog arrays:
      real(r8) u(im,jfirst-ng_d:jlast+ng_s,km)  ! Winds x
      real(r8) v(im,jfirst-ng_s:jlast+ng_d,km)  ! Winds y
      real(r8) pt(im,jfirst-ng_d:jlast+ng_d,km) ! Potential temperature

      real(r8) delp(im,jfirst:jlast,km)         ! Delta pressure
      real(r8) pkz(im, jfirst:jlast, km)
      real(r8) pe(im, km+1, jfirst:jlast)           ! Edge pressure
      real(r8) pk(im, jfirst:jlast, km+1)
      real(r8) phis(im,jfirst:jlast)

! !INPUT work arrays:
      real(r8) te(im, jfirst:jlast, km)              ! Work array (cache perf.)
      real(r8) dz(im, jfirst:jlast, km)              ! Work array (cache perf.)

! !OUTPUT PARAMETERS:
      real(r8) te0              ! globally integrated total energy

! !DESCRIPTION:
!    Determines the globally integrated total energy, if jfirst == 1
!    and jlast == jnp, otherwise it calculates the total energy in
!    the latitude slice jfirst:jlast. 
!
! !REVISION HISTORY:
!
! SJL 99.04.13 : Delivered as release 0.9.8
! WS  99.05.18 : Added im, jm, km, te, dz as arguments
! WS  99.05.25 : Replaced IMR by IM, JMR by JM-1; removed fvcore.h
! WS  99.10.11 : Ghosted U, now fully limited to jfirst:jlast
! WS  99.11.23 : Pruned te, additional cleaning
! WS  00.05.14 : Renamed ghost indices as per Kevin's definitions
! WS  00.07.13 : Changed PILGRIM API
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local
      real(r8) u2(im,jfirst:jlast+1)
      real(r8) v2(im,jfirst:jlast)
      real(r8) bte(im)
      real(r8) tte(jfirst:jlast)
      real(r8) te_sp
      real(r8) te_np
      real(r8) gztop(im)

      real(r8) dp, dl, xsum

      integer i, j, k, js2g0, jn2g0

      js2g0  = max(2,jfirst)
      jn2g0  = min(jm-1,jlast)
 
#if defined( SPMD )
      call mp_send_s(im, jm, jfirst, jlast, 1, km, ng_d, ng_s, u)
      call mp_barrier
      call mp_recv_n(im, jm, jfirst, jlast, 1, km, ng_d, ng_s, u)
#endif

!$omp parallel do private(i, j, k, u2, v2, te_sp, te_np)
      do 1000 k=1,km
      do j=js2g0,min(jlast+1,jm)
         do i=1,im
            u2(i,j) = u(i,j,k)**2
         enddo
      enddo

      do j=js2g0,jn2g0
         do i=1,im
            v2(i,j) = v(i,j,k)**2
         enddo
      enddo

      do j=js2g0,jn2g0
         do i=1,im-1
         te(i,j,k) = 0.25*(u2(i,j) + u2(i,j+1) + v2(i,j) + v2(i+1,j))
         enddo
         te(im,j,k) = 0.25*(u2(im,j) + u2(im,j+1) + v2(im,j) + v2(1,j))
      enddo

      do j=js2g0,jn2g0
         do i=1,im
            te(i,j,k) = delp(i,j,k) * ( te(i,j,k) +   &
                           cp*pt(i,j,k)*pkz(i,j,k)  ) 
         enddo
      enddo

      if ( jfirst == 1 ) then
           te_sp = 0.
        do i=1,im
           te_sp = te_sp + u2(i,2) + v2(i,2)
        enddo

        te_sp =  delp(1,1,k) * ( 0.5*te_sp/float(im) +    &
                                 cp*pt(1,1,k)*pkz(1,1,k) )
        do i=1,im
           te(i,1,k) = te_sp
        enddo
      endif

      if ( jlast == jm ) then
           te_np = 0.
        do i=1,im
           te_np = te_np + u2(i,jm) + v2(i,jm-1)
        enddo
           te_np = delp(1,jm,k) * ( 0.5*te_np/float(im) +       &
                                    cp*pt(1,jm,k)*pkz(1,jm,k) )
        do i=1,im
           te(i,jm,k) = te_np
        enddo
      endif

      do j=jfirst,jlast
         do i=1,im
            dz(i,j,k) = cp*pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))
         enddo
      enddo
1000  continue


!$omp parallel do private(i,j,k,bte,xsum,gztop)
      do 2000 j=jfirst,jlast
! Perform vertical integration
         do i=1,im
               gztop(i) = phis(i,j)
            do k=1,km
               gztop(i) = gztop(i) + dz(i,j,k)
            enddo
         enddo

      if (j == 1) then
! SP
            tte(1) = pe(1,km+1,1)*phis(1,1) - pe(1,1,1)*gztop(1)
            do k=1,km
               tte(1) = tte(1) + te(1,1,k)
            enddo
            tte(1)  = acap * tte(1)

      elseif (j == jm) then
! NP
            tte(jm) = pe(1,km+1,jm)*phis(1,jm) - pe(1,1,jm)*gztop(1)
            do k=1,km
               tte(jm) = tte(jm) + te(1,jm,k)
            enddo
            tte(jm) = acap * tte(jm)

      else
! Interior

         do i=1, im
            bte(i) = pe(i,km+1,j)*phis(i,j) - pe(i,1,j)*gztop(i)
         enddo

         do k=1,km
            do i=1,im
               bte(i) = bte(i) + te(i,j,k)
            enddo
         enddo

            xsum = 0.
         do i=1,im
            xsum = xsum + bte(i)
         enddo
         tte(j) = xsum*cosp(j)

      endif
2000  continue

      call par_vecsum(jm, jfirst, jlast, tte, te0)

      return
!EOC
      end subroutine benergy_do
!-----------------------------------------------------------------------


end module benergy             ! FastOpt
!-----------------------------------------------------------------------
