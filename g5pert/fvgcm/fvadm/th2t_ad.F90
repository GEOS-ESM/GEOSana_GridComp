!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: th2t_ad --- Convert (vir.) pot. temp to (vir.) temp (AD-only)
!
!
! !INTERFACE:
!
subroutine th2t_ad ( delp, pt, pt_ad, pkz_ad )

! !USES:

use precision
use mod_comm, only : numcpu
use stepon, only : im=>imr, &
                   jm=>jnp, &
                   km=>nl, &
                   jfirst, &
                   jlast, &
                   ng_d, &
                   ks, &
                   akap=>cappa, &
                   ptop

implicit none

! !INPUT PARAMETERS:

real(kind=r8), intent(in)    :: delp(im,jfirst:jlast,km)
real(kind=r8), intent(in)    :: pt(im,jfirst-ng_d:jlast+ng_d,km)    ! virtual potential temperature

! !INPUT/OUTPUT PARAMETERS:

real(kind=r8), intent(inout) :: pt_ad(im,jfirst-ng_d:jlast+ng_d,km) ! on  input: (vir.) potential temp. perturbation
                                                                    ! on output: (vir.)           temp. perturbation
real(kind=r8), intent(inout) :: pkz_ad(im,jfirst:jlast,km)

! !DESCRIPTION: Convert input pt_ad (virtual) potential temperature to 
!               (virtual) temperature. Notice this routine is unsual in the
!               sense that it overwrites the input array of perturbation pt_ad;
!               it also leave pt (the ref state) untouched. 
!
! !REVISION HISTORY:
!
!  31May2007  Todling    Initial code
!  01Jun2007  YZhu/RT    Bug fix in pt_ad calculation
!
!EOP
!-------------------------------------------------------------------------

! local parameters
! ----------------
real(kind=r8) :: pe  (im,km+1,jfirst:jlast)
real(kind=r8) :: peln(im,km+1,jfirst:jlast)
real(kind=r8) :: pk  (im,jfirst:jlast,km+1)
real(kind=r8) :: pkz (im,jfirst:jlast,km)

real(kind=r8) :: peln_ad(im,km+1,jfirst:jlast)
real(kind=r8) ptopk,peakap

integer nx,js2g0,jn2g0,i,j,k

if ( (jlast-jfirst+1)/numcpu >= 4 ) then
      nx = 1
else
      nx = 4
endif

! Fill in pressure edge arrays
! ----------------------------
ptopk = ptop ** akap
!$omp parallel do private(i,j,k)
do j=jfirst,jlast
   do i=1,im
      pe(i,1,j) = ptop
    enddo
    do k=2,km+1
       do i=1,im
           pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
       enddo
    enddo
enddo

!$omp parallel do default(shared) private(i,j,k)
do j=jfirst,jlast
   do i=1,im
      pk(i,j,1) = ptopk
   enddo
   do k=2,km+1
      do i=1,im
          pk(i,j,k) = pe(i,k,j)**akap
      enddo
   enddo
enddo

! Calculate p**kappa and it's tangent linear correspondent
! --------------------------------------------------------
call pkez   ( nx, im, km, jfirst, jlast, pe, pk, akap, ks, peln, pkz, .false.)


! Convert (vir.) potential temperature to (vir.) temperature
! ----------------------------------------------------------
js2g0 = max(2,jfirst)
jn2g0 = min(jm-1,jlast)
do k = 1, km
   if ( jlast==jm ) then
      j = jm
      do i = 1, im
!          pkz_ad(i,j,k) =  pt(i,j,k)*pt_ad(i,j,k)
            pt_ad(i,j,k) = pkz(i,j,k)*pt_ad(i,j,k)
      enddo
   endif
   if ( jfirst==1 ) then
      j = 1
      do i = 1, im
!          pkz_ad(i,j,k) =   pt(i,j,k)*pt_ad(i,j,k)
            pt_ad(i,j,k) =  pkz(i,j,k)*pt_ad(i,j,k)
      enddo
   endif
   do j = js2g0, jn2g0
      do i = 1, im
!          pkz_ad(i,j,k) =   pt(i,j,k)*pt_ad(i,j,k)
            pt_ad(i,j,k) =  pkz(i,j,k)*pt_ad(i,j,k)
      enddo
   enddo
enddo


end subroutine th2t_ad
