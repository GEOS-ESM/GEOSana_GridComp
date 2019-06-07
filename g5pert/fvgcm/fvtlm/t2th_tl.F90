!---------------------------------------------------------------------------
! NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3, GEOS/DAS !
!---------------------------------------------------------------------------
!BOP
!
! !ROUTINE: t2th_tl --- Convert (vir.) pot. temp to (vir.) temp (TL-only)
!
!
! !INTERFACE:
!
subroutine t2th_tl ( delp, delp_tl, pt, pt_tl )

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

real(kind=r8), intent(in)    :: pt(im,jfirst-ng_d:jlast+ng_d,km)    ! virtual potential temperature
real(kind=r8), intent(in)    :: delp(im,jfirst:jlast,km)
real(kind=r8), intent(in)    :: delp_tl(im,jfirst:jlast,km)

! !OUTPUT PARAMETERS:

real(kind=r8), intent(inout) :: pt_tl(im,jfirst-ng_d:jlast+ng_d,km) ! on  input: virtual pot. temp. perturbation
                                                                    ! on output: virtual      temp. perturbation

! !DESCRIPTION: Convert input pt_tl (virtual) potential temperature to 
!               (virtual) temperature. Notice this routine is unsual in the
!               sense that it overwrites the input array of perturbation pt_tl;
!               it also leave pt (the ref state) untouched. 
!
! !REMARKS: For whatever reason the ADM only verifies when I set pkz_tl
!           zero, here and in the reverse transformation. This needs
!           revision and correction.
!
! !REVISION HISTORY:
!
!  31May2007  Todling    Initial code
!
!EOP
!-------------------------------------------------------------------------

! local parameters
! ----------------
real(kind=r8) :: pe     (im,km+1,jfirst:jlast)
real(kind=r8) :: pe_tl  (im,km+1,jfirst:jlast)
real(kind=r8) :: peln   (im,km+1,jfirst:jlast)
real(kind=r8) :: peln_tl(im,km+1,jfirst:jlast)
real(kind=r8) :: pk     (im,jfirst:jlast,km+1)
real(kind=r8) :: pk_tl  (im,jfirst:jlast,km+1)
real(kind=r8) :: pkz    (im,jfirst:jlast,km)
real(kind=r8) :: pkz_tl (im,jfirst:jlast,km)

integer nx,js2g0,jn2g0,i,j,k
real(kind=r8) ptopk,peakap

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
      pe_tl(i,1,j) = 0.0
      pe(i,1,j) = ptop
    enddo
    do k=2,km+1
       do i=1,im
           pe_tl(i,k,j) = pe_tl(i,k-1,j) + delp_tl(i,j,k-1)
              pe(i,k,j) =    pe(i,k-1,j) +    delp(i,j,k-1)
       enddo
    enddo
enddo

!$omp parallel do default(shared) private(i,j,k)
do j=jfirst,jlast
   do i=1,im
      pk_tl(i,j,1) = 0.0
      pk(i,j,1) = ptopk
   enddo
   do k=2,km+1
      do i=1,im
          peakap = pe(i,k,j)**akap
          pk_tl(i,j,k) = akap*peakap * pe_tl(i,k,j) / pe(i,k,j)
             pk(i,j,k) =      peakap
      enddo
   enddo
enddo


! Calculate p**kappa and it's tangent linear correspondent
! --------------------------------------------------------
call pkez_tl( nx, im, km, jfirst, jlast, pe, pe_tl, pk, pk_tl, akap, ks, peln, peln_tl, pkz, pkz_tl, .false. )
pkz_tl = 0.

! Convert (vir.) potential temperature to (vir.) temperature
! ----------------------------------------------------------
js2g0 = max(2,jfirst)
jn2g0 = min(jm-1,jlast)
do k = 1, km
   do j = js2g0, jn2g0
      do i = 1, im
         pt_tl(i,j,k) = pt_tl(i,j,k)/pkz(i,j,k) - pkz_tl(i,j,k)*pt(i,j,k)/pkz(i,j,k)
      enddo
   enddo
   if ( jfirst==1 ) then
      j = 1
      do i = 1, im
         pt_tl(i,j,k) = pt_tl(i,j,k)/pkz(i,j,k) - pkz_tl(i,j,k)*pt(i,j,k)/pkz(i,j,k)
      enddo
   endif
   if ( jlast==jm ) then
      j = jm
      do i = 1, im
         pt_tl(i,j,k) = pt_tl(i,j,k)/pkz(i,j,k) - pkz_tl(i,j,k)*pt(i,j,k)/pkz(i,j,k)
      enddo
   endif
enddo

end subroutine t2th_tl
