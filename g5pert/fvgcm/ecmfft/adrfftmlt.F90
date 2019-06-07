!-----------------------------------------------------------------------
!BOP
! !IROUTINE: adrfftmlt --- adjoint fast Fourier transform
!
! !INTERFACE: 
 subroutine adrfftmlt(ada,work,trigs,ifax,inc,jump,n,lot,isign)

! !USES:
 use precision

! !DESCRIPTION:
!   The adjoint Fourier transform is the reverse Fourier transform
!   with appropriate normalization.
!
! !REVISION HISTORY:
!
!   01.09.25   Ralf Giering, FastOpt     Initial version
!
!EOP
!-----------------------------------------------------------------------
!BOC
      real(r8) ada(jump*lot)
      real(r8) work((n+1)*lot)
      real(r8) trigs(3*n/2+1)
      integer ifax(13)

      call rfftmlt(ada,work,trigs,ifax,inc,jump,n,lot,-isign)

!EOC
 end subroutine adrfftmlt
!-----------------------------------------------------------------------

