#include <misc.h>
#include <preproc.h>

subroutine v2xy (fldv, fldxyini, fldxy)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Perfrom grid-average from subgrid patch vector
! 
! Method: 
! Subgrid patch to grid average mapping: average a subgrid input vector 
! [fldv] of length [numpatch] to a 2-d [lsmlon] x [lsmlat] output array [fldxy]
! setting non-land points to [nonland]. Averaging is only done for points 
! that are not equal to "spval".
!
! Author: Gordon Bonan
! 
!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------

  use precision
  use clm_varder 
  use clm_varsur, only : landmask
  use clm_varpar, only : lsmlon, lsmlat
  use clm_varmap, only : numpatch, patchvec
  use clm_varcon, only : spval
  implicit none

! ------------------------ arguments----------------------------------
  real(r8), intent(in)  :: fldv(numpatch)       !subgrid vector input
  real(r8), intent(in)  :: fldxyini             !initial value of fldxy
  real(r8), intent(out) :: fldxy(lsmlon,lsmlat) !gridded output
! --------------------------------------------------------------------

! ------------------------ local variables ----------------------
  integer :: i                    !longitude index     
  integer :: j                    !latitude index      
  integer :: k                    !subgrid patch index 
  real(r8):: sumwt(lsmlon,lsmlat) !sum of wt
! ---------------------------------------------------------------

! Loop over subgrid patches to create grid average. 

  fldxy(:,:) = fldxyini
  sumwt(:,:) = 0.
  do k = 1, numpatch
     if (fldv(k) /= spval) then
        i = patchvec%ixy(k)      !longitude index for land point
        j = patchvec%jxy(k)      !latitude index for land point
        if (sumwt(i,j)==0.) fldxy(i,j) = 0.
        fldxy(i,j) = fldxy(i,j) + patchvec%wtxy(k)*fldv(k)
        sumwt(i,j) = sumwt(i,j) + patchvec%wtxy(k)
     endif
  end do
  where (landmask(:,:) == 1 .and. sumwt(:,:) /= 0.)
     fldxy(:,:) = fldxy(:,:)/sumwt(:,:)
  endwhere

  return
end subroutine v2xy
