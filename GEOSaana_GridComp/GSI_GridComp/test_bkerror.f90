subroutine test_bkerror

use constants, only: zero,one
use kinds, only: r_kind
use mpimod, only: mype
use gsi_4dvar, only: nsubwin
use gridmod, only: nlon,nlat,lon2,lat2
use state_vectors, only: allocate_state,deallocate_state
use control_vectors, only: control_vector
use control_vectors, only: allocate_cv,deallocate_cv
use control_vectors, only: dot_product,assignment(=)
use bias_predictors, only: predictors,allocate_preds,deallocate_preds,assignment(=)
use gsi_bundlemod, only: gsi_bundle
use gsi_bundlemod, only: gsi_bundlegetpointer
use gsi_bundlemod, only: assignment(=)
use gsi_4dcouplermod, only: gsi_4dcoupler_getpert_set
use gsi_4dcouplermod, only: gsi_4dcoupler_getpert
use gsi_4dcouplermod, only: gsi_4dcoupler_putpert_set
use gsi_4dcouplermod, only: gsi_4dcoupler_putpert
use gsi_4dcouplermod, only: gsi_4dcoupler_putpert_final
use hybrid_ensemble_parameters,only : l_hyb_ens,ntlevs_ens
use hybrid_ensemble_isotropic, only: bkerror_a_en
use obsmod, only: iadate

implicit none

type(gsi_bundle), allocatable :: fcgrad(:)
type(gsi_bundle) :: aux
type(control_vector) :: gradx,grady
type(predictors)     :: sbias
type(gsi_bundle),allocatable :: eval(:)

character(len=80):: ifname(nsubwin)
character(len=80):: ofname
integer :: ii, mydate(5),ier
integer :: nymd,nhms
logical :: dotest,fexist
real(r_kind), pointer :: ptr3d(:,:,:)

ifname = 'xinc.eta.nc4'
ofname = 'outvec_bkgtest'
inquire(file=ifname(1),exist=fexist)
if(.not.fexist) return
 
! start work space
allocate(fcgrad(nsubwin))
do ii=1,nsubwin
   call allocate_state(fcgrad(ii))
end do
call allocate_preds(sbias)
allocate(eval(ntlevs_ens))
do ii=1,ntlevs_ens
   call allocate_state(eval(ii))
end do

call allocate_cv(gradx)
call allocate_cv(grady)
gradx=zero
grady=zero
  
! get test vector (fcgrad)
call gsi_4dcoupler_getpert_set(nlat,nlon,lat2,lon2,nsubwin)
call gsi_4dcoupler_getpert(fcgrad,nsubwin,'inc',ifname)

if (l_hyb_ens) then
   eval(1)=fcgrad(1)
   fcgrad(1)=zero
   call ensctl2state_ad(eval,fcgrad(1),gradx)
endif
call control2state_ad(fcgrad,sbias,gradx)

! apply B to input (transformed) vector
call bkerror(gradx,grady)
if(l_hyb_ens) then
  call bkerror_a_en(gradx,grady)
endif

call control2state(grady,fcgrad,sbias)
if (l_hyb_ens) then
   call ensctl2state(grady,fcgrad(1),eval)
   fcgrad(1)=eval(1)
endif

! output column of covariance
nymd = iadate(1)*10000 + iadate(2)*100 + iadate(3)
nhms = iadate(4)*10000 + iadate(5)*100
call gsi_4dcoupler_putpert_set (nymd,nhms,ier)
call gsi_4dcoupler_putpert(fcgrad(1),nymd,nhms,'tlm',ofname)
call gsi_4dcoupler_putpert_final (ier)

! start work space
call deallocate_cv(gradx)
call deallocate_cv(grady)
do ii=1,ntlevs_ens
   call deallocate_state(eval(ii))
end do
deallocate(eval)
call deallocate_preds(sbias)
do ii=1,nsubwin
   call deallocate_state(fcgrad(ii))
end do
deallocate(fcgrad)

end subroutine test_bkerror
