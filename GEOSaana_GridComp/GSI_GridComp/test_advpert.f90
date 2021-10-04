subroutine test_advpert

use constants, only: zero
use mpimod, only: mype
use gsi_4dvar, only: nsubwin
use gridmod, only: nlon,nlat,lon2,lat2
use state_vectors, only: allocate_state,deallocate_state
use gsi_bundlemod, only: gsi_bundle
use gsi_4dcouplermod, only: gsi_4dcoupler_getpert_set
use gsi_4dcouplermod, only: gsi_4dcoupler_getpert
use gsi_4dcouplermod, only: gsi_4dcoupler_putpert_set
use gsi_4dcouplermod, only: gsi_4dcoupler_putpert
use gsi_4dcouplermod, only: gsi_4dcoupler_putpert_final

implicit none

type(gsi_bundle), allocatable :: fcgrad(:)

character(len=80):: ifname(nsubwin)
character(len=80):: ofname
integer, parameter :: tau_fcst = 24 ! hours
integer :: ii, mydate(5),ier
logical :: dotest,fexist

ifname = 'blob.nc4'
ofname = 'advected_blob.nc4'
inquire(file=ifname(1),exist=fexist)
if(.not.fexist) return
 
! start work space
allocate(fcgrad(nsubwin))
do ii=1,nsubwin
   call allocate_state(fcgrad(ii))
end do
  
! get test vector (fcgrad)
call gsi_4dcoupler_getpert_set(nlat,nlon,lat2,lon2,nsubwin)
call gsi_4dcoupler_getpert(fcgrad,nsubwin,'adm',ifname)

! propagate field
mydate(1)=1776;mydate(2)=7;mydate(3)=4;mydate(4)=12;mydate(5)=0
call advect_cv(mype,mydate,tau_fcst,fcgrad)

! after advection
call gsi_4dcoupler_putpert_set (17760704,120000,ier)
call gsi_4dcoupler_putpert(fcgrad(1),17760704,120000,'tlm',ofname)
call gsi_4dcoupler_putpert_final (ier)

! start work space
do ii=1,nsubwin
   call deallocate_state(fcgrad(ii))
end do
deallocate(fcgrad)

end subroutine test_advpert
