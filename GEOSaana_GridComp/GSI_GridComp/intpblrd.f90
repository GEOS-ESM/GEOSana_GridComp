module intpblrdmod
!$$$ module documentation block
!           .      .    .                                       .
! module:   intpblrdmod    module for intpblrd and its tangent linear intpblrd_tl
!   prgmmr:
!
! abstract: module for intpblrd and its tangent linear intpblrd_tl
!
! program history log:
!
!   2012-09-14  Syed RH Rizvi, NCAR/NESL/MMM/DAS  - implemented obs adjoint test  
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!
! subroutines included:
!   sub intpblrd
!
! variable definitions:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use m_obsNode, only: obsNode
use m_pblrdNode, only: pblrdNode
use m_pblrdNode, only: pblrdNode_typecast
use m_pblrdNode, only: pblrdNode_nextcast
use m_obsdiagNode, only: obsdiagNode_set
implicit none

PRIVATE
PUBLIC intpblrd

contains

subroutine intpblrd(pblrdhead,rval,sval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intpblrd      apply nonlin qc obs operator for conv. pblrd
!   prgmmr: zhu           org: np23                date: 2011-02-20
!
! abstract: apply observation operator and adjoint for conventional pblrd
!           observations with nonlinear qc operator
!
! program history log:
!
!   2012-09-14  Syed RH Rizvi, NCAR/NESL/MMM/DAS  - introduced ladtest_obs         
!   2014-12-03  derber  - modify so that use of obsdiags can be turned off
!
!   input argument list:
!     pblrdhead
!     spblrd    - increment in grid space
!     rpblrd
!
!   output argument list:
!     rpblrd    - results from observation operator (0 for no data)
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
  use constants, only: half,one,tiny_r_kind,cg_term
  use obsmod, only: lsaveobsens, l_do_adjoint,luse_obsdiag
  use qcmod, only: nlnqc_iter,varqc_iter
  use jfunc, only: jiter
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  use gsi_4dvar, only: ladtest_obs
  implicit none

! Declare passed variables
  class(obsNode  ),pointer, intent(in   ) :: pblrdhead
  type(gsi_bundle),         intent(in   ) :: sval
  type(gsi_bundle),         intent(inout) :: rval

! Declare local variables
  integer(i_kind) ier,istatus
  integer(i_kind) j1,j2,j3,j4
! real(r_kind) penalty
  real(r_kind) w1,w2,w3,w4
  real(r_kind) val
  real(r_kind) cg_pblrd,p0,grad,wnotgross,wgross,pg_pblrd
  real(r_kind),pointer,dimension(:) :: spblrd
  real(r_kind),pointer,dimension(:) :: rpblrd
  type(pblrdNode), pointer :: pblrdptr

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'pblrd',spblrd,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'pblrd',rpblrd,istatus);ier=istatus+ier
  if(ier/=0)return

  !pblrdptr => pblrdhead
  pblrdptr => pblrdNode_typecast(pblrdhead)
  do while (associated(pblrdptr))
     j1=pblrdptr%ij(1)
     j2=pblrdptr%ij(2)
     j3=pblrdptr%ij(3)
     j4=pblrdptr%ij(4)
     w1=pblrdptr%wij(1)
     w2=pblrdptr%wij(2)
     w3=pblrdptr%wij(3)
     w4=pblrdptr%wij(4)

!    Forward model
     val=w1*spblrd(j1)+w2*spblrd(j2)&
        +w3*spblrd(j3)+w4*spblrd(j4)

     if(luse_obsdiag)then
        if (lsaveobsens) then
           grad = val*pblrdptr%raterr2*pblrdptr%err2
           !-- pblrdptr%diags%obssen(jiter) = grad
           call obsdiagNode_set(pblrdptr%diags,jiter=jiter,obssen=grad)
        else
           !-- if (pblrdptr%luse) pblrdptr%diags%tldepart(jiter)=val
           if (pblrdptr%luse) call obsdiagNode_set(pblrdptr%diags,jiter=jiter,tldepart=val)
        endif
     endif

     if (l_do_adjoint) then
        if (.not. lsaveobsens) then
           if( .not. ladtest_obs)   val=val-pblrdptr%res

!          gradient of nonlinear operator
           if (nlnqc_iter .and. pblrdptr%pg > tiny_r_kind .and. &
                                pblrdptr%b  > tiny_r_kind) then
              pg_pblrd=pblrdptr%pg*varqc_iter
              cg_pblrd=cg_term/pblrdptr%b
              wnotgross= one-pg_pblrd
              wgross = pg_pblrd*cg_pblrd/wnotgross
              p0   = wgross/(wgross+exp(-half*pblrdptr%err2*val**2))
              val = val*(one-p0)
           endif
           if( ladtest_obs ) then
              grad = val
           else
              grad = val*pblrdptr%raterr2*pblrdptr%err2
           end if
        endif

!       Adjoint
        rpblrd(j1)=rpblrd(j1)+w1*grad
        rpblrd(j2)=rpblrd(j2)+w2*grad
        rpblrd(j3)=rpblrd(j3)+w3*grad
        rpblrd(j4)=rpblrd(j4)+w4*grad
     endif

     !pblrdptr => pblrdptr%llpoint
     pblrdptr => pblrdNode_nextcast(pblrdptr)

  end do

  return
end subroutine intpblrd

end module intpblrdmod
