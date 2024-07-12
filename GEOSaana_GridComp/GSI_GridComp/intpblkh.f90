module intpblkhmod
!$$$ module documentation block
!           .      .    .                                       .
! module:   intpblkhmod    module for intpblkh and its tangent linear intpblkh_tl
!   prgmmr:
!
! abstract: module for intpblkh and its tangent linear intpblkh_tl
!
! program history log:
!
!   2012-09-14  Syed RH Rizvi, NCAR/NESL/MMM/DAS  - implemented obs adjoint test  
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!
! subroutines included:
!   sub intpblkh
!
! variable definitions:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use m_obsNode, only: obsNode
use m_pblkhNode, only: pblkhNode
use m_pblkhNode, only: pblkhNode_typecast
use m_pblkhNode, only: pblkhNode_nextcast
use m_obsdiagNode, only: obsdiagNode_set
implicit none

PRIVATE
PUBLIC intpblkh

contains

subroutine intpblkh(pblkhhead,rval,sval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intpblkh      apply nonlin qc obs operator for conv. pblkh
!   prgmmr: zhu           org: np23                date: 2011-02-20
!
! abstract: apply observation operator and adjoint for conventional pblkh
!           observations with nonlinear qc operator
!
! program history log:
!
!   2012-09-14  Syed RH Rizvi, NCAR/NESL/MMM/DAS  - introduced ladtest_obs         
!   2014-12-03  derber  - modify so that use of obsdiags can be turned off
!
!   input argument list:
!     pblkhhead
!     spblkh    - increment in grid space
!     rpblkh
!
!   output argument list:
!     rpblkh    - results from observation operator (0 for no data)
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
  class(obsNode  ),pointer, intent(in   ) :: pblkhhead
  type(gsi_bundle),         intent(in   ) :: sval
  type(gsi_bundle),         intent(inout) :: rval

! Declare local variables
  integer(i_kind) ier,istatus
  integer(i_kind) j1,j2,j3,j4
! real(r_kind) penalty
  real(r_kind) w1,w2,w3,w4
  real(r_kind) val
  real(r_kind) cg_pblkh,p0,grad,wnotgross,wgross,pg_pblkh
  real(r_kind),pointer,dimension(:) :: spblkh
  real(r_kind),pointer,dimension(:) :: rpblkh
  type(pblkhNode), pointer :: pblkhptr

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'pblkh',spblkh,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'pblkh',rpblkh,istatus);ier=istatus+ier
  if(ier/=0)return

  !pblkhptr => pblkhhead
  pblkhptr => pblkhNode_typecast(pblkhhead)
  do while (associated(pblkhptr))
     j1=pblkhptr%ij(1)
     j2=pblkhptr%ij(2)
     j3=pblkhptr%ij(3)
     j4=pblkhptr%ij(4)
     w1=pblkhptr%wij(1)
     w2=pblkhptr%wij(2)
     w3=pblkhptr%wij(3)
     w4=pblkhptr%wij(4)

!    Forward model
     val=w1*spblkh(j1)+w2*spblkh(j2)&
        +w3*spblkh(j3)+w4*spblkh(j4)

     if(luse_obsdiag)then
        if (lsaveobsens) then
           grad = val*pblkhptr%raterr2*pblkhptr%err2
           !-- pblkhptr%diags%obssen(jiter) = grad
           call obsdiagNode_set(pblkhptr%diags,jiter=jiter,obssen=grad)
        else
           !-- if (pblkhptr%luse) pblkhptr%diags%tldepart(jiter)=val
           if (pblkhptr%luse) call obsdiagNode_set(pblkhptr%diags,jiter=jiter,tldepart=val)
        endif
     endif

     if (l_do_adjoint) then
        if (.not. lsaveobsens) then
           if( .not. ladtest_obs)   val=val-pblkhptr%res

!          gradient of nonlinear operator
           if (nlnqc_iter .and. pblkhptr%pg > tiny_r_kind .and. &
                                pblkhptr%b  > tiny_r_kind) then
              pg_pblkh=pblkhptr%pg*varqc_iter
              cg_pblkh=cg_term/pblkhptr%b
              wnotgross= one-pg_pblkh
              wgross = pg_pblkh*cg_pblkh/wnotgross
              p0   = wgross/(wgross+exp(-half*pblkhptr%err2*val**2))
              val = val*(one-p0)
           endif
           if( ladtest_obs ) then
              grad = val
           else
              grad = val*pblkhptr%raterr2*pblkhptr%err2
           end if
        endif

!       Adjoint
        rpblkh(j1)=rpblkh(j1)+w1*grad
        rpblkh(j2)=rpblkh(j2)+w2*grad
        rpblkh(j3)=rpblkh(j3)+w3*grad
        rpblkh(j4)=rpblkh(j4)+w4*grad
     endif

     !pblkhptr => pblkhptr%llpoint
     pblkhptr => pblkhNode_nextcast(pblkhptr)

  end do

  return
end subroutine intpblkh

end module intpblkhmod
