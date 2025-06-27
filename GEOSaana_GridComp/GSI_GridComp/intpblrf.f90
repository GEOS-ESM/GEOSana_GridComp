module intpblrfmod
!$$$ module documentation block
!           .      .    .                                       .
! module:   intpblrfmod    module for intpblrf and its tangent linear intpblrf_tl
!   prgmmr:
!
! abstract: module for intpblrf and its tangent linear intpblrf_tl
!
! program history log:
!
!   2012-09-14  Syed RH Rizvi, NCAR/NESL/MMM/DAS  - implemented obs adjoint test  
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!
! subroutines included:
!   sub intpblrf
!
! variable definitions:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use m_obsNode, only: obsNode
use m_pblrfNode, only: pblrfNode
use m_pblrfNode, only: pblrfNode_typecast
use m_pblrfNode, only: pblrfNode_nextcast
use m_obsdiagNode, only: obsdiagNode_set
implicit none

PRIVATE
PUBLIC intpblrf

contains

subroutine intpblrf(pblrfhead,rval,sval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intpblrf      apply nonlin qc obs operator for conv. pblrf
!   prgmmr: zhu           org: np23                date: 2011-02-20
!
! abstract: apply observation operator and adjoint for conventional pblrf
!           observations with nonlinear qc operator
!
! program history log:
!
!   2012-09-14  Syed RH Rizvi, NCAR/NESL/MMM/DAS  - introduced ladtest_obs         
!   2014-12-03  derber  - modify so that use of obsdiags can be turned off
!
!   input argument list:
!     pblrfhead
!     spblrf    - increment in grid space
!     rpblrf
!
!   output argument list:
!     rpblrf    - results from observation operator (0 for no data)
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
  class(obsNode  ),pointer, intent(in   ) :: pblrfhead
  type(gsi_bundle),         intent(in   ) :: sval
  type(gsi_bundle),         intent(inout) :: rval

! Declare local variables
  integer(i_kind) ier,istatus
  integer(i_kind) j1,j2,j3,j4
! real(r_kind) penalty
  real(r_kind) w1,w2,w3,w4
  real(r_kind) val
  real(r_kind) cg_pblrf,p0,grad,wnotgross,wgross,pg_pblrf
  real(r_kind),pointer,dimension(:) :: spblrf
  real(r_kind),pointer,dimension(:) :: rpblrf
  type(pblrfNode), pointer :: pblrfptr

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'pblrf',spblrf,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'pblrf',rpblrf,istatus);ier=istatus+ier
  if(ier/=0)return

  !pblrfptr => pblrfhead
  pblrfptr => pblrfNode_typecast(pblrfhead)
  do while (associated(pblrfptr))
     j1=pblrfptr%ij(1)
     j2=pblrfptr%ij(2)
     j3=pblrfptr%ij(3)
     j4=pblrfptr%ij(4)
     w1=pblrfptr%wij(1)
     w2=pblrfptr%wij(2)
     w3=pblrfptr%wij(3)
     w4=pblrfptr%wij(4)

!    Forward model
     val=w1*spblrf(j1)+w2*spblrf(j2)&
        +w3*spblrf(j3)+w4*spblrf(j4)

     if(luse_obsdiag)then
        if (lsaveobsens) then
           grad = val*pblrfptr%raterr2*pblrfptr%err2
           !-- pblrfptr%diags%obssen(jiter) = grad
           call obsdiagNode_set(pblrfptr%diags,jiter=jiter,obssen=grad)
        else
           !-- if (pblrfptr%luse) pblrfptr%diags%tldepart(jiter)=val
           if (pblrfptr%luse) call obsdiagNode_set(pblrfptr%diags,jiter=jiter,tldepart=val)
        endif
     endif

     if (l_do_adjoint) then
        if (.not. lsaveobsens) then
           if( .not. ladtest_obs)   val=val-pblrfptr%res

!          gradient of nonlinear operator
           if (nlnqc_iter .and. pblrfptr%pg > tiny_r_kind .and. &
                                pblrfptr%b  > tiny_r_kind) then
              pg_pblrf=pblrfptr%pg*varqc_iter
              cg_pblrf=cg_term/pblrfptr%b
              wnotgross= one-pg_pblrf
              wgross = pg_pblrf*cg_pblrf/wnotgross
              p0   = wgross/(wgross+exp(-half*pblrfptr%err2*val**2))
              val = val*(one-p0)
           endif
           if( ladtest_obs ) then
           grad = val
           else
           grad = val*pblrfptr%raterr2*pblrfptr%err2
           end if
        endif

!       Adjoint
        rpblrf(j1)=rpblrf(j1)+w1*grad
        rpblrf(j2)=rpblrf(j2)+w2*grad
        rpblrf(j3)=rpblrf(j3)+w3*grad
        rpblrf(j4)=rpblrf(j4)+w4*grad
     endif

     !pblrfptr => pblrfptr%llpoint
     pblrfptr => pblrfNode_nextcast(pblrfptr)

  end do

  return
end subroutine intpblrf

end module intpblrfmod
