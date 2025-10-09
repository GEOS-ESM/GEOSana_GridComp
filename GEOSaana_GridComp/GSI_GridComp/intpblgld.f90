module intpblgldmod
!$$$ module documentation block
!           .      .    .                                       .
! module:   intpblgldmod    module for intpblgld and its tangent linear intpblgld_tl
!   prgmmr:
!
! abstract: module for intpblgld and its tangent linear intpblgld_tl
!
! program history log:
!
!   2012-09-14  Syed RH Rizvi, NCAR/NESL/MMM/DAS  - implemented obs adjoint test  
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!
! subroutines included:
!   sub intpblgld
!
! variable definitions:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use m_obsNode, only: obsNode
use m_pblgldNode, only: pblgldNode
use m_pblgldNode, only: pblgldNode_typecast
use m_pblgldNode, only: pblgldNode_nextcast
use m_obsdiagNode, only: obsdiagNode_set
implicit none

PRIVATE
PUBLIC intpblgld

contains

subroutine intpblgld(pblgldhead,rval,sval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intpblgld      apply nonlin qc obs operator for conv. pblgld
!   prgmmr: zhu           org: np23                date: 2011-02-20
!
! abstract: apply observation operator and adjoint for conventional pblgld
!           observations with nonlinear qc operator
!
! program history log:
!
!   2012-09-14  Syed RH Rizvi, NCAR/NESL/MMM/DAS  - introduced ladtest_obs         
!   2014-12-03  derber  - modify so that use of obsdiags can be turned off
!
!   input argument list:
!     pblgldhead
!     spblgld    - increment in grid space
!     rpblgld
!
!   output argument list:
!     rpblgld    - results from observation operator (0 for no data)
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
  class(obsNode  ),pointer, intent(in   ) :: pblgldhead
  type(gsi_bundle),         intent(in   ) :: sval
  type(gsi_bundle),         intent(inout) :: rval

! Declare local variables
  integer(i_kind) ier,istatus
  integer(i_kind) j1,j2,j3,j4
! real(r_kind) penalty
  real(r_kind) w1,w2,w3,w4
  real(r_kind) val
  real(r_kind) cg_pblgld,p0,grad,wnotgross,wgross,pg_pblgld
  real(r_kind),pointer,dimension(:) :: spblgld
  real(r_kind),pointer,dimension(:) :: rpblgld
  type(pblgldNode), pointer :: pblgldptr

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'pblgld',spblgld,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'pblgld',rpblgld,istatus);ier=istatus+ier
  if(ier/=0)return

  !pblgldptr => pblgldhead
  pblgldptr => pblgldNode_typecast(pblgldhead)
  do while (associated(pblgldptr))
     j1=pblgldptr%ij(1)
     j2=pblgldptr%ij(2)
     j3=pblgldptr%ij(3)
     j4=pblgldptr%ij(4)
     w1=pblgldptr%wij(1)
     w2=pblgldptr%wij(2)
     w3=pblgldptr%wij(3)
     w4=pblgldptr%wij(4)

!    Forward model
     val=w1*spblgld(j1)+w2*spblgld(j2)&
        +w3*spblgld(j3)+w4*spblgld(j4)

     if(luse_obsdiag)then
        if (lsaveobsens) then
           grad = val*pblgldptr%raterr2*pblgldptr%err2
           !-- pblgldptr%diags%obssen(jiter) = grad
           call obsdiagNode_set(pblgldptr%diags,jiter=jiter,obssen=grad)
        else
           !-- if (pblgldptr%luse) pblgldptr%diags%tldepart(jiter)=val
           if (pblgldptr%luse) call obsdiagNode_set(pblgldptr%diags,jiter=jiter,tldepart=val)
        endif
     endif

     if (l_do_adjoint) then
        if (.not. lsaveobsens) then
           if( .not. ladtest_obs)   val=val-pblgldptr%res

!          gradient of nonlinear operator
           if (nlnqc_iter .and. pblgldptr%pg > tiny_r_kind .and. &
                                pblgldptr%b  > tiny_r_kind) then
              pg_pblgld=pblgldptr%pg*varqc_iter
              cg_pblgld=cg_term/pblgldptr%b
              wnotgross= one-pg_pblgld
              wgross = pg_pblgld*cg_pblgld/wnotgross
              p0   = wgross/(wgross+exp(-half*pblgldptr%err2*val**2))
              val = val*(one-p0)
           endif
           if( ladtest_obs ) then
              grad = val
           else
              grad = val*pblgldptr%raterr2*pblgldptr%err2
           end if
        endif

!       Adjoint
        rpblgld(j1)=rpblgld(j1)+w1*grad
        rpblgld(j2)=rpblgld(j2)+w2*grad
        rpblgld(j3)=rpblgld(j3)+w3*grad
        rpblgld(j4)=rpblgld(j4)+w4*grad
     endif

     !pblgldptr => pblgldptr%llpoint
     pblgldptr => pblgldNode_nextcast(pblgldptr)

  end do

  return
end subroutine intpblgld

end module intpblgldmod
