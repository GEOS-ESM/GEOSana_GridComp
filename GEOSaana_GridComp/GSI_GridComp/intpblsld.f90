module intpblsldmod
!$$$ module documentation block
!           .      .    .                                       .
! module:   intpblsldmod    module for intpblsld and its tangent linear intpblsld_tl
!   prgmmr:
!
! abstract: module for intpblsld and its tangent linear intpblsld_tl
!
! program history log:
!
!   2012-09-14  Syed RH Rizvi, NCAR/NESL/MMM/DAS  - implemented obs adjoint test  
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!
! subroutines included:
!   sub intpblsld
!
! variable definitions:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use m_obsNode, only: obsNode
use m_pblsldNode, only: pblsldNode
use m_pblsldNode, only: pblsldNode_typecast
use m_pblsldNode, only: pblsldNode_nextcast
use m_obsdiagNode, only: obsdiagNode_set
implicit none

PRIVATE
PUBLIC intpblsld

contains

subroutine intpblsld(pblsldhead,rval,sval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intpblsld      apply nonlin qc obs operator for conv. pblsld
!   prgmmr: zhu           org: np23                date: 2011-02-20
!
! abstract: apply observation operator and adjoint for conventional pblsld
!           observations with nonlinear qc operator
!
! program history log:
!
!   2012-09-14  Syed RH Rizvi, NCAR/NESL/MMM/DAS  - introduced ladtest_obs         
!   2014-12-03  derber  - modify so that use of obsdiags can be turned off
!
!   input argument list:
!     pblsldhead
!     spblsld    - increment in grid space
!     rpblsld
!
!   output argument list:
!     rpblsld    - results from observation operator (0 for no data)
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
  class(obsNode  ),pointer, intent(in   ) :: pblsldhead
  type(gsi_bundle),         intent(in   ) :: sval
  type(gsi_bundle),         intent(inout) :: rval

! Declare local variables
  integer(i_kind) ier,istatus
  integer(i_kind) j1,j2,j3,j4
! real(r_kind) penalty
  real(r_kind) w1,w2,w3,w4
  real(r_kind) val
  real(r_kind) cg_pblsld,p0,grad,wnotgross,wgross,pg_pblsld
  real(r_kind),pointer,dimension(:) :: spblsld
  real(r_kind),pointer,dimension(:) :: rpblsld
  type(pblsldNode), pointer :: pblsldptr

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'pblsld',spblsld,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'pblsld',rpblsld,istatus);ier=istatus+ier
  if(ier/=0)return

  !pblsldptr => pblsldhead
  pblsldptr => pblsldNode_typecast(pblsldhead)
  do while (associated(pblsldptr))
     j1=pblsldptr%ij(1)
     j2=pblsldptr%ij(2)
     j3=pblsldptr%ij(3)
     j4=pblsldptr%ij(4)
     w1=pblsldptr%wij(1)
     w2=pblsldptr%wij(2)
     w3=pblsldptr%wij(3)
     w4=pblsldptr%wij(4)

!    Forward model
     val=w1*spblsld(j1)+w2*spblsld(j2)&
        +w3*spblsld(j3)+w4*spblsld(j4)

     if(luse_obsdiag)then
        if (lsaveobsens) then
           grad = val*pblsldptr%raterr2*pblsldptr%err2
           !-- pblsldptr%diags%obssen(jiter) = grad
           call obsdiagNode_set(pblsldptr%diags,jiter=jiter,obssen=grad)
        else
           !-- if (pblsldptr%luse) pblsldptr%diags%tldepart(jiter)=val
           if (pblsldptr%luse) call obsdiagNode_set(pblsldptr%diags,jiter=jiter,tldepart=val)
        endif
     endif

     if (l_do_adjoint) then
        if (.not. lsaveobsens) then
           if( .not. ladtest_obs)   val=val-pblsldptr%res

!          gradient of nonlinear operator
           if (nlnqc_iter .and. pblsldptr%pg > tiny_r_kind .and. &
                                pblsldptr%b  > tiny_r_kind) then
              pg_pblsld=pblsldptr%pg*varqc_iter
              cg_pblsld=cg_term/pblsldptr%b
              wnotgross= one-pg_pblsld
              wgross = pg_pblsld*cg_pblsld/wnotgross
              p0   = wgross/(wgross+exp(-half*pblsldptr%err2*val**2))
              val = val*(one-p0)
           endif
           if( ladtest_obs ) then
              grad = val
           else
              grad = val*pblsldptr%raterr2*pblsldptr%err2
           end if
        endif

!       Adjoint
        rpblsld(j1)=rpblsld(j1)+w1*grad
        rpblsld(j2)=rpblsld(j2)+w2*grad
        rpblsld(j3)=rpblsld(j3)+w3*grad
        rpblsld(j4)=rpblsld(j4)+w4*grad
     endif

     !pblsldptr => pblsldptr%llpoint
     pblsldptr => pblsldNode_nextcast(pblsldptr)

  end do

  return
end subroutine intpblsld

end module intpblsldmod
