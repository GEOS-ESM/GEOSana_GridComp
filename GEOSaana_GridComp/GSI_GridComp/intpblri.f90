module intpblrimod
!$$$ module documentation block
!           .      .    .                                       .
! module:   intpblrimod    module for intpblri and its tangent linear intpblri_tl
!   prgmmr:
!
! abstract: module for intpblri and its tangent linear intpblri_tl
!
! program history log:
!
!   2012-09-14  Syed RH Rizvi, NCAR/NESL/MMM/DAS  - implemented obs adjoint test  
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!
! subroutines included:
!   sub intpblri
!
! variable definitions:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

use m_obsNode, only: obsNode
use m_pblriNode, only: pblriNode
use m_pblriNode, only: pblriNode_typecast
use m_pblriNode, only: pblriNode_nextcast
use m_obsdiagNode, only: obsdiagNode_set
implicit none

PRIVATE
PUBLIC intpblri

contains

subroutine intpblri(pblrihead,rval,sval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intpblri      apply nonlin qc obs operator for conv. pblri
!   prgmmr: zhu           org: np23                date: 2011-02-20
!
! abstract: apply observation operator and adjoint for conventional pblri
!           observations with nonlinear qc operator
!
! program history log:
!
!   2012-09-14  Syed RH Rizvi, NCAR/NESL/MMM/DAS  - introduced ladtest_obs         
!   2014-12-03  derber  - modify so that use of obsdiags can be turned off
!
!   input argument list:
!     pblrihead
!     spblri    - increment in grid space
!     rpblri
!
!   output argument list:
!     rpblri    - results from observation operator (0 for no data)
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
  class(obsNode  ),pointer, intent(in   ) :: pblrihead
  type(gsi_bundle),         intent(in   ) :: sval
  type(gsi_bundle),         intent(inout) :: rval

! Declare local variables
  integer(i_kind) ier,istatus
  integer(i_kind) j1,j2,j3,j4
! real(r_kind) penalty
  real(r_kind) w1,w2,w3,w4
  real(r_kind) val
  real(r_kind) cg_pblri,p0,grad,wnotgross,wgross,pg_pblri
  real(r_kind),pointer,dimension(:) :: spblri
  real(r_kind),pointer,dimension(:) :: rpblri
  type(pblriNode), pointer :: pblriptr

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'pblri',spblri,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'pblri',rpblri,istatus);ier=istatus+ier
  if(ier/=0)return

  !pblriptr => pblrihead
  pblriptr => pblriNode_typecast(pblrihead)
  do while (associated(pblriptr))
     j1=pblriptr%ij(1)
     j2=pblriptr%ij(2)
     j3=pblriptr%ij(3)
     j4=pblriptr%ij(4)
     w1=pblriptr%wij(1)
     w2=pblriptr%wij(2)
     w3=pblriptr%wij(3)
     w4=pblriptr%wij(4)

!    Forward model
     val=w1*spblri(j1)+w2*spblri(j2)&
        +w3*spblri(j3)+w4*spblri(j4)

     if(luse_obsdiag)then
        if (lsaveobsens) then
           grad = val*pblriptr%raterr2*pblriptr%err2
           !-- pblriptr%diags%obssen(jiter) = grad
           call obsdiagNode_set(pblriptr%diags,jiter=jiter,obssen=grad)
        else
           !-- if (pblriptr%luse) pblriptr%diags%tldepart(jiter)=val
           if (pblriptr%luse) call obsdiagNode_set(pblriptr%diags,jiter=jiter,tldepart=val)
        endif
     endif

     if (l_do_adjoint) then
        if (.not. lsaveobsens) then
           if( .not. ladtest_obs)   val=val-pblriptr%res

!          gradient of nonlinear operator
           if (nlnqc_iter .and. pblriptr%pg > tiny_r_kind .and. &
                                pblriptr%b  > tiny_r_kind) then
              pg_pblri=pblriptr%pg*varqc_iter
              cg_pblri=cg_term/pblriptr%b
              wnotgross= one-pg_pblri
              wgross = pg_pblri*cg_pblri/wnotgross
              p0   = wgross/(wgross+exp(-half*pblriptr%err2*val**2))
              val = val*(one-p0)
           endif
           if( ladtest_obs ) then
              grad = val
           else
              grad = val*pblriptr%raterr2*pblriptr%err2
           end if
        endif

!       Adjoint
        rpblri(j1)=rpblri(j1)+w1*grad
        rpblri(j2)=rpblri(j2)+w2*grad
        rpblri(j3)=rpblri(j3)+w3*grad
        rpblri(j4)=rpblri(j4)+w4*grad
     endif

     !pblriptr => pblriptr%llpoint
     pblriptr => pblriNode_nextcast(pblriptr)

  end do

  return
end subroutine intpblri

end module intpblrimod
