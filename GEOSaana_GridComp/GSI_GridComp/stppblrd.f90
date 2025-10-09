module stppblrdmod

!$$$ module documentation block
!           .      .    .                                       .
! module:   stppblrdmod    module for stppblrd
!  prgmmr:
!
! abstract: module for stppblrd
!
! program history log:
!   2009-02-24  zhu
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!
! subroutines included:
!   sub stppblrd
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

PRIVATE
PUBLIC stppblrd

contains

subroutine stppblrd(pblrdhead,rval,sval,out,sges,nstep)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    stppblrd      calculate penalty and contribution to stepsize
!   prgmmr: derber           org: np23                date: 2004-07-20
!
! abstract: calculate penalty and contribution to stepsize for surface pressure
!            with addition of nonlinear qc
!
! program history log:
!   2011-02-23  zhu  - update
!
!   input argument list:
!     pblrdhead
!     rpblrd     - search direction for pblrd
!     spblrd     - analysis increment for pblrd
!     sges     - step size estimate (nstep)
!     nstep    - number of stepsizes  (==0 means use outer iteration values)
!                                         
!   output argument list:         
!     out(1:nstep)   - contribution to penalty for conventional pblrd - sges(1:nstep)
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind,r_quad
  use qcmod, only: nlnqc_iter,varqc_iter
  use constants, only: half,one,two,tiny_r_kind,cg_term,zero_quad
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  use m_obsNode, only: obsNode
  use m_pblrdNode, only: pblrdNode
  use m_pblrdNode, only: pblrdNode_typecast
  use m_pblrdNode, only: pblrdNode_nextcast
  implicit none

! Declare passed variables
  class(obsNode), pointer             ,intent(in   ) :: pblrdhead
  integer(i_kind)                     ,intent(in   ) :: nstep
  real(r_quad),dimension(max(1,nstep)),intent(inout) :: out
  type(gsi_bundle)                    ,intent(in   ) :: rval,sval
  real(r_kind),dimension(max(1,nstep)),intent(in   ) :: sges

! Declare local variables  
  integer(i_kind) j1,j2,j3,j4,kk,ier,istatus
  real(r_kind) w1,w2,w3,w4
  real(r_kind) val,val2
  real(r_kind) cg_pblrd,pblrd,wgross,wnotgross
  real(r_kind),dimension(max(1,nstep)):: pen
  real(r_kind) pg_pblrd
  real(r_kind),pointer,dimension(:) :: spblrd
  real(r_kind),pointer,dimension(:) :: rpblrd
  type(pblrdNode), pointer :: pblrdptr

  out=zero_quad

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'pblrd',spblrd,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'pblrd',rpblrd,istatus);ier=istatus+ier
  if(ier/=0)return

  pblrdptr => pblrdNode_typecast(pblrdhead)
  do while (associated(pblrdptr))
     if(pblrdptr%luse)then
        if(nstep > 0)then
           j1=pblrdptr%ij(1)
           j2=pblrdptr%ij(2)
           j3=pblrdptr%ij(3)
           j4=pblrdptr%ij(4)
           w1=pblrdptr%wij(1)
           w2=pblrdptr%wij(2)
           w3=pblrdptr%wij(3)
           w4=pblrdptr%wij(4)

           val =w1*rpblrd(j1)+w2*rpblrd(j2)+w3*rpblrd(j3)+w4*rpblrd(j4)
           val2=w1*spblrd(j1)+w2*spblrd(j2)+w3*spblrd(j3)+w4*spblrd(j4)-pblrdptr%res

           do kk=1,nstep
              pblrd=val2+sges(kk)*val
              pen(kk)= pblrd*pblrd*pblrdptr%err2
           end do
        else
           pen(1)=pblrdptr%res*pblrdptr%res*pblrdptr%err2
        end if
 
!  Modify penalty term if nonlinear QC
        if (nlnqc_iter .and. pblrdptr%pg > tiny_r_kind .and.  &
                             pblrdptr%b  > tiny_r_kind) then
           pg_pblrd=pblrdptr%pg*varqc_iter
           cg_pblrd=cg_term/pblrdptr%b
           wnotgross= one-pg_pblrd
           wgross = pg_pblrd*cg_pblrd/wnotgross
           do kk=1,max(1,nstep)
              pen(kk)= -two*log((exp(-half*pen(kk)) + wgross)/(one+wgross))
           end do
        endif

        out(1) = out(1)+pen(1)*pblrdptr%raterr2
        do kk=2,nstep
           out(kk) = out(kk)+(pen(kk)-pen(1))*pblrdptr%raterr2
        end do
     end if

     pblrdptr => pblrdNode_nextcast(pblrdptr)

  end do
  
  return
end subroutine stppblrd

end module stppblrdmod
