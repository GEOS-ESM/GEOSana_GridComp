module stppblkhmod

!$$$ module documentation block
!           .      .    .                                       .
! module:   stppblkhmod    module for stppblkh
!  prgmmr:
!
! abstract: module for stppblkh
!
! program history log:
!   2009-02-24  zhu
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!
! subroutines included:
!   sub stppblkh
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

PRIVATE
PUBLIC stppblkh

contains

subroutine stppblkh(pblkhhead,rval,sval,out,sges,nstep)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    stppblkh      calculate penalty and contribution to stepsize
!   prgmmr: derber           org: np23                date: 2004-07-20
!
! abstract: calculate penalty and contribution to stepsize for surface pressure
!            with addition of nonlinear qc
!
! program history log:
!   2011-02-23  zhu  - update
!
!   input argument list:
!     pblkhhead
!     rpblkh     - search direction for pblkh
!     spblkh     - analysis increment for pblkh
!     sges     - step size estimate (nstep)
!     nstep    - number of stepsizes  (==0 means use outer iteration values)
!                                         
!   output argument list:         
!     out(1:nstep)   - contribution to penalty for conventional pblkh - sges(1:nstep)
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind,r_quad
  use qcmod, only: nlnqc_iter,varqc_iter,vqc
  use constants, only: half,one,two,tiny_r_kind,cg_term,zero_quad
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  use m_obsNode, only: obsNode
  use m_pblkhNode, only: pblkhNode
  use m_pblkhNode, only: pblkhNode_typecast
  use m_pblkhNode, only: pblkhNode_nextcast
  implicit none

! Declare passed variables
  class(obsNode), pointer             ,intent(in   ) :: pblkhhead
  integer(i_kind)                     ,intent(in   ) :: nstep
  real(r_quad),dimension(max(1,nstep)),intent(inout) :: out
  type(gsi_bundle)                    ,intent(in   ) :: rval,sval
  real(r_kind),dimension(max(1,nstep)),intent(in   ) :: sges

! Declare local variables  
  integer(i_kind) j1,j2,j3,j4,kk,ier,istatus
  real(r_kind) w1,w2,w3,w4
  real(r_kind) val,val2
  real(r_kind) cg_pblkh,pblkh,wgross,wnotgross
  real(r_kind),dimension(max(1,nstep)):: pen
  real(r_kind) pg_pblkh
  real(r_kind),pointer,dimension(:) :: spblkh
  real(r_kind),pointer,dimension(:) :: rpblkh
  type(pblkhNode), pointer :: pblkhptr

  out=zero_quad

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'pblkh',spblkh,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'pblkh',rpblkh,istatus);ier=istatus+ier
  if(ier/=0)return

  pblkhptr => pblkhNode_typecast(pblkhhead)
  do while (associated(pblkhptr))
     if(pblkhptr%luse)then
        if(nstep > 0)then
           j1=pblkhptr%ij(1)
           j2=pblkhptr%ij(2)
           j3=pblkhptr%ij(3)
           j4=pblkhptr%ij(4)
           w1=pblkhptr%wij(1)
           w2=pblkhptr%wij(2)
           w3=pblkhptr%wij(3)
           w4=pblkhptr%wij(4)

           val =w1*rpblkh(j1)+w2*rpblkh(j2)+w3*rpblkh(j3)+w4*rpblkh(j4)
           val2=w1*spblkh(j1)+w2*spblkh(j2)+w3*spblkh(j3)+w4*spblkh(j4)-pblkhptr%res

           do kk=1,nstep
              pblkh=val2+sges(kk)*val
              pen(kk)= pblkh*pblkh*pblkhptr%err2
           end do
        else
           pen(1)=pblkhptr%res*pblkhptr%res*pblkhptr%err2
        end if
 
!  Modify penalty term if nonlinear QC
        if (vqc .and. nlnqc_iter .and. pblkhptr%pg > tiny_r_kind .and.  &
                             pblkhptr%b  > tiny_r_kind) then
           pg_pblkh=pblkhptr%pg*varqc_iter
           cg_pblkh=cg_term/pblkhptr%b
           wnotgross= one-pg_pblkh
           wgross = pg_pblkh*cg_pblkh/wnotgross
           do kk=1,max(1,nstep)
              pen(kk)= -two*log((exp(-half*pen(kk)) + wgross)/(one+wgross))
           end do
        endif

        out(1) = out(1)+pen(1)*pblkhptr%raterr2
        do kk=2,nstep
           out(kk) = out(kk)+(pen(kk)-pen(1))*pblkhptr%raterr2
        end do
     end if

     pblkhptr => pblkhNode_nextcast(pblkhptr)

  end do
  
  return
end subroutine stppblkh

end module stppblkhmod
