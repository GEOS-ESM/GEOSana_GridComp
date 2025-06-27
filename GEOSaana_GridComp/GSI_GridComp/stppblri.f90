module stppblrimod

!$$$ module documentation block
!           .      .    .                                       .
! module:   stppblrimod    module for stppblri
!  prgmmr:
!
! abstract: module for stppblri
!
! program history log:
!   2009-02-24  zhu
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!
! subroutines included:
!   sub stppblri
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

PRIVATE
PUBLIC stppblri

contains

subroutine stppblri(pblrihead,rval,sval,out,sges,nstep)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    stppblri      calculate penalty and contribution to stepsize
!   prgmmr: derber           org: np23                date: 2004-07-20
!
! abstract: calculate penalty and contribution to stepsize for surface pressure
!            with addition of nonlinear qc
!
! program history log:
!   2011-02-23  zhu  - update
!
!   input argument list:
!     pblrihead
!     rpblri     - search direction for pblri
!     spblri     - analysis increment for pblri
!     sges     - step size estimate (nstep)
!     nstep    - number of stepsizes  (==0 means use outer iteration values)
!                                         
!   output argument list:         
!     out(1:nstep)   - contribution to penalty for conventional pblri - sges(1:nstep)
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
  use m_pblriNode, only: pblriNode
  use m_pblriNode, only: pblriNode_typecast
  use m_pblriNode, only: pblriNode_nextcast
  implicit none

! Declare passed variables
  class(obsNode), pointer             ,intent(in   ) :: pblrihead
  integer(i_kind)                     ,intent(in   ) :: nstep
  real(r_quad),dimension(max(1,nstep)),intent(inout) :: out
  type(gsi_bundle)                    ,intent(in   ) :: rval,sval
  real(r_kind),dimension(max(1,nstep)),intent(in   ) :: sges

! Declare local variables  
  integer(i_kind) j1,j2,j3,j4,kk,ier,istatus
  real(r_kind) w1,w2,w3,w4
  real(r_kind) val,val2
  real(r_kind) cg_pblri,pblri,wgross,wnotgross
  real(r_kind),dimension(max(1,nstep)):: pen
  real(r_kind) pg_pblri
  real(r_kind),pointer,dimension(:) :: spblri
  real(r_kind),pointer,dimension(:) :: rpblri
  type(pblriNode), pointer :: pblriptr

  out=zero_quad

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'pblri',spblri,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'pblri',rpblri,istatus);ier=istatus+ier
  if(ier/=0)return

  pblriptr => pblriNode_typecast(pblrihead)
  do while (associated(pblriptr))
     if(pblriptr%luse)then
        if(nstep > 0)then
           j1=pblriptr%ij(1)
           j2=pblriptr%ij(2)
           j3=pblriptr%ij(3)
           j4=pblriptr%ij(4)
           w1=pblriptr%wij(1)
           w2=pblriptr%wij(2)
           w3=pblriptr%wij(3)
           w4=pblriptr%wij(4)

           val =w1*rpblri(j1)+w2*rpblri(j2)+w3*rpblri(j3)+w4*rpblri(j4)
           val2=w1*spblri(j1)+w2*spblri(j2)+w3*spblri(j3)+w4*spblri(j4)-pblriptr%res

           do kk=1,nstep
              pblri=val2+sges(kk)*val
              pen(kk)= pblri*pblri*pblriptr%err2
           end do
        else
           pen(1)=pblriptr%res*pblriptr%res*pblriptr%err2
        end if
 
!  Modify penalty term if nonlinear QC
        if (vqc .and. nlnqc_iter .and. pblriptr%pg > tiny_r_kind .and.  &
                             pblriptr%b  > tiny_r_kind) then
           pg_pblri=pblriptr%pg*varqc_iter
           cg_pblri=cg_term/pblriptr%b
           wnotgross= one-pg_pblri
           wgross = pg_pblri*cg_pblri/wnotgross
           do kk=1,max(1,nstep)
              pen(kk)= -two*log((exp(-half*pen(kk)) + wgross)/(one+wgross))
           end do
        endif

        out(1) = out(1)+pen(1)*pblriptr%raterr2
        do kk=2,nstep
           out(kk) = out(kk)+(pen(kk)-pen(1))*pblriptr%raterr2
        end do
     end if

     pblriptr => pblriNode_nextcast(pblriptr)

  end do
  
  return
end subroutine stppblri

end module stppblrimod
