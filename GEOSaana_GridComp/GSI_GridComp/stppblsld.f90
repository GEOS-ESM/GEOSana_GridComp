module stppblsldmod

!$$$ module documentation block
!           .      .    .                                       .
! module:   stppblsldmod    module for stppblsld
!  prgmmr:
!
! abstract: module for stppblsld
!
! program history log:
!   2009-02-24  zhu
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!
! subroutines included:
!   sub stppblsld
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

PRIVATE
PUBLIC stppblsld

contains

subroutine stppblsld(pblsldhead,rval,sval,out,sges,nstep)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    stppblsld      calculate penalty and contribution to stepsize
!   prgmmr: derber           org: np23                date: 2004-07-20
!
! abstract: calculate penalty and contribution to stepsize for surface pressure
!            with addition of nonlinear qc
!
! program history log:
!   2011-02-23  zhu  - update
!
!   input argument list:
!     pblsldhead
!     rpblsld     - search direction for pblsld
!     spblsld     - analysis increment for pblsld
!     sges     - step size estimate (nstep)
!     nstep    - number of stepsizes  (==0 means use outer iteration values)
!                                         
!   output argument list:         
!     out(1:nstep)   - contribution to penalty for conventional pblsld - sges(1:nstep)
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
  use m_pblsldNode, only: pblsldNode
  use m_pblsldNode, only: pblsldNode_typecast
  use m_pblsldNode, only: pblsldNode_nextcast
  implicit none

! Declare passed variables
  class(obsNode), pointer             ,intent(in   ) :: pblsldhead
  integer(i_kind)                     ,intent(in   ) :: nstep
  real(r_quad),dimension(max(1,nstep)),intent(inout) :: out
  type(gsi_bundle)                    ,intent(in   ) :: rval,sval
  real(r_kind),dimension(max(1,nstep)),intent(in   ) :: sges

! Declare local variables  
  integer(i_kind) j1,j2,j3,j4,kk,ier,istatus
  real(r_kind) w1,w2,w3,w4
  real(r_kind) val,val2
  real(r_kind) cg_pblsld,pblsld,wgross,wnotgross
  real(r_kind),dimension(max(1,nstep)):: pen
  real(r_kind) pg_pblsld
  real(r_kind),pointer,dimension(:) :: spblsld
  real(r_kind),pointer,dimension(:) :: rpblsld
  type(pblsldNode), pointer :: pblsldptr

  out=zero_quad

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'pblsld',spblsld,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'pblsld',rpblsld,istatus);ier=istatus+ier
  if(ier/=0)return

  pblsldptr => pblsldNode_typecast(pblsldhead)
  do while (associated(pblsldptr))
     if(pblsldptr%luse)then
        if(nstep > 0)then
           j1=pblsldptr%ij(1)
           j2=pblsldptr%ij(2)
           j3=pblsldptr%ij(3)
           j4=pblsldptr%ij(4)
           w1=pblsldptr%wij(1)
           w2=pblsldptr%wij(2)
           w3=pblsldptr%wij(3)
           w4=pblsldptr%wij(4)

           val =w1*rpblsld(j1)+w2*rpblsld(j2)+w3*rpblsld(j3)+w4*rpblsld(j4)
           val2=w1*spblsld(j1)+w2*spblsld(j2)+w3*spblsld(j3)+w4*spblsld(j4)-pblsldptr%res

           do kk=1,nstep
              pblsld=val2+sges(kk)*val
              pen(kk)= pblsld*pblsld*pblsldptr%err2
           end do
        else
           pen(1)=pblsldptr%res*pblsldptr%res*pblsldptr%err2
        end if
 
!  Modify penalty term if nonlinear QC
        if (nlnqc_iter .and. pblsldptr%pg > tiny_r_kind .and.  &
                             pblsldptr%b  > tiny_r_kind) then
           pg_pblsld=pblsldptr%pg*varqc_iter
           cg_pblsld=cg_term/pblsldptr%b
           wnotgross= one-pg_pblsld
           wgross = pg_pblsld*cg_pblsld/wnotgross
           do kk=1,max(1,nstep)
              pen(kk)= -two*log((exp(-half*pen(kk)) + wgross)/(one+wgross))
           end do
        endif

        out(1) = out(1)+pen(1)*pblsldptr%raterr2
        do kk=2,nstep
           out(kk) = out(kk)+(pen(kk)-pen(1))*pblsldptr%raterr2
        end do
     end if

     pblsldptr => pblsldNode_nextcast(pblsldptr)

  end do
  
  return
end subroutine stppblsld

end module stppblsldmod
