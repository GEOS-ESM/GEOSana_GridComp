module stppblgldmod

!$$$ module documentation block
!           .      .    .                                       .
! module:   stppblgldmod    module for stppblgld
!  prgmmr:
!
! abstract: module for stppblgld
!
! program history log:
!   2009-02-24  zhu
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!
! subroutines included:
!   sub stppblgld
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

PRIVATE
PUBLIC stppblgld

contains

subroutine stppblgld(pblgldhead,rval,sval,out,sges,nstep)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    stppblgld      calculate penalty and contribution to stepsize
!   prgmmr: derber           org: np23                date: 2004-07-20
!
! abstract: calculate penalty and contribution to stepsize for surface pressure
!            with addition of nonlinear qc
!
! program history log:
!   2011-02-23  zhu  - update
!
!   input argument list:
!     pblgldhead
!     rpblgld     - search direction for pblgld
!     spblgld     - analysis increment for pblgld
!     sges     - step size estimate (nstep)
!     nstep    - number of stepsizes  (==0 means use outer iteration values)
!                                         
!   output argument list:         
!     out(1:nstep)   - contribution to penalty for conventional pblgld - sges(1:nstep)
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
  use m_pblgldNode, only: pblgldNode
  use m_pblgldNode, only: pblgldNode_typecast
  use m_pblgldNode, only: pblgldNode_nextcast
  implicit none

! Declare passed variables
  class(obsNode), pointer             ,intent(in   ) :: pblgldhead
  integer(i_kind)                     ,intent(in   ) :: nstep
  real(r_quad),dimension(max(1,nstep)),intent(inout) :: out
  type(gsi_bundle)                    ,intent(in   ) :: rval,sval
  real(r_kind),dimension(max(1,nstep)),intent(in   ) :: sges

! Declare local variables  
  integer(i_kind) j1,j2,j3,j4,kk,ier,istatus
  real(r_kind) w1,w2,w3,w4
  real(r_kind) val,val2
  real(r_kind) cg_pblgld,pblgld,wgross,wnotgross
  real(r_kind),dimension(max(1,nstep)):: pen
  real(r_kind) pg_pblgld
  real(r_kind),pointer,dimension(:) :: spblgld
  real(r_kind),pointer,dimension(:) :: rpblgld
  type(pblgldNode), pointer :: pblgldptr

  out=zero_quad

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'pblgld',spblgld,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'pblgld',rpblgld,istatus);ier=istatus+ier
  if(ier/=0)return

  pblgldptr => pblgldNode_typecast(pblgldhead)
  do while (associated(pblgldptr))
     if(pblgldptr%luse)then
        if(nstep > 0)then
           j1=pblgldptr%ij(1)
           j2=pblgldptr%ij(2)
           j3=pblgldptr%ij(3)
           j4=pblgldptr%ij(4)
           w1=pblgldptr%wij(1)
           w2=pblgldptr%wij(2)
           w3=pblgldptr%wij(3)
           w4=pblgldptr%wij(4)

           val =w1*rpblgld(j1)+w2*rpblgld(j2)+w3*rpblgld(j3)+w4*rpblgld(j4)
           val2=w1*spblgld(j1)+w2*spblgld(j2)+w3*spblgld(j3)+w4*spblgld(j4)-pblgldptr%res

           do kk=1,nstep
              pblgld=val2+sges(kk)*val
              pen(kk)= pblgld*pblgld*pblgldptr%err2
           end do
        else
           pen(1)=pblgldptr%res*pblgldptr%res*pblgldptr%err2
        end if
 
!  Modify penalty term if nonlinear QC
        if (nlnqc_iter .and. pblgldptr%pg > tiny_r_kind .and.  &
                             pblgldptr%b  > tiny_r_kind) then
           pg_pblgld=pblgldptr%pg*varqc_iter
           cg_pblgld=cg_term/pblgldptr%b
           wnotgross= one-pg_pblgld
           wgross = pg_pblgld*cg_pblgld/wnotgross
           do kk=1,max(1,nstep)
              pen(kk)= -two*log((exp(-half*pen(kk)) + wgross)/(one+wgross))
           end do
        endif

        out(1) = out(1)+pen(1)*pblgldptr%raterr2
        do kk=2,nstep
           out(kk) = out(kk)+(pen(kk)-pen(1))*pblgldptr%raterr2
        end do
     end if

     pblgldptr => pblgldNode_nextcast(pblgldptr)

  end do
  
  return
end subroutine stppblgld

end module stppblgldmod
