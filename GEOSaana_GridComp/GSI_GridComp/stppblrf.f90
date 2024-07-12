module stppblrfmod

!$$$ module documentation block
!           .      .    .                                       .
! module:   stppblrfmod    module for stppblrf
!  prgmmr:
!
! abstract: module for stppblrf
!
! program history log:
!   2009-02-24  zhu
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!
! subroutines included:
!   sub stppblrf
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

PRIVATE
PUBLIC stppblrf

contains

subroutine stppblrf(pblrfhead,rval,sval,out,sges,nstep)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    stppblrf      calculate penalty and contribution to stepsize
!   prgmmr: derber           org: np23                date: 2004-07-20
!
! abstract: calculate penalty and contribution to stepsize for surface pressure
!            with addition of nonlinear qc
!
! program history log:
!   2011-02-23  zhu  - update
!
!   input argument list:
!     pblrfhead
!     rpblrf     - search direction for pblrf
!     spblrf     - analysis increment for pblrf
!     sges     - step size estimate (nstep)
!     nstep    - number of stepsizes  (==0 means use outer iteration values)
!                                         
!   output argument list:         
!     out(1:nstep)   - contribution to penalty for conventional pblrf - sges(1:nstep)
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
  use m_pblrfNode, only: pblrfNode
  use m_pblrfNode, only: pblrfNode_typecast
  use m_pblrfNode, only: pblrfNode_nextcast
  implicit none

! Declare passed variables
  class(obsNode), pointer             ,intent(in   ) :: pblrfhead
  integer(i_kind)                     ,intent(in   ) :: nstep
  real(r_quad),dimension(max(1,nstep)),intent(inout) :: out
  type(gsi_bundle)                    ,intent(in   ) :: rval,sval
  real(r_kind),dimension(max(1,nstep)),intent(in   ) :: sges

! Declare local variables  
  integer(i_kind) j1,j2,j3,j4,kk,ier,istatus
  real(r_kind) w1,w2,w3,w4
  real(r_kind) val,val2
  real(r_kind) cg_pblrf,pblrf,wgross,wnotgross
  real(r_kind),dimension(max(1,nstep)):: pen
  real(r_kind) pg_pblrf
  real(r_kind),pointer,dimension(:) :: spblrf
  real(r_kind),pointer,dimension(:) :: rpblrf
  type(pblrfNode), pointer :: pblrfptr

  out=zero_quad

! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'pblrf',spblrf,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'pblrf',rpblrf,istatus);ier=istatus+ier
  if(ier/=0)return

  pblrfptr => pblrfNode_typecast(pblrfhead)
  do while (associated(pblrfptr))
     if(pblrfptr%luse)then
        if(nstep > 0)then
           j1=pblrfptr%ij(1)
           j2=pblrfptr%ij(2)
           j3=pblrfptr%ij(3)
           j4=pblrfptr%ij(4)
           w1=pblrfptr%wij(1)
           w2=pblrfptr%wij(2)
           w3=pblrfptr%wij(3)
           w4=pblrfptr%wij(4)

           val =w1*rpblrf(j1)+w2*rpblrf(j2)+w3*rpblrf(j3)+w4*rpblrf(j4)
           val2=w1*spblrf(j1)+w2*spblrf(j2)+w3*spblrf(j3)+w4*spblrf(j4)-pblrfptr%res

           do kk=1,nstep
              pblrf=val2+sges(kk)*val
              pen(kk)= pblrf*pblrf*pblrfptr%err2
              print*, "kk=",kk,",pen(kk)=",pen(kk), ",pblrf=",pblrf,",sges(kk)=",sges(kk)
           end do
        else
           pen(1)=pblrfptr%res*pblrfptr%res*pblrfptr%err2
           print*, "pen(1)=",pen
        end if
 
!  Modify penalty term if nonlinear QC
        if (vqc .and. nlnqc_iter .and. pblrfptr%pg > tiny_r_kind .and.  &
                             pblrfptr%b  > tiny_r_kind) then
           pg_pblrf=pblrfptr%pg*varqc_iter
           cg_pblrf=cg_term/pblrfptr%b
           wnotgross= one-pg_pblrf
           wgross = pg_pblrf*cg_pblrf/wnotgross
           do kk=1,max(1,nstep)
              pen(kk)= -two*log((exp(-half*pen(kk)) + wgross)/(one+wgross))
              print*, "L150: kk=",kk,",pen(kk)=",pen(kk)
           end do
        endif

        out(1) = out(1)+pen(1)*pblrfptr%raterr2
        do kk=2,nstep
           out(kk) = out(kk)+(pen(kk)-pen(1))*pblrfptr%raterr2
           print*, "L150: kk=",kk,",out(kk)=",out(kk)
        end do
     end if

     pblrfptr => pblrfNode_nextcast(pblrfptr)

  end do
  
  return
end subroutine stppblrf

end module stppblrfmod
