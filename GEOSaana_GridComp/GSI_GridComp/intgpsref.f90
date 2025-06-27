module intgpsrefmod

!$$$ module documentation block
!           .      .    .                                       .
! module:   intrefmod    module for intref and its tangent linear intref_tl
!   prgmmr:
!
! abstract: module for intref and its tangent linear intref_tl
!
! program history log:
!   2005-05-13  Yanqiu zhu - wrap intref and its tangent linear intref_tl into one module
!   2005-11-16  Derber - remove interfaces
!   2008-11-28  Todling - add interface back
!   2009-08-13  lueken - update documentation
!   2013-10-28  todling - rename p3d to prse
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!
! subroutines included:
!   sub intgpsref_
!
! variable definitions:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
use m_obsNode, only: obsNode
use m_gpsrefNode, only: gpsrefNode
use m_gpsrefNode, only: gpsrefNode_typecast
use m_gpsrefNode, only: gpsrefNode_nextcast
use m_obsdiagNode, only: obsdiagNode_set
implicit none

PRIVATE
PUBLIC intgpsref

interface intgpsref; module procedure &
          intgpsref_
end interface

contains

subroutine intgpsref_(gpshead,rval,sval)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    intref      apply nonlinqc obs operator refractivity
!   prgmmr: cucurull, l.     org: JCSDA/NCEP          date: 2004-04-29
!
! abstract: apply gps local refractivity operator and adjoint with
!           addition of nonlinear qc.
!
! program history log:
!   2004-04-29  cucurull- original code
!   2004-06-21  treadon - update documentation
!   2004-08-02  treadon - add only to module use, add intent in/out
!   2004-10-08  parrish - add nonlinear qc option
!   2004-11-19  cucurull- add increments for surface pressure and temperature at levels
!                          below observation. Install non-linear forward operator.
!   2005-01-26  cucurull- Implement local GPS RO operator
!   2005-03-01  parrish - nonlinear qc change as above; correct bug in zeroing of tl_AD
!   2005-03-23  cucurull- correct bounds for obs below the second level; place 
!                         bounds for k1 and k2
!   2005-04-11  treadon - merge intref and intref_qc into single routine
!   2005-08-02  derber  - modify for variational qc parameters for each ob
!   2005-09-28  derber  - consolidate location and weight arrays
!   2005-12-02  cucurull - fix bug for dimensions of sp and rp
!   2006-01-03  treadon - include r_kind type in w1,w2,...,w12 declaration
!   2006-07-28  derber  - modify to use new inner loop obs data structure
!                       - unify NL qc
!   2006-09-06  cucurull - generalize code to hybrid vertical coordinate and modify to use
!                          surface pressure
!   2007-01-13  derber - clean up code and use coding standards
!   2007-03-28  derber - turn intref into generalized intgps
!   2007-07-26  cucurull - in/out 3d pressure to update code to generalized vertical coordinate
!   2008-06-02  safford - rm unused vars
!   2008-11-26  todling  - add 4dvar and GSI adjoint capability (obs binning, obsens, etc) 
!   2008-11-26  todling - turned FOTO optional; changed ptr%time handle
!   2010-05-13  todling - update to use gsi_bundle; update interface
!   2014-12-03  derber  - modify so that use of obsdiags can be turned off
!   2025-02-07  eyang   - add jac, wij, ij, res for lower and upper obs to assimilate refractivity gradient
!   
!   input argument list:
!     gpshead  - obs type pointer to obs structure
!     st       - input temperature correction field
!     sq       - input q correction field
!     sp       - input (3D) p correction field
!
!   output argument list:
!     gpshead  - obs type pointer to obs structure
!     rt       - output t vector after inclusion of gps local refractivity
!     rq       - output q vector after inclusion of gps local refractivity
!     rp       - output p vector after inclusion of gps local refractivity
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind
  use obsmod, only: lsaveobsens,l_do_adjoint,luse_obsdiag
  use qcmod, only: nlnqc_iter,varqc_iter
  use gridmod, only: nsig
  use constants, only: zero,one,half,tiny_r_kind,cg_term,r3600
  use jfunc, only: jiter
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  use gsi_4dvar, only: ladtest_obs
  implicit none

! Declare passed variables
  class(obsNode),  pointer, intent(in   ) :: gpshead
  type(gsi_bundle),         intent(in   ) :: sval
  type(gsi_bundle),         intent(inout) :: rval

! Declare local variables
  integer(i_kind) j,ier,istatus
  integer(i_kind),dimension(nsig):: i1,i2,i3,i4
  integer(i_kind),dimension(nsig):: i1_b,i2_b,i3_b,i4_b
  integer(i_kind),dimension(nsig):: i1_a,i2_a,i3_a,i4_a
  real(r_kind) :: w1,w2,w3,w4
  real(r_kind) :: w1_b,w2_b,w3_b,w4_b
  real(r_kind) :: w1_a,w2_a,w3_a,w4_a
  real(r_kind) :: p_TL,p_AD,t_TL,t_AD,q_TL,q_AD
  real(r_kind) :: p_TL_b,t_TL_b,q_TL_b
  real(r_kind) :: p_TL_a,t_TL_a,q_TL_a
  real(r_kind) :: val,pg_gps
  real(r_kind) :: val_b,val_a
  real(r_kind) :: val_b_new,val_a_new
  real(r_kind) :: dz_b,dz_a
  real(r_kind) :: d_b,d_a
  real(r_kind) ::cg_gps,grad,p0,wnotgross,wgross
  real(r_kind) ::grad_a,grad_b
  real(r_kind),pointer,dimension(:) :: st,sq
  real(r_kind),pointer,dimension(:) :: rt,rq
  real(r_kind),pointer,dimension(:) :: sp
  real(r_kind),pointer,dimension(:) :: rp
  type(gpsrefNode), pointer :: gpsrefptr

!  If no gps obs return
  if(.not. associated(gpshead))return
! Retrieve pointers
! Simply return if any pointer not found
  ier=0
  call gsi_bundlegetpointer(sval,'tv'  ,st,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(sval,'q'   ,sq,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(sval,'prse',sp,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'tv'  ,rt,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'q'   ,rq,istatus);ier=istatus+ier
  call gsi_bundlegetpointer(rval,'prse',rp,istatus);ier=istatus+ier
  if(ier/=0)return

  !gpsptr => gpshead
  gpsrefptr => gpsrefNode_typecast(gpshead)
  do while (associated(gpsrefptr))

! Load location information into local variables
     do j=1,nsig
        i1(j)= gpsrefptr%ij(1,j)
        i2(j)= gpsrefptr%ij(2,j)
        i3(j)= gpsrefptr%ij(3,j)
        i4(j)= gpsrefptr%ij(4,j)
        ! lower obs (eyang)
        if (gpsrefptr%dhgt_b/=zero) then
           i1_b(j)= gpsrefptr%ij_b(1,j)
           i2_b(j)= gpsrefptr%ij_b(2,j)
           i3_b(j)= gpsrefptr%ij_b(3,j)
           i4_b(j)= gpsrefptr%ij_b(4,j)
        end if
        ! upper obs (eyang)
        if (gpsrefptr%dhgt/=zero) then
           i1_a(j)= gpsrefptr%ij_a(1,j)
           i2_a(j)= gpsrefptr%ij_a(2,j)
           i3_a(j)= gpsrefptr%ij_a(3,j)
           i4_a(j)= gpsrefptr%ij_a(4,j)
        end if
     enddo
     if (gpsrefptr%dhgt_b==zero) then
        print*, 'yeg intgpsref L179 no lower obs, i1_b(1)=',i1_b(1)
     else
        print*, 'yeg intgpsref L181 with lower obs, i1_b(1)=',i1_b(1)
     end if
     w1=gpsrefptr%wij(1)
     w2=gpsrefptr%wij(2)
     w3=gpsrefptr%wij(3)
     w4=gpsrefptr%wij(4)
     ! lower obs (eyang)
     if (gpsrefptr%dhgt_b/=zero) then
        w1_b=gpsrefptr%wij_b(1)
        w2_b=gpsrefptr%wij_b(2)
        w3_b=gpsrefptr%wij_b(3)
        w4_b=gpsrefptr%wij_b(4)
     end if
     ! upper obs (eyang)
     if (gpsrefptr%dhgt/=zero) then
        w1_a=gpsrefptr%wij_a(1)
        w2_a=gpsrefptr%wij_a(2)
        w3_a=gpsrefptr%wij_a(3)
        w4_a=gpsrefptr%wij_a(4)
     end if
     if (gpsrefptr%dhgt==zero) then
        print*, 'yeg intgpsref L202 no upper obs, w1_a=',w1_a
     else
        print*, 'yeg intgpsref L204 with upper obs, w1_a=',w1_a
     end if

     val=zero
     ! lower and upper obs (eyang)
     val_b=zero
     val_a=zero
     val_b_new=zero
     val_a_new=zero
     dz_b=zero
     dz_a=zero
     d_b=zero
     d_a=zero

!  local refractivity (linear operator)

     do j=1,nsig
        ! val = H_lin H_intp dx
        t_TL=w1* st(i1(j))+w2* st(i2(j))+w3* st(i3(j))+w4* st(i4(j))
        q_TL=w1* sq(i1(j))+w2* sq(i2(j))+w3* sq(i3(j))+w4* sq(i4(j))
        p_TL=w1* sp(i1(j))+w2* sp(i2(j))+w3* sp(i3(j))+w4* sp(i4(j))
        val = val + p_TL*gpsrefptr%jac_p(j) + t_TL*gpsrefptr%jac_t(j)+q_TL*gpsrefptr%jac_q(j)

        ! lower obs (eyang)
        ! val_b = H_lin_b H_intp_b dx
        if (gpsrefptr%dhgt_b/=zero) then
           t_TL_b=w1_b* st(i1_b(j))+w2_b* st(i2_b(j))+w3_b* st(i3_b(j))+w4_b* st(i4_b(j))
           q_TL_b=w1_b* sq(i1_b(j))+w2_b* sq(i2_b(j))+w3_b* sq(i3_b(j))+w4_b* sq(i4_b(j))
           p_TL_b=w1_b* sp(i1_b(j))+w2_b* sp(i2_b(j))+w3_b* sp(i3_b(j))+w4_b* sp(i4_b(j))
           val_b = val_b + p_TL_b*gpsrefptr%jac_p_b(j) + t_TL_b*gpsrefptr%jac_t_b(j)+q_TL_b*gpsrefptr%jac_q_b(j)
        end if

        ! upper obs (eyang)
        ! val_a = H_lin_a H_intp_a dx
        if (gpsrefptr%dhgt/=zero) then
           t_TL_a=w1_a* st(i1_a(j))+w2_a* st(i2_a(j))+w3_a* st(i3_a(j))+w4_a* st(i4_a(j))
           q_TL_a=w1_a* sq(i1_a(j))+w2_a* sq(i2_a(j))+w3_a* sq(i3_a(j))+w4_a* sq(i4_a(j))
           p_TL_a=w1_a* sp(i1_a(j))+w2_a* sp(i2_a(j))+w3_a* sp(i3_a(j))+w4_a* sp(i4_a(j))
           val_a = val_a + p_TL_a*gpsrefptr%jac_p_a(j) + t_TL_a*gpsrefptr%jac_t_a(j)+q_TL_a*gpsrefptr%jac_q_a(j)
        end if
     end do

     print*, 'yeg intgpsref L244 w1_b,i1_b(1)=',w1_b,i1_b(1)

     if (luse_obsdiag)then
        if (lsaveobsens) then
           grad = val*gpsrefptr%raterr2*gpsrefptr%err2
           !-- gpsptr%diags%obssen(jiter) = grad
           call obsdiagNode_set(gpsrefptr%diags,jiter=jiter,obssen=grad)
        else
           !-- if (gpsptr%luse) gpsptr%diags%tldepart(jiter)=val
           !if (gpsrefptr%luse) call obsdiagNode_set(gpsrefptr%diags,jiter=jiter,tldepart=val) ! commented out (eyang)

           ! val should be val_a minus val for gradient between ith and (i+1)th obs at ith obs (eyang) 
           ! ( (d2-d1) - (H2 dx - H1 dx) )*1/dz
           !-----------
           ! ( (H2 dx - H1 dx) )*1/dz
           !-----------
           if (gpsrefptr%dhgt/=zero) then ! if upper obs (i.e., there is ref grad. with upper obs)
              if (gpsrefptr%luse) call obsdiagNode_set(gpsrefptr%diags,jiter=jiter,tldepart=(val_a-val)/gpsrefptr%dhgt ) 
           else
              if (gpsrefptr%luse) call obsdiagNode_set(gpsrefptr%diags,jiter=jiter,tldepart=zero) 
           end if
        endif
     endif

!    Do adjoint
     if (l_do_adjoint) then

        if (.not. lsaveobsens) then
!           if( .not. ladtest_obs)  val=val-gpsptr%res
           !-------------------------------------
           print*, 'yeg intgpsref L263 ladtest_obs=',ladtest_obs
           if( .not. ladtest_obs) then ! gradient with below and above obs (eyang)
           ! val = val - gpsptr%res

              !---------------------------------------------
              ! val_a (upper obs) = H_lin3 H_intp3 dx3
              ! val   (middle obs)= H_lin2 H_intp2 dx2
              ! val_b (lower obs) = H_lin1 H_intp1 dx1
              !---------------------------------------------
              ! (Upper obs) val_a
              ! 1/(dz2)^2 ( H_lin2 H_intp2)^T [ (d3-d2) - (H_lin3 H_intp3 dx3 - H_lin2 H_intp2 dx2) ] sig_g2^(-2)
              !---------------------------------------------
              ! dz2 = (z3-z2) -> gpsptr%dhgt
              ! (d3 - d2) -> gpsptr%res
              if ( gpsrefptr%dhgt==zero ) then ! no gradient available (w/ upper obs)
                 val_a_new = zero
                 print*, "yeg intgpsref L278 wo upper obs, dhgt=",gpsrefptr%dhgt
              else
                 val_a_new = (one/( gpsrefptr%dhgt )**2) * ( gpsrefptr%res - (val_a - val) )
              end if
 
              !---------------------------------------------
              ! (Lower obs) val_b
              ! 1/(dz1)^2 (-H_lin2 H_intp2)^T [ (d2-d1) - (H_lin2 H_intp2 dx2 - H_lin1 H_intp1 dx1) ] sig_g1^(-2)
              !---------------------------------------------
              ! dz1 = (z2-z1) -> gpsptr%dhgt_b
              ! (d2 - d1) -> gpsptr%res_b
              ! val_b_new = (one/( gpsptr%dhgt_b )**2) * ( gpsptr%res_b - (val - val_b) )*(-1)
              if ( gpsrefptr%dhgt_b==zero ) then ! no gradient available (w/ lower obs)
                 val_b_new = zero
                 print*, "yeg intgpsref L292 wo lower obs, dhgt=",gpsrefptr%dhgt_b
              else
                 val_b_new = (one/( gpsrefptr%dhgt_b )**2) * ( (val - val_b) - gpsrefptr%res_b ) 
              end if

           end if
           !-------------------------------------
 
!          needed for gradient of nonlinear qc operator
           if (nlnqc_iter .and. gpsrefptr%pg > tiny_r_kind .and.  &
                                gpsrefptr%b  > tiny_r_kind) then
              pg_gps=gpsrefptr%pg*varqc_iter
              cg_gps=cg_term/gpsrefptr%b
              wnotgross= one-pg_gps
              wgross = pg_gps*cg_gps/wnotgross
              !p0   = wgross/(wgross+exp(-half*gpsrefptr%err2*val**2)) ! comment out (eyang)
              p0   = wgross/(wgross+exp(-half*gpsrefptr%err2*(val_a_new*gpsrefptr%dhgt)**2)) ! val_a_new (grad with upper obs)
              !val = val*(one-p0) ! commented out (eyang)
              val = (val_a_new*gpsrefptr%dhgt)*(one-p0)
              print*, 'yeg intgpsref L322, nlnqc_iter,val=',nlnqc_iter,val
           endif
       
           if( ladtest_obs) then
              grad = val
           else
!              grad = val*gpsptr%raterr2*gpsptr%err2
              !-------------------------------------
              ! For refractivity gradient (eyang)
              !-------------------------------------
              grad_a = val_a_new*gpsrefptr%raterr2*gpsrefptr%err2
              grad_b = val_b_new*gpsrefptr%raterr2_b*gpsrefptr%err2_b
              grad = grad_a + grad_b
              print*, 'yeg intgpsref L319 grad_a,grad_b,grad=',grad_a,grad_b,grad

           end if
        endif


!       adjoint 

        ! (H_lin2 H_intp2)^T = (H_intp2)^T (H_lin2)^T
        do j=1,nsig

           t_AD = grad*gpsrefptr%jac_t(j)
           rt(i1(j))=rt(i1(j))+w1*t_AD
           rt(i2(j))=rt(i2(j))+w2*t_AD
           rt(i3(j))=rt(i3(j))+w3*t_AD
           rt(i4(j))=rt(i4(j))+w4*t_AD
           q_AD = grad*gpsrefptr%jac_q(j)
           rq(i1(j))=rq(i1(j))+w1*q_AD
           rq(i2(j))=rq(i2(j))+w2*q_AD
           rq(i3(j))=rq(i3(j))+w3*q_AD
           rq(i4(j))=rq(i4(j))+w4*q_AD
           p_AD = grad*gpsrefptr%jac_p(j)
           rp(i1(j))=rp(i1(j))+w1*p_AD
           rp(i2(j))=rp(i2(j))+w2*p_AD
           rp(i3(j))=rp(i3(j))+w3*p_AD
           rp(i4(j))=rp(i4(j))+w4*p_AD
        
        enddo

     endif

     !gpsptr => gpsptr%llpoint
     gpsrefptr => gpsrefNode_nextcast(gpsrefptr)

  end do

  return
end subroutine intgpsref_


end module intgpsrefmod
