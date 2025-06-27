module stpgpsrefmod

!$$$ module documentation block
!           .      .    .                                       .
! module:   stpgpsrefmod    module for stpref and its tangent linear stpref_tl
!  prgmmr:
!
! abstract: module for stpref and its tangent linear stpref_tl
!
! program history log:
!   2005-05-19  Yanqiu zhu - wrap stpref and its tangent linear stpref_tl into one module
!   2005-11-16  Derber - remove interfaces
!   2009-08-12  lueken - updated documentation
!   2010-05-13  todling - uniform interface across stp routines
!   2013-10-28  todling - rename p3d to prse
!   2016-05-18  guo     - replaced ob_type with polymorphic obsNode through type casting
!   2025-02-18  eyang   - modified stpgps to stpgpsref
!
! subroutines included:
!   sub stpgpsref
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

implicit none

PRIVATE
PUBLIC stpgpsref

contains

subroutine stpgpsref(gpshead,rval,sval,out,sges,nstep)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram: stpref    compute contribution to penalty and stepsize
!                       from ref, using nonlinear qc    
!   prgmmr: cucurull,l.     org: JCSDA/NCEP           date: 2004-05-06
!
! abstract:  This routine applies the (linear) operator for local 
!            refractivity and linear linear estimate for step size. 
!            This version includes nonlinear qc.
!
! program history log:
!   2004-05-06  cucurull 
!   2004-06-21  treadon - update documentation
!   2004-08-02  treadon - add only to module use, add intent in/out
!   2004=10-08  parrish - add nonlinear qc option
!   2004-11-30  cucurull- add increments for surface pressure and temperature at levels
!                         below observation. Install non-linear forward operator.
!   2005-01-26  cucurull- Implement local GPS RO operator
!   2005-03-23  cucurull- correct bouds for obs below the second level
!   2005-04-11  treadon - merge stpref and stpref_qc into single routine
!   2005-08-02  derber  - modify for variational qc parameters for each ob
!   2005-09-28  derber  - consolidate location and weight arrays
!   2005-12-02  cucurull - fix bug for dimensions of sp and rp
!   2007-03-19  tremolet - binning of observations
!   2007-07-28  derber  - modify to use new inner loop obs data structure
!                       - unify NL qc
!   2006-09-06  cucurull - generalize code to hybrid vertical coordinate and modify to use 
!                          surface pressure
!   2007-01-13  derber - clean up code and use coding standards
!   2007-06-04  derber  - use quad precision to get reproducability over number of processors
!   2007-07-26  cucurull - input 3d pressure to update code to generalized vertical coordinate
!   2008-04-11  safford - rm unused vars
!   2008-12-03  todling - changed handling of ptr%time
!   2010-01-04  zhang,b - bug fix: accumulate penalty for multiple obs bins
!   2010-05-13  todling - update to use gsi_bundle
!
!   input argument list:
!     gpshead
!     rt    - search direction (gradxJ) for virtual temperature
!     rq    - search direction (gradxJ) for specific humidity
!     rp    - search direction (gradxJ) for (3D) pressure
!     st    - analysis increment (correction) for virtual temperature
!     sq    - analysis increment (correction) for specific humidity
!     sp    - analysis increment (correction) for (3D) pressure
!     sges  - stepsize estimates (nstep)
!     nstep - number of stepsize estimates (==0 means use outer iteration values)
!                                         
!   output argument list:
!     out(1:nstep)- contribution to penalty from local gps refractivity - sges(1:nstep)
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
  use kinds, only: r_kind,i_kind,r_quad
  use qcmod, only: nlnqc_iter,varqc_iter
  use constants, only: zero,one,two,half,tiny_r_kind,cg_term,zero_quad,r3600
  use gridmod, only: nsig
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlegetpointer
  use m_obsNode, only: obsNode
  use m_gpsrefNode, only: gpsrefNode
  use m_gpsrefNode, only: gpsrefNode_typecast
  use m_gpsrefNode, only: gpsrefNode_nextcast
  implicit none

! Declare passed variables
  class(obsNode   ),pointer           ,intent(in   ) :: gpshead
  integer(i_kind)                     ,intent(in   ) :: nstep
  real(r_quad),dimension(max(1,nstep)),intent(inout) :: out
  type(gsi_bundle)                    ,intent(in   ) :: rval,sval
  real(r_kind),dimension(max(1,nstep)),intent(in   ) :: sges

! Declare local variables
  integer(i_kind) j,kk,ier,istatus
  integer(i_kind),dimension(nsig):: i1,i2,i3,i4
  integer(i_kind),dimension(nsig):: i1_a,i2_a,i3_a,i4_a
  real(r_kind) :: val,val2
  real(r_kind) :: w1,w2,w3,w4
  real(r_kind) :: w1_a,w2_a,w3_a,w4_a
  real(r_kind) :: q_TL,p_TL,t_TL
  real(r_kind) :: rq_TL,rp_TL,rt_TL
  real(r_kind) :: q_TL_a,p_TL_a,t_TL_a
  real(r_kind) :: rq_TL_a,rp_TL_a,rt_TL_a
  real(r_kind),pointer,dimension(:) :: st,sq
  real(r_kind),pointer,dimension(:) :: rt,rq
  real(r_kind),pointer,dimension(:) :: sp
  real(r_kind),pointer,dimension(:) :: rp
  type(gpsrefNode), pointer :: gpsrefptr

  real(r_kind) cg_gps,wgross,wnotgross
  real(r_kind) pg_gps,nref
  real(r_kind),dimension(max(1,nstep))::pen

! Initialize penalty, b1, and b3 to zero
  out=zero_quad

!  If no gps data return
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


! Loop over observations
  gpsrefptr => gpsrefNode_typecast(gpshead)
  do while (associated(gpsrefptr))

     if(gpsrefptr%luse)then

        !val2=-gpsptr%res
        !d2-d1
        val2=-gpsrefptr%res
        if(nstep > 0)then
           do j=1,nsig
              i1(j)= gpsrefptr%ij(1,j)
              i2(j)= gpsrefptr%ij(2,j)
              i3(j)= gpsrefptr%ij(3,j)
              i4(j)= gpsrefptr%ij(4,j)
              ! upper obs (eyang)
              i1_a(j)= gpsrefptr%ij_a(1,j)
              i2_a(j)= gpsrefptr%ij_a(2,j)
              i3_a(j)= gpsrefptr%ij_a(3,j)
              i4_a(j)= gpsrefptr%ij_a(4,j)
           enddo
           w1=gpsrefptr%wij(1)
           w2=gpsrefptr%wij(2)
           w3=gpsrefptr%wij(3)
           w4=gpsrefptr%wij(4)
           ! upper obs (eyang)
           w1_a=gpsrefptr%wij_a(1)
           w2_a=gpsrefptr%wij_a(2)
           w3_a=gpsrefptr%wij_a(3)
           w4_a=gpsrefptr%wij_a(4)

           val=zero

           do j=1,nsig
              t_TL =w1* st(i1(j))+w2* st(i2(j))+w3* st(i3(j))+w4* st(i4(j))
              rt_TL=w1* rt(i1(j))+w2* rt(i2(j))+w3* rt(i3(j))+w4* rt(i4(j))
              q_TL =w1* sq(i1(j))+w2* sq(i2(j))+w3* sq(i3(j))+w4* sq(i4(j))
              rq_TL=w1* rq(i1(j))+w2* rq(i2(j))+w3* rq(i3(j))+w4* rq(i4(j))
              p_TL =w1* sp(i1(j))+w2* sp(i2(j))+w3* sp(i3(j))+w4* sp(i4(j))
              rp_TL=w1* rp(i1(j))+w2* rp(i2(j))+w3* rp(i3(j))+w4* rp(i4(j))
              ! upper obs (eyang)
              t_TL_a =w1_a* st(i1_a(j))+w2_a* st(i2_a(j))+w3_a* st(i3_a(j))+w4_a* st(i4_a(j))
              rt_TL_a=w1_a* rt(i1_a(j))+w2_a* rt(i2_a(j))+w3_a* rt(i3_a(j))+w4_a* rt(i4_a(j))
              q_TL_a =w1_a* sq(i1_a(j))+w2_a* sq(i2_a(j))+w3_a* sq(i3_a(j))+w4_a* sq(i4_a(j))
              rq_TL_a=w1_a* rq(i1_a(j))+w2_a* rq(i2_a(j))+w3_a* rq(i3_a(j))+w4_a* rq(i4_a(j))
              p_TL_a =w1_a* sp(i1_a(j))+w2_a* sp(i2_a(j))+w3_a* sp(i3_a(j))+w4_a* sp(i4_a(j))
              rp_TL_a=w1_a* rp(i1_a(j))+w2_a* rp(i2_a(j))+w3_a* rp(i3_a(j))+w4_a* rp(i4_a(j))

              !----------------------
              ! val2 and val (eyang)
              ! why not for lower obs?
              ! how can we incorporate lower obs in this format? (EYANG)
              !----------------------
!              val2 = val2 + t_tl*gpsptr%jac_t(j)+ q_tl*gpsptr%jac_q(j)+p_tl*gpsptr%jac_p(j) 
!              val  = val + rt_tl*gpsptr%jac_t(j)+rq_tl*gpsptr%jac_q(j)+rp_tl*gpsptr%jac_p(j)

              ! val2 = -gpsptr%res = -(d2-d1)
              ! val2 <=  -(d2-d1) + (H_lin2 H_intp2 dx2 - H_lin1 H_intp1 dx1)
              if (gpsrefptr%dhgt==zero) then ! no gradient available (w/ upper obs)
                 !===================NEED TO TAKE CARE OF (yzhu)
                 val2 = zero ! is it right? yzhu
                 val = zero
                 print*, 'yeg stpgpsref L208 NoGradientUpper dght=',gpsrefptr%dhgt
              else
                 val2 = val2 + t_tl_a*gpsrefptr%jac_t_a(j)+ q_tl_a*gpsrefptr%jac_q_a(j)+p_tl_a*gpsrefptr%jac_p_a(j) & 
                             - t_tl  *gpsrefptr%jac_t(j)  - q_tl  *gpsrefptr%jac_q(j)  -p_tl  *gpsrefptr%jac_p(j) 

                 val  = val + rt_tl_a*gpsrefptr%jac_t_a(j)+rq_tl_a*gpsrefptr%jac_q_a(j)+rp_tl_a*gpsrefptr%jac_p_a(j) &
                            - rt_tl  *gpsrefptr%jac_t(j)  -rq_tl  *gpsrefptr%jac_q(j)  -rp_tl  *gpsrefptr%jac_p(j)
                 print*, 'yeg stpgpsref L216: gpsrefptr%dhgt=',gpsrefptr%dhgt
                 print*, 'yeg stpgpsref L217: val,val2=',val,val2
              end if

           enddo


!          penalty and gradient

           do kk=1,nstep
              nref=val2+sges(kk)*val
              pen(kk)=nref*nref*gpsrefptr%err2
              print*, 'stpgpsref L228: kk,val,val2=',kk,val,val2
              print*, 'stpgpsref L229: nref,err2=',nref,gpsrefptr%err2
              print*, 'stpgpsref L230: pen(kk)=',pen(kk)
           end do
        else
           pen(1)=val2*val2*gpsrefptr%err2
        end if

!       Modify penalty term if nonlinear QC
        if (nlnqc_iter .and. gpsrefptr%pg > tiny_r_kind .and. gpsrefptr%b > tiny_r_kind) then
           pg_gps=gpsrefptr%pg*varqc_iter
           cg_gps=cg_term/gpsrefptr%b
           wnotgross= one-pg_gps
           wgross = pg_gps*cg_gps/wnotgross
           do kk=1,max(1,nstep)
              pen(kk) = -two*log((exp(-half*pen(kk)) + wgross)/(one+wgross))
           end do
        endif
  
!       Cost function
        out(1) = out(1)+pen(1)*gpsrefptr%raterr2
        do kk=2,nstep
           out(kk) = out(kk)+(pen(kk)-pen(1))*gpsrefptr%raterr2
        end do

     endif
  
     gpsrefptr => gpsrefNode_nextcast(gpsrefptr)

  end do

 return
end subroutine stpgpsref


end module stpgpsrefmod
