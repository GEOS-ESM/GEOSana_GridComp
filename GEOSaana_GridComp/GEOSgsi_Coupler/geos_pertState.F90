!#define DEBUG_TRACE
!#define DEBUG_VERBOSE
!#define DBLE_STATE
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  geos_pertState --- redistributions between GSI and aGCM
!
! !INTERFACE:
!
  module geos_pertState

! !USES:


  use kinds,       only : r_kind,r_single,i_kind
#ifdef DBLE_STATE
  use ESMF,    only : r_ESkind => ESMF_KIND_R8
#else
  use ESMF,    only : r_ESkind => ESMF_KIND_R4
#endif
  use mpimod,      only : mype,mpi_rtype,mpi_real16,mpi_comm_world
  use gridmod,     only : strip
  use gridmod,     only : displs_s,ijn_s
  use gridmod,     only : nlat, nlon     ! no. lat/lon
  use gridmod,     only : lat1, lon1     ! no. lat/lon on subdomain (no buffer)
  use gridmod,     only : lat2, lon2     ! no. lat/lon on subdomain (buffer pnts on ends)
  use gridmod,     only : nsig           ! no. levels
  use gridmod,     only : iglobal        ! no. of horizontal points on global grid
  use gridmod,     only : ijn            ! no. of horiz. pnts for each subdomain (no buffer)
  use gridmod,     only : displs_g       ! comm. array, displacement for receive on global grid
  use gridmod,     only : itotsub        ! no. of horizontal points of all subdomains combined
  use gridmod,     only : bk5         
  use gsi_bias,    only : reorder21,reorder12
  use constants,   only : zero,one,tiny_r_kind
  use state_vectors,only: dot_product

  use m_interpack,    only : interpack_terpv
  use m_interpack_ad, only : interpack_terpv_ad

  use mpeu_util, only: tell,perr,die

      implicit none

! !PUBLIC MEMBER FUNCTIONS:

      PRIVATE

  public:: pertState_set	! set a scalar or vector to a pertState
  public:: pertState_show	! display contents of a pertState

  public:: pertState_setNeeded	! set variables in a state (of export) needed
  public:: pertState_getIFx	! get a pertState from a gsiPert vector, only if a
  				! variable is required by the gsiPert vector.

  public:: pertState_putIFx	! put a pertState into a gsiPert vector, only if a
  				! variable is present in the gsiPert vector.

  public:: pertState_dot_product, &
	   dot_product		! dot product of 2 states or of 1 or 2 variables

  public:: TLMvect		! "which" values used internally to flag the
  public:: ADMvect		! vector class, thus the operation mode.

  character(len=*),parameter:: TLMvect = "tlm"
  character(len=*),parameter:: ADMvect = "adm"

!  public:: pertState_create	! (s,x)		create s based on x
!  public:: pertState_add	! (s,d,a)	add values of d to s, s+=[a*]d
!  public:: pertState_copy	! (x,s)		copy values of x to s
!  public:: pertState_destroy	! (s)		destroy s

! !PUBLIC DATA:

!
! !DESCRIPTION: Maps gsi increments on to gcm perturbations and vice-versa.
!
! !REMARKS:
!
!   1) This package assumes xpert to have A-grid winds in its u/v slot
!   2) This package assumes xpert to have TV in the its pt slot, unless when
!      perturbation vector in read in from file.
!
! !TO DO:
!
!   1) Need to work on O3 (ADM/TLM not ready for this)
!   2) This interface should ultimately become an ESMF-Grided Component
!      and therefore become non-specific, i.e., applicable to an GCM
!      TL and AD models. For now, this is specific to GEOS-5.  
!   3) Allow vectype_ to be reset from RC file specific to this code to
!      handle GEOS-4 perturbation vector.
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  01Apr2010  Treadon   move strip to gridmod
!  15May2010  Todling   Update to use GSI_Bundle
!  31Aug2010  Guo       Added internal_state variables.
!  20May2011  Guo       This version is workable with stub_GEOSagcmPert,
!			and it has be verified against a GSI built-in
!			adjoint-test in GSI evaljgrad().
!  28Oct2013  Todling   Rename p3d to prse
!
!EOP
!-------------------------------------------------------------------------

  interface pertState_set; module procedure &
  	set2avalue_,	&
	set3d0_    ,	&
  	set3dv_    ,	&
  	set2astate_; end interface

  interface pertState_show; module procedure &
  	show_; end interface

  interface pertState_setNeeded; module procedure &
  	setNeeded_; end interface

  interface pertState_dot_product; module procedure &
	dotprod_st_    , &
	dotprod_2st_   , &
	dotprod_vmask_ , &
	dotprod_2vmask_, &
	dotprod_var_   , &
  	dotprod_2var_  ; end interface
  interface dot_product; module procedure &
	dotprod_st_    , &
	dotprod_2st_   , &
	dotprod_vmask_ , &
	dotprod_2vmask_, &
	dotprod_var_   , &
  	dotprod_2var_  ; end interface

  interface pertState_getIFx; module procedure &
  	getIFx_; end interface
  interface pertState_putIFx; module procedure &
  	putIFx_; end interface

  character(len=*), parameter:: myname = 'geos_pertState'

  	! These are expected from a GSI_bundle
  character(len=*),parameter::    u_GSIpert='u'		! = U
  character(len=*),parameter::    v_GSIpert='v'		! = V
  character(len=*),parameter::   ps_GSIpert='ps'	! = sum(DP)_all + p_top
  character(len=*),parameter:: prse_GSIpert='prse'	! = sum(DP)_k
  character(len=*),parameter::   tv_GSIpert='tv'	! = TV
  character(len=*),parameter:: tsen_GSIpert='tsen'	! = ?
  character(len=*),parameter::    q_GSIpert='q'		! = QV
  character(len=*),parameter::   oz_GSIpert='oz'	! = O3, in GG=gram/gram
  character(len=*),parameter::   cw_GSIpert='cw'	! = QL + QI
  character(len=*),parameter::  sst_GSIpert='sst'	! = ?

  character(len=*),dimension(8),parameter:: var3dList_GSIpert = &
    (/ u_GSIpert,   v_GSIpert,prse_GSIpert,  tv_GSIpert, &
    tsen_GSIpert,   q_GSIpert,  oz_GSIpert,  cw_GSIpert /)

  character(len=*),dimension(2),parameter:: var2dList_GSIpert = &
    (/ps_GSIpert, sst_GSIpert /)

  	! These are expected from a PGCM state
  character(len=*),parameter::    u_GCMpert='U'		! = u
  character(len=*),parameter::    v_GCMpert='V'		! = v
  character(len=*),parameter::   tv_GCMpert='TV'	! = tv
  character(len=*),parameter::   dp_GCMpert='DP'	! = d(prse,ps)/dk
  character(len=*),parameter::   qv_GCMpert='QV'	! = q
  character(len=*),parameter::   ql_GCMpert='QL'	! = cw - qi
  character(len=*),parameter::   qi_GCMpert='QI'	! = qi (=0)
  character(len=*),parameter::   o3_GCMpert='O3'	! = oz, in PPMV

  character(len=*),dimension(8),parameter:: var3dList_GCMpert = &
    (/ u_GCMpert,  v_GCMpert, tv_GCMpert, dp_GCMpert, &
      qv_GCMpert, ql_GCMpert, qi_GCMpert, o3_GCMpert  /)

   integer(i_kind),parameter:: ROOT =0
   real(r_kind), parameter :: GG_per_PPMV= 1.657E-6_r_kind ! GSI Ozone is in GG
   real(r_kind), parameter :: kPa_per_Pa = 0.001_r_kind
   real(r_kind), parameter :: Pa_per_kPa = 1000._r_kind
   integer(i_kind),parameter:: GENERIC_ERROR = -9999

#include "MAPL_ErrLog.h"
#include "mytrace.H"

contains
function dotprod_st_(xst,hw) result(dot_)
  use ESMF, only: ESMF_State
  use kinds, only: r_quad,i_kind
  implicit none
  type(ESMF_State),intent(in) :: xst
  integer(i_kind),optional,intent(in):: hw	! halowidth
  real(r_quad):: dot_
  dot_=0._r_quad
  dot_=dot_+dotprod_var_(xst,xst, u_GCMpert,hw=hw)
  dot_=dot_+dotprod_var_(xst,xst, v_GCMpert,hw=hw)
  dot_=dot_+dotprod_var_(xst,xst,tv_GCMpert,hw=hw)
  dot_=dot_+dotprod_var_(xst,xst,dp_GCMpert,hw=hw)
  dot_=dot_+dotprod_var_(xst,xst,qv_GCMpert,hw=hw)
  dot_=dot_+dotprod_var_(xst,xst,ql_GCMpert,hw=hw)
  dot_=dot_+dotprod_var_(xst,xst,qi_GCMpert,hw=hw)
  dot_=dot_+dotprod_var_(xst,xst,o3_GCMpert,hw=hw)
end function dotprod_st_
function dotprod_2st_(xst,yst,hw) result(dot_)
  use ESMF, only: ESMF_State
  use kinds, only: r_quad,i_kind
  implicit none
  type(ESMF_State),intent(in) :: xst
  type(ESMF_State),intent(in) :: yst
  integer(i_kind),optional,intent(in):: hw	! halowidth
  real(r_quad):: dot_
  dot_=0._r_quad
  dot_=dot_+dotprod_var_(xst,yst, u_GCMpert,hw=hw)
  dot_=dot_+dotprod_var_(xst,yst, v_GCMpert,hw=hw)
  dot_=dot_+dotprod_var_(xst,yst,tv_GCMpert,hw=hw)
  dot_=dot_+dotprod_var_(xst,yst,dp_GCMpert,hw=hw)
  dot_=dot_+dotprod_var_(xst,yst,qv_GCMpert,hw=hw)
  dot_=dot_+dotprod_var_(xst,yst,ql_GCMpert,hw=hw)
  dot_=dot_+dotprod_var_(xst,yst,qi_GCMpert,hw=hw)
  dot_=dot_+dotprod_var_(xst,yst,o3_GCMpert,hw=hw)
end function dotprod_2st_

function dotprod_var_(xst,yst,vnm,hw) result(dot_)
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_SUCCESS
  use MAPL_Mod, only: ESMFL_StateGetPointerToData
  use kinds, only: r_single,r_quad,i_kind
  implicit none
  type(ESMF_State),target,intent(in):: xst
  type(ESMF_State),target,intent(in):: yst
  character(len=*),intent(in):: vnm
  integer(i_kind),optional,intent(in):: hw	! halowidth
  real(r_quad):: dot_

  type(ESMF_State),pointer:: p_xst,p_yst
  real(r_ESkind),pointer,dimension(:,:,:):: p_xvar,p_yvar
  real(r_quad):: x,y
  logical:: sameshape
  integer:: ier
  character(len=*), parameter :: myname_ = myname//'::dotprod_var_'
#ifdef _MYNAME_
#undef _MYNAME_
#endif
#define _MYNAME_ myname_//'("'//trim(vnm)//'")'
_ENTRY_(_MYNAME_)

  p_xst => xst
  p_xvar=> null()
  call ESMFL_StateGetPointerToData(p_xst,p_xvar,vnm,rc=ier)
		if(ier/=ESMF_SUCCESS) &
		  call die(myname_,'ESMFL_StateGetPointerToData(xst,"'//trim(vnm)//'"), rc =',ier)

		if(.not.associated(p_xvar)) &
		  call die(myname_,'.not.associated(ptr => p_xst%"'//trim(vnm)//'")')

  p_yst => yst
  p_yvar => null()
  call ESMFL_StateGetPointerToData(p_yst,p_yvar,vnm,rc=ier)
		if(ier/=ESMF_SUCCESS) &
		  call die(myname_,'ESMFL_StateGetPointerToData(yst,"'//trim(vnm)//'"), rc =',ier)

		if(.not.associated(p_yvar)) &
		  call die(myname_,'.not.associated(ptr => p_yst%"'//trim(vnm)//'")')


  sameshape = all(shape(p_xvar)==shape(p_yvar))
  if(.not.sameshape) then
    call perr(myname_,'mismatched shape, xst("'//trim(vnm)//'").vs.yst("'//trim(vnm)//'"')
    call perr(myname_,'            shape(xst("'//trim(vnm)//'") =',shape(p_xvar))
    call perr(myname_,'            shape(yst("'//trim(vnm)//'") =',shape(p_yvar))
    call die(myname_)
  endif

  dot_=par_dotprod3d_(p_xvar,p_yvar,hw=hw)
  nullify(p_xst,p_yst)
  nullify(p_xvar,p_yvar)
_EXIT_(_MYNAME_)
end function dotprod_var_

function dotprod_vmask_(xst,yst,vmask,hw) result(dot_)
  use ESMF, only: ESMF_State
  use kinds, only: r_quad,i_kind
  use mpeu_util, only: tell,perr,die
  implicit none
  type(ESMF_State),target,intent(in):: xst
  type(ESMF_State),target,intent(in):: yst
  character(len=*),dimension(:),intent(in):: vmask
  integer(i_kind),optional,intent(in):: hw	! halowidth
  real(r_quad):: dot_

  character(len=*), parameter :: myname_ = myname//'::dotprod_vmask_'
  integer(i_kind):: iv
_ENTRY_(myname_)

  dot_=0.
  do iv=1,size(vmask)
    dot_=dot_+dotprod_var_(xst,yst,vmask(iv),hw=hw)
  enddo

_EXIT_(myname_)
end function dotprod_vmask_

function dotprod_2var_(xst,xnm,yst,ynm,hw) result(dot_)
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_SUCCESS
  use MAPL_Mod, only: ESMFL_StateGetPointerToData
  use kinds, only: r_single,r_quad,i_kind
  implicit none
  type(ESMF_State),target,intent(in):: xst
  character(len=*),intent(in):: xnm
  type(ESMF_State),target,intent(in):: yst
  character(len=*),intent(in):: ynm
  integer(i_kind),optional,intent(in):: hw	! halowidth
  real(r_quad):: dot_

  real(r_ESkind),pointer,dimension(:,:,:):: p_xvar,p_yvar
  type(ESMF_State),pointer:: p_xst,p_yst
  logical:: sameshape
  integer(i_kind):: ier
  character(len=*), parameter :: myname_ = myname//'::dotprod_2var_'
#ifdef _MYNAME_
#undef _MYNAME_
#endif
#define _MYNAME_ myname_//'("'//trim(xnm)//'", "'//trim(ynm)//'")'
_ENTRY_(_MYNAME_)

  p_xst => xst
  p_xvar=> null()
  call ESMFL_StateGetPointerToData(p_xst,p_xvar,xnm,rc=ier)
  	if(ier/=ESMF_SUCCESS) &
	  call die(myname_,'ESMFL_StateGetPointerToData(xst,"'//trim(xnm)//'"), rc =',ier)
  p_yst => yst
  p_yvar => null()
  call ESMFL_StateGetPointerToData(p_yst,p_yvar,ynm,rc=ier)
  	if(ier/=ESMF_SUCCESS) &
	  call die(myname_,'ESMFL_StateGetPointerToData(yst,"'//trim(ynm)//'"), rc =',ier)

  sameshape = all(shape(p_xvar)==shape(p_yvar))
  if(.not.sameshape) then
    call perr(myname_,'mismatched shape, xst("'//trim(xnm)//'").vs.yst("'//trim(ynm)//'"')
    call perr(myname_,'            shape(xst("'//trim(xnm)//'") =',shape(p_xvar))
    call perr(myname_,'            shape(yst("'//trim(ynm)//'") =',shape(p_yvar))
    call die(myname_)
  endif

  dot_=par_dotprod3d_(p_xvar,p_yvar,hw=hw)
  nullify(p_xst,p_yst)
  nullify(p_xvar,p_yvar)
_EXIT_(_MYNAME_)
end function dotprod_2var_

function dotprod_2vmask_(xst,xmask,yst,ymask,hw) result(dot_)
  use ESMF, only: ESMF_State
  use kinds, only: r_quad,i_kind
  use mpeu_util, only: tell,perr,die
  implicit none
  type(ESMF_State),target,intent(in):: xst
  character(len=*),dimension(:),intent(in):: xmask
  type(ESMF_State),target,intent(in):: yst
  character(len=*),dimension(:),intent(in):: ymask
  integer(i_kind),optional,intent(in):: hw	! halowidth
  real(r_quad):: dot_

  character(len=*), parameter :: myname_ = myname//'::dotprod_2vmask_'
  integer(i_kind):: iv
_ENTRY_(myname_)

  if(size(xmask)/=size(ymask)) then
    call perr(myname_,'size(xmask)/=size(ymask)')
    call perr(myname_,'size(xmask) =',size(xmask))
    call perr(myname_,'size(ymask) =',size(ymask))
    call die(myname_)
  endif

  dot_=0.
  do iv=1,size(xmask)
    dot_=dot_+dotprod_2var_(xst,xmask(iv),yst,ymask(iv),hw=hw)
  enddo

_EXIT_(myname_)
end function dotprod_2vmask_

function par_dotprod3d_(x,y,hw) result(dot_)
  use kinds, only: r_single, r_quad
  use mpeu_mpif, only: mpi_real16
  use mpimod, only: mpi_comm_world
  use mpeu_util, only: die
  use mpl_allreducemod, only: mpl_allreduce
  implicit none
  real(r_ESkind),dimension(:,:,:),intent(in):: x,y
  integer(i_kind),optional,intent(in):: hw	! halowidth
  real(r_quad  ):: dot_

  integer(i_kind):: i,j,k, hw_
  real(r_quad):: x_,y_,dp_(1)
  character(len=*), parameter :: myname_ = myname//'::par_dotprod3d_'
_ENTRY_(myname_)

  hw_=0
  if(present(hw)) hw_=hw

  dp_(1)=0._r_quad
  do k=1,size(x,3)
  do j=1+hw_,size(x,2)-hw_
  do i=1+hw_,size(x,1)-hw_
    x_=real(x(i,j,k),kind=r_quad)
    y_=real(y(i,j,k),kind=r_quad)
    dp_(1)=dp_(1)+x_*y_
  enddo
  enddo
  enddo

  call mpl_allreduce(1,dp_(:))
  dot_=dp_(1)
_EXIT_(myname_)
end function par_dotprod3d_

subroutine set2avalue_(pert,v,rc)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 subroutine set2avalue_
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:	 2010-11-04
!
! abstract: set a value (i.e. 0.0) to a pertState
!
! program history log:
!   2010-11-04  j guo   - added this document block
!
!   input argument list: see Fortran 90 style document below
!
!   output argument list: see Fortran 90 style document below
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:
!
!$$$  end subprogram documentation block

! subroutine interface:

  use ESMF, only: ESMF_State
  use gsi_bundlemod, only: gsi_bundle
!  use kinds, only: i_kind,r_kind
!  use mpeu_util, only: tell,perr,die
  implicit none

  type(ESMF_State) ,intent(inout):: pert	! output to a pre-configged state
  real(kind=r_kind),intent(in   ):: v	 	! a value
  integer(kind=i_kind),optional,intent(  out):: rc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname_=myname//'::set2avalue_'
  integer(i_kind):: ier

_ENTRY_(myname_)
  if(present(rc)) rc=0

	! Get import variables, to initialize them

! pert%u = v
  call set3d0_(pert,u_GCMpert,v,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3d0_("'//	&
		trim(u_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

! pert%v = v
  call set3d0_(pert,v_GCMpert,v,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3d0_("'//	&
		trim(v_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

! pert%tv = v
  call set3d0_(pert,tv_GCMpert,v,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3d0_("'//	&
		trim(tv_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

! pert%qv = v
  call set3d0_(pert,qv_GCMpert,v,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3d0_("'//	&
		trim(qv_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

! pert%o3 = v
  call set3d0_(pert,o3_GCMpert,v,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3d0_("'//	&
		trim(o3_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

! pert%ql = v
  call set3d0_(pert,ql_GCMpert,v,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3d0_("'//	&
		trim(ql_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

! pert%qi = v
  call set3d0_(pert,qi_GCMpert,v,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3d0_("'//	&
		trim(qi_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

! pert%dp = v
  call set3d0_(pert,dp_GCMpert,v,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3d0_("'//	&
		trim(dp_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

_EXIT_(myname_)
end subroutine set2avalue_
subroutine set2astate_(spert,rpert,rc)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 subroutine set2astate_
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:	 2010-11-04
!
! abstract: set a value (i.e. 0.0) to a pertState
!
! program history log:
!   2010-11-04  j guo   - added this document block
!
!   input argument list: see Fortran 90 style document below
!
!   output argument list: see Fortran 90 style document below
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:
!
!$$$  end subprogram documentation block

! subroutine interface:

  use ESMF, only: ESMF_State
  use gsi_bundlemod, only: gsi_bundle
  use kinds, only: i_kind
  use mpeu_util, only: tell,perr,die
  implicit none

  type(ESMF_State) ,target,intent(inout):: spert	! output to a pre-configged state
  type(ESMF_State) ,target,intent(in   ):: rpert	! the RHS input state
  integer(i_kind),optional,intent(  out):: rc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname_=myname//'::set2astate_'
  integer(i_kind):: ier

_ENTRY_(myname_)
  if(present(rc)) rc=0

	! Get import variables, to initialize them

! pert%u = v
  call set3dv_(spert,u_GCMpert,rpert,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3dv_("'//	&
		trim(u_GCMpert)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

! pert%v = v
  call set3dv_(spert,v_GCMpert,rpert,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3dv_("'//	&
		trim(v_GCMpert)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

! pert%tv = v
  call set3dv_(spert,tv_GCMpert,rpert,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3dv_("'//	&
		trim(tv_GCMpert)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

! pert%qv = v
  call set3dv_(spert,qv_GCMpert,rpert,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3dv_("'//	&
		trim(qv_GCMpert)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

! pert%o3 = v
  call set3dv_(spert,o3_GCMpert,rpert,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3dv_("'//	&
		trim(o3_GCMpert)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

! pert%ql = v
  call set3dv_(spert,ql_GCMpert,rpert,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3dv_("'//	&
		trim(ql_GCMpert)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

! pert%qi = v
  call set3dv_(spert,qi_GCMpert,rpert,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3dv_("'//	&
		trim(qi_GCMpert)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

! pert%dp = v
  call set3dv_(spert,dp_GCMpert,rpert,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'set3dv_("'//	&
		trim(dp_GCMpert)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

_EXIT_(myname_)
end subroutine set2astate_

subroutine show_(spert,who,which,rc)
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_SUCCESS
  use mpeu_util, only: perr,die
  use kinds    , only: i_kind
  implicit none
  type(ESMF_State),intent(inout):: spert
  character(len=*),intent(in):: who,which
  integer(i_kind),optional,intent(out):: rc

  integer(i_kind):: ier
  character(len=*),parameter:: myname_=myname//"::show_"
_ENTRY_(myname_)

  if(present(rc)) rc=0

  call show3dv_(spert, u_GCMpert, who, which, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'show3dv_("'//trim( u_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  call show3dv_(spert, v_GCMpert, who, which, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'show3dv_("'//trim( v_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  call show3dv_(spert,tv_GCMpert, who, which, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'show3dv_("'//trim(tv_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  call show3dv_(spert,dp_GCMpert, who, which, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'show3dv_("'//trim(dp_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  call show3dv_(spert,qv_GCMpert, who, which, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'show3dv_("'//trim(qv_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  call show3dv_(spert,ql_GCMpert, who, which, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'show3dv_("'//trim(ql_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  call show3dv_(spert,qi_GCMpert, who, which, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'show3dv_("'//trim(qi_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  call show3dv_(spert,o3_GCMpert, who, which, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'show3dv_("'//trim(o3_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

_EXIT_(myname_)
end subroutine show_
subroutine show3dv_(spert, snm, who, which, rc)
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_SUCCESS
  use MAPL_Mod, only: ESMFL_StateGetPointerToData
  use mpeu_util, only: tell,perr,die
  use kinds    , only: i_kind

  implicit none
  type(ESMF_State),target,intent(inout):: spert
  character(len=*),intent(in):: snm
  character(len=*),intent(in):: who,which
  integer(i_kind) ,intent(out):: rc

  integer(i_kind):: n,ier
  real(r_kind):: a,d
  real(r_ESkind),pointer,dimension(:,:,:):: sptr3
  character(len=*),parameter :: myname_=myname//'::show3dv_'
#ifdef _MYNAME_
#undef _MYNAME_
#endif
#define _MYNAME_ myname_//'('//trim(which)//'%"'//trim(snm)//'")'
_ENTRY_(_MYNAME_)

  sptr3 => null()
  call ESMFL_StateGetPointerToData(spert,sptr3,snm,rc=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'ESMFL_StateGetPointerToData('//trim(which)//'%"'//trim(snm)//'"), rc =',rc)
	  _EXIT_(_MYNAME_)
	  return
	endif

  	if(.not.associated(sptr3)) then
	  call perr(myname_,'.not.associated(ptr => '//trim(which)//'%"'//trim(snm)//'")')
	  rc=GENERIC_ERROR
	  _EXIT_(_MYNAME_)
	  return
	endif

  n=size(sptr3)
  a=sum(sptr3)/max(1,n)
  d=sum((sptr3-a)*(sptr3-a))/max(1,n-1)
  call tell(who,'local-size('//trim(which)//'%"'//trim(snm)//'") =',n)
  call tell(who,'local-mean('//trim(which)//'%"'//trim(snm)//'") =',a)
  call tell(who,'local-stdv('//trim(which)//'%"'//trim(snm)//'") =',d)

  nullify(sptr3)
_EXIT_(_MYNAME_)
end subroutine show3dv_

subroutine setNeeded_(spert,rc)
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_StateGet
  use ESMF, only: ESMF_StateIntent_Flag
  use ESMF, only: ESMF_StateIntent_EXPORT
  use ESMF, only: ESMF_SUCCESS
  use ESMF, only: operator(/=)
  use mpeu_util, only: perr,die
  use kinds    , only: i_kind
  implicit none
  type(ESMF_State),intent(inout):: spert
  integer(i_kind),optional,intent(out):: rc

  type(ESMF_StateIntent_Flag):: stype
  integer(i_kind):: ier
  character(len=*),parameter:: myname_=myname//"::setNeeded_"
_ENTRY_(myname_)

  if(present(rc)) rc=0

  call ESMF_StateGet(spert, stateIntent=stype, rc=ier)
  		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_StateGet(), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  _EXIT_(myname_)
		  return
		endif

		if(stype /= ESMF_STATEINTENT_EXPORT) then
		  call perr(myname_,'unexpected state-type, .not. ESMF_STATE_EXPORT')
		  if(.not.present(rc)) call die(myname_)
		  rc=GENERIC_ERROR
		  _EXIT_(myname_)
		  return
		endif

  call mySetNeeded_(spert, u_GCMpert, rc=ier)
	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'mySetNeeded_("'//trim( u_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  call mySetNeeded_(spert, v_GCMpert, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'mySetNeeded_("'//trim( v_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  call mySetNeeded_(spert,tv_GCMpert, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'mySetNeeded_("'//trim(tv_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  call mySetNeeded_(spert,dp_GCMpert, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'mySetNeeded_("'//trim(dp_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  call mySetNeeded_(spert,qv_GCMpert, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'mySetNeeded_("'//trim(qv_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  call mySetNeeded_(spert,ql_GCMpert, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'mySetNeeded_("'//trim(ql_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  call mySetNeeded_(spert,qi_GCMpert, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'mySetNeeded_("'//trim(qi_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  call mySetNeeded_(spert,o3_GCMpert, rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'mySetNeeded_("'//trim(o3_GCMpert)//'", rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

_EXIT_(myname_)
end subroutine setNeeded_

subroutine mySetNeeded_(s,snm,rc)
  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_SUCCESS
  use MAPL_Mod, only: ESMFL_StateGetPointerToData
  use kinds, only: i_kind,r_kind
  use mpeu_util, only: tell,perr,die
  implicit none
  type(ESMF_State), target,intent(inout):: s
  character(len=*),        intent(in   ):: snm	! variable name in s
  integer(i_kind) ,        intent(  out):: rc
  
  character(len=*),parameter:: myname_=myname//"::mySetNeeded_"
  real(r_ESkind),pointer,dimension(:,:,:):: sptr3
#ifdef _MYNAME_
#undef _MYNAME_
#endif
#define _MYNAME_ myname_//'%("'//trim(snm)//'")'
_ENTRY_(_MYNAME_)
  sptr3 => null()
  rc=ESMF_SUCCESS
  call ESMFL_StateGetPointerToData(s,sptr3,snm,alloc=.true.,rc=rc)
  	if(rc/=ESMF_SUCCESS) then
	  call perr(myname_,'ESMFL_StateGetPointerToData("'//trim(snm)//'",alloc), rc =',rc)
	  _EXIT_(_MYNAME_)
	  return
	endif
  nullify(sptr3)
_EXIT_(_MYNAME_)
end subroutine mySetNeeded_

subroutine getIFx_(pert,x,which,rc,teston)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 subroutine getIFx_
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:	 2010-11-04
!
! abstract: conditionally get from "x", a state in gsi_bundle
!
! program history log:
!   2010-11-04  j guo   - added this document block
!
!   input argument list: see Fortran 90 style document below
!
!   output argument list: see Fortran 90 style document below
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:
!
!$$$  end subprogram documentation block

! subroutine interface:

  use ESMF, only: ESMF_State
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlePutVar
  use kinds, only: i_kind,r_kind,r_quad
  use mpeu_util, only: tell,perr,die
  implicit none

  type(ESMF_State),        intent(inout):: pert	! output to a pre-configged state
  type(gsi_bundle),pointer              :: x	! input if associated().
  character(len=*),        intent(in   ):: which	! TLM or ADM
  integer(i_kind),optional,intent(  out):: rc
  logical        ,optional,intent(in   ):: teston

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname_=myname//'::getIFx_'
  integer(i_kind):: ier
  real(r_quad):: d
  type(gsi_bundle):: a,b
  logical:: teston_

_ENTRY_(myname_)
  if(present(rc)) rc=0
  teston_=.false.
  if(present(teston)) teston_=teston

	! Get import variables, to initialize them
  select case(associated(x))
  case(.true.)

    	! Adjoint of extending basic variables x%(tv,q,ps) to x%(prse, tsen)
    select case(which)
    case(ADMvect)
      if(teston_) call adtest_copy_(x,b,myname_)	! b = x -- deep_copy x b

      ! x = M''x -- x is the workspace
      call tv2tsen_ad_( ptr3d_to_(x,  tv_GSIpert), &
			ptr3d_to_(x,   q_GSIpert), &
			ptr3d_to_(x,tsen_GSIpert)  )
      call ps2prse_ad_( ptr2d_to_(x,  ps_GSIpert), &
			ptr3d_to_(x,  tv_GSIpert), &
			ptr3d_to_(x,prse_GSIpert)  )

      if(teston_) then
        call adtest_copy_(x,a,myname_)			! a = x -- deep_copy x a

	! a = Ma -- a is the workspace
	call ps2prse_tl_( ptr2d_to_(a,  ps_GSIpert), &
			  ptr3d_to_(a,  tv_GSIpert), &
			  ptr3d_to_(a,prse_GSIpert)  )
	call tv2tsen_tl_( ptr3d_to_(a,  tv_GSIpert), &
			  ptr3d_to_(a,   q_GSIpert), &
			  ptr3d_to_(a,tsen_GSIpert)  )

	!call adtest_show_(b,x,x,a,myname_)		! (x,x) .vs. (b,a)
	call adtest_dstr_(b,myname_)
	call adtest_dstr_(a,myname_)
      endif

!    case(TLMvect)
!      call gsi_bundlePutVar(x, tsen_GSIpert, 0._r_kind, ier)
!		if(ier/=0) then
!		  call perr(myname_,'gsi_bundlePutVar("'//tsen_GSIpert//'",0.), istatus =',ier)
!		  if(.not.present(rc)) call die(myname_)
!		  rc=ier
!		  _EXIT_(myname_)
!		  return
!		endif
!      call gsi_bundlePutVar(x, prse_GSIpert, 0._r_kind, ier)
!		if(ier/=0) then
!		  call perr(myname_,'gsi_bundlePutVar("'//prse_GSIpert//'",0.), istatus =',ier)
!		  if(.not.present(rc)) call die(myname_)
!		  rc=ier
!		  _EXIT_(myname_)
!		  return
!		endif
    endselect

    call gsi2pgcmIFx_(x,pert,which, stat=ier) !, teston=teston)
		if(ier/=0) then
		  call perr(myname_,'gsi2pgcmIFx_("'//which//'"), stat =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  _EXIT_(myname_)
		  return
		endif

  case default
    call pertState_set(pert,0._r_kind,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'pertState_set("'//which//'",0.), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(myname_)
	  return
	endif

  endselect

#ifdef DEBUG_VERBOSE
  call pertState_show(pert, myname_, which)
#endif

_EXIT_(myname_)
end subroutine getIFx_

subroutine adtest_copy_(xi,xo,where)
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundleCreate
  use gsi_bundlemod, only: assignment(=)
  use mpeu_util, only: perr,die
  implicit none
  type(gsi_bundle),intent(in ):: xi
  type(gsi_bundle),intent(out):: xo
  character(len=*),intent(in ):: where

  integer(i_kind):: ier
  character(len=*),parameter:: myname_='/adtest_copy_'
  call gsi_bundleCreate(xo,xi,'adtest_'//trim(xi%name), istatus=ier)
		if(ier/=0) then
		  call perr(where//myname_, &
		  	'gsi_bundleCreate(adtest_"'//trim(xi%name)//'"), istatus =',ier)
		  call  die(where//myname_)
		endif
  xo=xi
end subroutine adtest_copy_

subroutine adtest_dstr_(x,where)
  use kinds, only: i_kind
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundleDestroy
  use mpeu_util, only: perr,die
  implicit none
  type(gsi_bundle),intent(inout):: x
  character(len=*),intent(in   ):: where

  integer(i_kind):: ier
  character(len=*),parameter:: myname_='/adtest_dstr_'
  call gsi_bundleDestroy(x,istatus=ier)
  		if(ier/=0) then
		  call perr(where//myname_, &
		  	'gsi_bundleDestroy("'//trim(x%name)//'"), istatus =',ier)
		  call  die(where//myname_)
		endif
end subroutine adtest_dstr_

subroutine tv2tsen_tl_(p_tv,p_qv,p_tsen)
  use kinds, only: r_kind
  implicit none
  real(r_kind),pointer,dimension(:,:,:):: p_tv,p_qv,p_tsen
  call tv_to_tsen(p_tv,p_qv,p_tsen)
end subroutine tv2tsen_tl_
subroutine tv2tsen_ad_(p_tv,p_qv,p_tsen)
  use kinds, only: r_kind
  implicit none
  real(r_kind),pointer,dimension(:,:,:):: p_tv,p_qv,p_tsen
  call tv_to_tsen_ad(p_tv,p_qv,p_tsen)
end subroutine tv2tsen_ad_
subroutine ps2prse_tl_(p_ps,p_tv,p_prse)
  use kinds, only: r_kind
  implicit none
  real(r_kind),pointer,dimension(:,:  ):: p_ps
  real(r_kind),pointer,dimension(:,:,:):: p_tv,p_prse
  call getprs_tl(p_ps,p_tv,p_prse)
end subroutine ps2prse_tl_
subroutine ps2prse_ad_(p_ps,p_tv,p_prse)
  use kinds, only: r_kind
  implicit none
  real(r_kind),pointer,dimension(:,:  ):: p_ps
  real(r_kind),pointer,dimension(:,:,:):: p_tv,p_prse
  call getprs_ad(p_ps,p_tv,p_prse)
end subroutine ps2prse_ad_
function ptr3d_to_(x,var) result(ptr)
  use kinds, only: r_kind
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundleGetPointer
  use mpeu_util, only: perr,die
  implicit none
  type(gsi_bundle),target,intent(in):: x
  character(len=*),intent(in):: var
  real(r_kind),pointer,dimension(:,:,:):: ptr
  integer(i_kind):: rc
  character(len=*),parameter:: myname_=myname//"::ptr3d_to_"
  call gsi_bundleGetPointer(x,var,ptr,istatus=rc)
  		if(rc/=0) then
		  call perr(myname_,'gsi_bundleGetPointer("'// &
		  	trim(var)//'"), istatus =',rc)
		  call die(myname_)
		endif
end function ptr3d_to_
function ptr2d_to_(x,var) result(ptr)
  use kinds, only: r_kind
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundleGetPointer
  use mpeu_util, only: perr,die
  implicit none
  type(gsi_bundle),target,intent(in):: x
  character(len=*),intent(in):: var
  real(r_kind),pointer,dimension(:,:):: ptr
  integer(i_kind):: rc
  character(len=*),parameter:: myname_=myname//"::ptr2d_to_"
  call gsi_bundleGetPointer(x,var,ptr,istatus=rc)
  		if(rc/=0) then
		  call perr(myname_,'gsi_bundleGetPointer("'// &
		  	trim(var)//'"), istatus =',rc)
		  call die(myname_)
		endif
end function ptr2d_to_

subroutine putIFx_(pert,x,which,rc,teston)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 subroutine putIFx_
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:	 2010-11-04
!
! abstract: put into gsi_bundle variables from a pertState
!
! program history log:
!   2010-11-04  j guo   - added this document block
!
!   input argument list: see Fortran 90 style document below
!
!   output argument list: see Fortran 90 style document below
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:
!
!$$$  end subprogram documentation block

! subroutine interface:

  use ESMF, only: ESMF_State
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundlePutVar
  use kinds, only: i_kind,r_kind,r_quad
  use mpeu_util, only: tell,perr,die
  implicit none

  type(ESMF_State),        intent(inout):: pert	! input of a pre-configured state
  type(gsi_bundle),pointer              :: x	! output if associated().
  character(len=*),        intent(in   ):: which	! TLM or ADM
  integer(i_kind),optional,intent(  out):: rc
  logical        ,optional,intent(in   ):: teston

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname_=myname//'::putIFx_'
  integer(i_kind):: ier
  type(gsi_bundle):: a,b
  real(r_quad):: d
  logical:: teston_

_ENTRY_(myname_)
  if(present(rc)) rc=0
  teston_=.false.
  if(present(teston)) teston_=teston

#ifdef DEBUG_VERBOSE
  call pertState_show(pert, myname_, which)
#endif

	! Get import variables, to initialize them
  select case(associated(x))
  case(.true.)

    call pgcm2gsiIFx_(pert,x,which, stat=ier) !, teston=teston)
		if(ier/=0) then
		  call perr(myname_,'pgcm2gsiIFx_("'//which//'"), stat =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  _EXIT_(myname_)
		  return
		endif

    call gsi_bundlePutVar(x, tsen_GSIpert, 0._r_kind, ier)
		if(ier/=0) then
		  call perr(myname_,'gsi_bundlePutVar("'//tsen_GSIpert//'",0.), istatus =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  _EXIT_(myname_)
		  return
		endif
    call gsi_bundlePutVar(x, prse_GSIpert, 0._r_kind, ier)
		if(ier/=0) then
		  call perr(myname_,'gsi_bundlePutVar("'//prse_GSIpert//'",0.), istatus =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  _EXIT_(myname_)
		  return
		endif

    	! Extend variables x%(tv,q,ps) to x%(prse,tsen)
    select case(which)
    case(TLMvect)
      if(teston_) call adtest_copy_(x,b,myname_)	! b=x -- deep_copy x b

      ! x = Mx -- x is the workspace
      call ps2prse_tl_( ptr2d_to_(x,  ps_GSIpert), &
			ptr3d_to_(x,  tv_GSIpert), &
			ptr3d_to_(x,prse_GSIpert)  )
      call tv2tsen_tl_( ptr3d_to_(x,  tv_GSIpert), &
			ptr3d_to_(x,   q_GSIpert), &
			ptr3d_to_(x,tsen_GSIpert)  )

      if(teston_) then
        call adtest_copy_(x,a,myname_)		! a=x -- deep_copy x a

	! a = M''a -- a is the workspace
	call tv2tsen_ad_( ptr3d_to_(a,  tv_GSIpert), &
			  ptr3d_to_(a,   q_GSIpert), &
			  ptr3d_to_(a,tsen_GSIpert)  )
	call ps2prse_ad_( ptr2d_to_(a,  ps_GSIpert), &
			  ptr3d_to_(a,  tv_GSIpert), &
			  ptr3d_to_(a,prse_GSIpert)  )

	!call adtest_show_(b,x,x,a,myname_)		! (b,a) vs (x,x)
	call adtest_dstr_(b,myname_)
	call adtest_dstr_(a,myname_)
      endif

    endselect
  endselect

_EXIT_(myname_)
end subroutine putIFx_

subroutine set3d0_(s,snm,v,rc)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 subroutine set3d0_
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:	 2010-08-10
!
! abstract: Let s%("snm") = (scalor)v
!
! program history log:
!   2010-08-10  j guo   - added this document block
!
!   input argument list: see Fortran 90 style document below
!
!   output argument list: see Fortran 90 style document below
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:
!
!$$$  end subprogram documentation block

! subroutine interface:

  use ESMF, only: ESMF_State
  use MAPL_Mod, only: ESMFL_StateGetPointerToData
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundleGetPointer
  use kinds, only: i_kind,r_kind
  use mpeu_util, only: tell,perr,die
  implicit none
  type(ESMF_State), target,intent(inout):: s	! s%v = x%v
  character(len=*)        ,intent(in   ):: snm	! variable name in s
  real   (r_kind)         ,intent(in   ):: v	! a real value
  integer(i_kind),optional,intent(  out):: rc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname_=myname//"::set3d0_"
  integer(i_kind):: ier,n
  real(r_ESkind),pointer,dimension(:,:,:):: sptr3
  real(r_kind):: a,d

#ifdef _MYNAME_
#undef _MYNAME_
#endif
#define _MYNAME_ myname_//'%("'//trim(snm)//'")'
_ENTRY_(_MYNAME_)
  if(present(rc)) rc=0

! spert("x") = v
  sptr3 => null()
  call ESMFL_StateGetPointerToData(s,sptr3,snm,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'ESMFL_StateGetPointerToData("'//trim(snm)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(_MYNAME_)
	  return
	endif

  	if(.not.associated(sptr3)) then
	  call perr(myname_,'.not.associated(ptr => "'//trim(snm)//'")')
	  rc=GENERIC_ERROR
	  _EXIT_(_MYNAME_)
	  return
	endif

  sptr3(:,:,:)=v
  nullify(sptr3)

_EXIT_(_MYNAME_)
end subroutine set3d0_

subroutine set3dv_(s,snm,r,rc)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 subroutine set3dv_
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:	 2010-08-10
!
! abstract: Let s%("snm") = r%("snm")
!
! program history log:
!   2010-08-10  j guo   - added this document block
!
!   input argument list: see Fortran 90 style document below
!
!   output argument list: see Fortran 90 style document below
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:
!
!$$$  end subprogram documentation block

! subroutine interface:

  use ESMF, only: ESMF_State
  use ESMF, only: ESMF_SUCCESS
  use MAPL_Mod, only: ESMFL_StateGetPointerToData
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundleGetPointer
  use kinds, only: i_kind,r_kind
  use mpeu_util, only: tell,perr,die
  implicit none
  type(ESMF_State),target,intent(inout):: s	! s%v = r%v
  character(len=*),       intent(in   ):: snm	! variable name in s
  type(ESMF_State),target,intent(in   ):: r	!
  integer(kind=i_kind),optional,intent(  out):: rc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname_=myname//"::set3dv_"
  integer(i_kind):: ier,n
  type(ESMF_State),pointer:: r_
  real(r_ESkind),pointer,dimension(:,:,:):: rptr3
  real(r_ESkind),pointer,dimension(:,:,:):: sptr3

#ifdef _MYNAME_
#undef _MYNAME_
#endif
#define _MYNAME_ myname_//'%("'//trim(snm)//'")'
_ENTRY_(_MYNAME_)
  if(present(rc)) rc=0

  sptr3 => null()
  call ESMFL_StateGetPointerToData(s,sptr3,snm,rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'ESMFL_StateGetPointerToData(s%"'//trim(snm)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(_MYNAME_)
	  return
	endif


  	if(.not.associated(sptr3)) then
	  call perr(myname_,'.not.associated(ptr_exp => "'//trim(snm)//'")')
	  rc=GENERIC_ERROR
	  _EXIT_(_MYNAME_)
	  return
	endif

  r_ => r	! The first dummy argument (r) of ESMFL_StateGetPointerToData()
  		! is defined as INTENT(INOUT), which does not make the sense
		! to me.  So this (r_) is just trying to fool the compiler.
  rptr3 => null()
  call ESMFL_StateGetPointerToData(r_,rptr3,snm,rc=ier)
  	if(ier/=ESMF_SUCCESS) then
	  call perr(myname_,'ESMFL_StateGetPointerToData(r%"'//trim(snm)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(_MYNAME_)
	  return
	endif
  nullify(r_)

  	if(.not.associated(rptr3)) then
	  call perr(myname_,'.not.associated(ptr_imp => "'//trim(snm)//'")')
	  rc=GENERIC_ERROR
	  _EXIT_(_MYNAME_)
	  return
	endif

  sptr3(:,:,:)=rptr3(:,:,:)

  nullify(sptr3)
  nullify(rptr3)
_EXIT_(_MYNAME_)
end subroutine set3dv_

subroutine set2d0_(s,snm,v,rc)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 subroutine set2d0_
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:	 2010-08-10
!
! abstract: Let s%("snm") = v
!
! program history log:
!   2010-08-10  j guo   - added this document block
!
!   input argument list: see Fortran 90 style document below
!
!   output argument list: see Fortran 90 style document below
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:
!
!$$$  end subprogram documentation block

! subroutine interface:

  use ESMF, only: ESMF_State
  use MAPL_Mod, only: ESMFL_StateGetPointerToData
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundleGetPointer
  use kinds, only: i_kind,r_kind
  use mpeu_util, only: tell,perr,die
  implicit none
  type(ESMF_State) ,target,intent(inout):: s	! s%v = x%v
  character(len=*) ,       intent(in   ):: snm	! variable name in s
  real(kind=r_kind),       intent(in   ):: v	! a real value
  integer(kind=i_kind),optional,intent(out):: rc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname_=myname//"::set2d0_"
  integer(i_kind):: ier
  real(r_ESkind),pointer,dimension(:,:):: sptr2

#ifdef _MYNAME_
#undef _MYNAME_
#endif
#define _MYNAME_ myname_//'%("'//trim(snm)//'")'
_ENTRY_(_MYNAME_)
  if(present(rc)) rc=0

! s%("snm")=v
  sptr2 => null()
  call ESMFL_StateGetPointerToData(s,sptr2,snm,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'ESMFL_StateGetPointerToData("'//trim(snm)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(_MYNAME_)
	  return
	endif

  	if(.not.associated(sptr2)) then
	  call perr(myname_,'.not.associated(ptr => "'//trim(snm)//'")')
	  rc=GENERIC_ERROR
	  _EXIT_(_MYNAME_)
	  return
	endif

  sptr2(:,:)=v

  nullify(sptr2)
_EXIT_(_MYNAME_)
end subroutine set2d0_

subroutine set2dv_(s,snm,r,rc)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 subroutine set2dv_
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:	 2010-08-10
!
! abstract: Let s%("snm") = r%("rnm")
!
! program history log:
!   2010-08-10  j guo   - added this document block
!
!   input argument list: see Fortran 90 style document below
!
!   output argument list: see Fortran 90 style document below
!
! attributes:
!   language: Fortran 90 and/or above
!   machine:
!
!$$$  end subprogram documentation block

! subroutine interface:

  use ESMF, only: ESMF_State
  use MAPL_Mod, only: ESMFL_StateGetPointerToData
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundleGetPointer
  use kinds, only: i_kind,r_kind
  use mpeu_util, only: tell,perr,die
  implicit none
  type(ESMF_State),target,intent(inout):: s	! s%v = r%v
  character(len=*),       intent(in   ):: snm	! variable name in s and r
  type(ESMF_State),target,intent(in   ):: r	!
  integer(kind=i_kind),optional,intent(  out):: rc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname_=myname//"::set2dv_"
  integer(i_kind):: ier
  type(ESMF_State),pointer:: r_
  real(r_ESkind),pointer,dimension(:,:):: rptr2
  real(r_ESkind),pointer,dimension(:,:):: sptr2

#ifdef _MYNAME_
#undef _MYNAME_
#endif
#define _MYNAME_ myname_//'%("'//trim(snm)//'")'
_ENTRY_(_MYNAME_)
  if(present(rc)) rc=0

  sptr2 => null()
  call ESMFL_StateGetPointerToData(s,sptr2,snm,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'ESMFL_StateGetPointerToData("'//trim(snm)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(_MYNAME_)
	  return
	endif

  	if(.not.associated(sptr2)) then
	  call perr(myname_,'.not.associated(ptr_exp => "'//trim(snm)//'")')
	  rc=GENERIC_ERROR
	  _EXIT_(_MYNAME_)
	  return
	endif

  r_ => r	! The first dummy argument (r) of ESMFL_StateGetPointerToData()
  		! is defined as INTENT(INOUT), which does not make the sense
		! to me.  So this (r_) is just trying to fool the compiler.
  rptr2 => null()
  call ESMFL_StateGetPointerToData(r_,rptr2,snm,rc=ier)
  	if(ier/=0) then
	  call perr(myname_,'ESMFL_StateGetPointerToData("'//trim(snm)//'"), rc =',ier)
	  if(.not.present(rc)) call die(myname_)
	  rc=ier
	  _EXIT_(_MYNAME_)
	  return
	endif
  nullify(r_)

  	if(.not.associated(rptr2)) then
	  call perr(myname_,'.not.associated(ptr_imp => "'//trim(snm)//'")')
	  rc=GENERIC_ERROR
	  _EXIT_(_MYNAME_)
	  return
	endif

  sptr2(:,:)=rptr2(:,:)

  nullify(sptr2)
  nullify(rptr2)
_EXIT_(_MYNAME_)
end subroutine set2dv_

!------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1    !
!------------------------------------------------------------------------
!BOP
!
! !ROUTINE: gsi2pgcmIFx_:  Convert GSI increments to GCM perturbations
!
! !INTERFACE:

      subroutine gsi2pgcmIFx_ ( xx, spert, which, stat )

! !USES:

      use ESMF, only: ESMF_State
      use ESMF, only: ESMF_SUCCESS
      use MAPL_Mod, only: ESMFL_StateGetPointerToData

      use gsi_bundlemod, only: gsi_bundle
      use gsi_bundlemod, only: gsi_bundlegetpointer
      use gsi_bundlemod, only: assignment(=)

      use kinds, only: i_kind
      use kinds, only: r_kind
      implicit none

! !INPUT PARAMETERS:

      type(gsi_bundle), intent(in) :: xx    ! GSI increment vector
      character(len=*), intent(in) :: which ! adm or tlm

! !OUTPUT PARAMETERS:

      type(ESMF_State), intent(inout) :: spert ! GCM perturbation vector
      integer(i_kind) , intent(out)   :: stat  ! return error code

! !DESCRIPTION: Converts GSI increments in to ADM/TLM perturbations.
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  21Sep2007  Todling   Handles for O3 and CW.
!  19Nov2008  Todling   Update to use gsi-3d pressure instead of ps.
!  15May2010  Todling   Update to use GSI_Bundle
!  21Feb2011  Todling   Adapt to work with MAPL-SimpleBundle
!  30Mar2011  Guo       Reworked to apply to an ESMF_STATE:pert
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'::gsi2pgcmIFx_'

      	! temporary work space arrays in AGCM perturbation space.
      real(r_kind), allocatable, dimension(:,:,:) :: sub_tv,sub_u,sub_v,sub_q,sub_dp,sub_p3
      real(r_kind), allocatable, dimension(:,:  ) :: sub_ps
      real(r_kind), allocatable, dimension(:,:,:) :: sub_oz		  ! nc>=1
      real(r_kind), allocatable, dimension(:,:,:) :: sub_ql,sub_qi,sub_cw ! nc>=2

      	! indices to fields in input gsi_bundle:xx
      integer(i_kind)  i_u,i_v,i_t,i_q,i_p,i_oz,i_cw

      	! pointers to fields in output ESMF_State:spert
      real(r_ESkind),pointer,dimension(:,:,:):: p_u,p_v,p_tv,p_dp,p_qv,p_o3,p_ql,p_qi

      integer(i_kind)  i,j,k,ijk,ij
      integer(i_kind)  ierr,istatus,rc,status
      character(len=255) :: whatin
      integer(i_kind), parameter :: nc=2

      _ENTRY_(myname_)
      stat = ESMF_SUCCESS
      select case(which)
      case(ADMvect,TLMvect)
      case default
        call perr(myname_,'unknown mode, which =',which)
	stat=GENERIC_ERROR
        _EXIT_(myname_)
	return
      endselect

!     Initializes internal arrays
!     ---------------------------
      allocate (sub_u (lat2,lon2,nsig), sub_v (lat2,lon2,nsig),&
                sub_tv(lat2,lon2,nsig), sub_q (lat2,lon2,nsig),&
                sub_ps(lat2,lon2     ), sub_dp(lat2,lon2,nsig),&
		sub_p3(lat2,lon2,nsig+1), stat=ierr )
        if ( ierr/=0 ) then
	    call perr(myname_,'allocate(sub_u , ...), stat =',ierr)
	    stat=ierr
	    _EXIT_(myname_)
            return
        end if
      if ( nc>=1 ) then
          allocate (sub_oz(lat2,lon2,nsig), stat=ierr )
            if ( ierr/=0 ) then
	        call perr(myname_,'allocate(sub_oz), stat =',ierr)
	        stat=ierr
		_EXIT_(myname_)
                return
            end if
      endif
      if ( nc>=2 ) then
          allocate (sub_cw(lat2,lon2,nsig), &
	            sub_ql(lat2,lon2,nsig), &
                    sub_qi(lat2,lon2,nsig), stat=ierr )
            if ( ierr/=0 ) then
	        call perr(myname_,'allocate(sub_cw, sub_ql, sub_qi), stat =',ierr)
	        stat=ierr
		_EXIT_(myname_)
                return
            end if
      endif

!     Get poiners to GSI state-vector
!     -------------------------------
      ierr=0
      call gsi_bundlegetpointer(xx,   u_GSIpert, i_u   , istatus); if(istatus/=0) i_u   =-1
      call gsi_bundlegetpointer(xx,   v_GSIpert, i_v   , istatus); if(istatus/=0) i_v   =-1
      call gsi_bundlegetpointer(xx,  tv_GSIpert, i_t   , istatus); if(istatus/=0) i_t   =-1
      call gsi_bundlegetpointer(xx,   q_GSIpert, i_q   , istatus); if(istatus/=0) i_q   =-1
      call gsi_bundlegetpointer(xx,  ps_GSIpert, i_p   , istatus); if(istatus/=0) i_p   =-1

      if(min(i_u, i_v, i_t, i_q, i_p)<0) then
        call perr(myname_,'gsi_bundlegetpointer(xx), missing basic variable(s)')
	if(i_u   <0) call perr(myname_,'i_u   =',i_u   )
	if(i_v   <0) call perr(myname_,'i_v   =',i_v   )
	if(i_t   <0) call perr(myname_,'i_t   =',i_t   )
	if(i_q   <0) call perr(myname_,'i_q   =',i_q   )
	if(i_p   <0) call perr(myname_,'i_p   =',i_p   )
	stat=GENERIC_ERROR
	_EXIT_(myname_)
	return
      endif

      call gsi_bundlegetpointer(xx,  oz_GSIpert, i_oz  , istatus); if(istatus/=0) i_oz  =-1
      call gsi_bundlegetpointer(xx,  cw_GSIpert, i_cw  , istatus); if(istatus/=0) i_cw  =-1

      sub_u (:,:,:) = xx%r3(i_u )%q(:,:,:)
      sub_v (:,:,:) = xx%r3(i_v )%q(:,:,:)
      sub_tv(:,:,:) = xx%r3(i_t )%q(:,:,:)
      sub_q (:,:,:) = xx%r3(i_q )%q(:,:,:)
      sub_ps(:,:  ) = xx%r2(i_p )%q(:,:  )
      sub_dp(:,:,:) = 0.
      sub_p3(:,:,:) = 0.

      if(nc>=1) then
        sub_oz(:,:,:) = 0.
        if(i_oz>0) sub_oz(:,:,:) = xx%r3(i_oz)%q(:,:,:)
      endif

      if(nc>=2) then
        sub_cw(:,:,:) = 0.
        if(i_cw>0) sub_cw(:,:,:) = xx%r3(i_cw)%q(:,:,:)
        sub_ql(:,:,:) = 0.
        sub_qi(:,:,:) = 0.
      endif

!     convert perturbations from GSIpert space to GCMpert space
!     ---------------------------
      select case(which)
      case(ADMvect)
          do k=1,size(sub_dp,3)
            sub_dp(:,:,k) = sub_dp(:,:,k) + kPa_per_Pa * sub_ps(:,:)
          enddo
          sub_ps(:,:)=0.

	  	! adjoint-convert ozone unit in PPMV to GG
	  if(nc>=1) sub_oz(:,:,:) = sub_oz(:,:,:) * GG_per_PPMV

		! The forward operation is
		!    cw   = (1 1)(QL) = QL+QI
		!		 (QI)
		!
		! Its adjoint operation is
		!   (QL*) = (1)(cw*)  = (cw*)
		!   (QI*) = (1)         (cw*)

	  if(nc>=2) then
            sub_ql(:,:,:) = sub_cw(:,:,:)
            sub_qi(:,:,:) = sub_cw(:,:,:)
	    sub_cw(:,:,:) = 0.
	  endif

      case(TLMvect)
	  call getprs_tl(sub_ps,sub_tv,sub_p3)
	  sub_p3(:,:,:) = Pa_per_kPa * sub_p3(:,:,:)
	  do k=1,size(sub_dp,3)
            sub_dp(:,:,k) = sub_p3(:,:,k) - sub_p3(:,:,k+1)
	  enddo

	  	! convert ozone unit in GG to PPMV
	  if(nc>=1) sub_oz(:,:,:) = sub_oz(:,:,:) / GG_per_PPMV

		! The forward operation is
		!   (QL) = (1)(cw)
		!   (QI) = (0)
		!
		! Its adjoint operation is
		!    cw* = (1 0)(QL*) = QL*
		!		(QI*)
	  if(nc>=2) then
            sub_ql(:,:,:) = sub_cw(:,:,:)
	    sub_qi(:,:,:) = 0.
	  endif
      endselect

!     Gather from GSI subdomains/Scatter to GCM
!     -----------------------------------------
	  ierr=0
	  p_u => null()
	  call ESMFL_StateGetPointerToData(spert, p_u ,  u_GCMpert, rc=istatus)
	  	if(istatus/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMFL_StateGetPointerToData("'//trim(u_GCMpert)//'"), rc =',istatus)
		  ierr=GENERIC_ERROR
		endif
	  p_v => null()
	  call ESMFL_StateGetPointerToData(spert, p_v ,  v_GCMpert, rc=istatus)
	  	if(istatus/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMFL_StateGetPointerToData("'//trim(v_GCMpert)//'"), rc =',istatus)
		  ierr=GENERIC_ERROR
		endif
	  p_tv => null()
	  call ESMFL_StateGetPointerToData(spert, p_tv, tv_GCMpert, rc=istatus)
	  	if(istatus/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMFL_StateGetPointerToData("'//trim(tv_GCMpert)//'"), rc =',istatus)
		  ierr=GENERIC_ERROR
		endif
	  p_dp => null()
	  call ESMFL_StateGetPointerToData(spert, p_dp, dp_GCMpert, rc=istatus)
	  	if(istatus/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMFL_StateGetPointerToData("'//trim(dp_GCMpert)//'"), rc =',istatus)
		  ierr=GENERIC_ERROR
		endif
	  p_qv => null()
	  call ESMFL_StateGetPointerToData(spert, p_qv, qv_GCMpert, rc=istatus)
	  	if(istatus/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMFL_StateGetPointerToData("'//trim(qv_GCMpert)//'"), rc =',istatus)
		  ierr=GENERIC_ERROR
		endif

	  if(ierr/=0) then
	    stat=GENERIC_ERROR
	    _EXIT_(myname_)
	    return
	  endif

          call gsi2pert_ ( sub_u , p_u , ierr )
          call gsi2pert_ ( sub_v , p_v , ierr )
          call gsi2pert_ ( sub_tv, p_tv, ierr )
          call gsi2pert_ ( sub_dp, p_dp, ierr )
          call gsi2pert_ ( sub_q , p_qv, ierr )

	  p_o3 => null()
	  call ESMFL_StateGetPointerToData(spert, p_o3, o3_GCMpert, rc=istatus)
	  	if(istatus/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMFL_StateGetPointerToData("'//trim(o3_GCMpert)//'"), rc =',istatus)
		  ierr=GENERIC_ERROR
		endif

	  if(ierr/=0) then
	    stat=GENERIC_ERROR
	    _EXIT_(myname_)
	    return
	  endif

	  p_o3(:,:,:) = 0.
	  if(nc>=1) then
            call gsi2pert_ ( sub_oz, p_o3, ierr )
	    		if(ierr/=0) then
	  		  call perr(myname_,'unfinished O3 convertion, ierr =', ierr)
			  stat = ierr
			  _EXIT_(myname_)
			  return
			endif
	  endif

	  p_ql => null()
	  call ESMFL_StateGetPointerToData(spert, p_ql, ql_GCMpert, rc=istatus)
	  	if(istatus/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMFL_StateGetPointerToData("'//trim(ql_GCMpert)//'"), rc =',istatus)
		  ierr=GENERIC_ERROR
		endif
	  p_qi => null()
	  call ESMFL_StateGetPointerToData(spert, p_qi, qi_GCMpert, rc=istatus)
	  	if(istatus/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMFL_StateGetPointerToData("'//trim(qi_GCMpert)//'"), rc =',istatus)
		  ierr=GENERIC_ERROR
		endif

	  if(ierr/=0) then
	    stat=GENERIC_ERROR
	    _EXIT_(myname_)
	    return
	  endif

	  p_ql(:,:,:) = 0.
	  p_qi(:,:,:) = 0.
          if (nc>=2) then
            call gsi2pert_ ( sub_ql, p_ql, ierr )
	    		if(ierr/=0) then
	  		  call perr(myname_,'unfinished QL convertion, ierr =', ierr)
			  stat = ierr
			  _EXIT_(myname_)
			  return
			endif
            call gsi2pert_ ( sub_qi, p_qi, ierr )
	    		if(ierr/=0) then
	  		  call perr(myname_,'unfinished QI convertion, ierr =', ierr)
			  stat = ierr
			  _EXIT_(myname_)
			  return
			endif
	  endif

      deallocate (sub_u, sub_v, sub_tv, sub_q, sub_ps, sub_dp, sub_p3, stat=ierr )
		if ( ierr/=0 ) then
		  call perr(myname_,'deallocate(sub_u, ...), stat =', ierr)
		  stat=ierr
		  _EXIT_(myname_)
		  return
		endif
      if ( nc>=1 ) then
          deallocate (sub_oz, stat=ierr )
		if ( ierr/=0 ) then
		  call perr(myname_,'deallocate(sub_oz), stat =', ierr)
		  stat=ierr
		  _EXIT_(myname_)
		  return
		endif
      endif
      if ( nc>=2 ) then
          deallocate (sub_cw, sub_ql, sub_qi, stat=ierr )
		if ( ierr/=0 ) then
		  call perr(myname_,'deallocate(sub_cw, sub_ql, sub_qi), stat =', ierr)
		  stat=ierr
		  _EXIT_(myname_)
		  return
		end if
      endif

      _EXIT_(myname_)
      CONTAINS

      subroutine gsi2pert_ ( sub, fld, stat_ )

      implicit none
      real(r_kind),   intent(in)  :: sub(:,:,:)
      real(r_ESkind), intent(inout) :: fld(:,:,:)
      integer(i_kind),intent(out) :: stat_

      character(len=*), parameter :: myname_ = myname//'*gsi2pert_'

      real(r_kind), allocatable :: fldsm(:,:)
      real(r_kind), allocatable :: work4d(:,:,:,:)   ! auxliar 4d array
      real(r_kind), allocatable :: work3d(:,:,:)     ! auxliar 3d array
      real(r_kind), allocatable :: work(:)

      integer(i_kind) imr,jnp,nl
      integer(i_kind) i,j,k,mm1,ierr

      mm1 = mype+1 
      stat_ = 0

      allocate ( work3d(nlon,nlat,nsig), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work3d)'
            return
        end if
      allocate ( work(max(iglobal,itotsub)), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work)'
            return
        end if
      allocate ( fldsm(lat1*lon1,nsig), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(fldsm)'
            return
        end if

!     Strip off boundary points from subdomains
!     -----------------------------------------
      call strip(sub,fldsm,nsig)

!     Gather GSI perturbations to root processor
!     ------------------------------------------
      do k=1,nsig
         call mpi_gatherv(fldsm(1,k),ijn(mm1),mpi_rtype,&
              work,ijn,displs_g,mpi_rtype,&
              ROOT,mpi_comm_world,ierr)
         if (mype==ROOT) then
            call reorder12(work,work3d(:,:,k))
         endif
      end do

      deallocate ( fldsm, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': delloc(fldsm)'
            return
        end if
      deallocate ( work, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': delloc(work)'
            return
        end if

      call mygetdims_(imr,jnp,nl)
      allocate ( work4d(imr,jnp,nl,1), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work4d)'
            return
        end if

!     Flip horizontal and vertical
!     ----------------------------
      if ( mype==ROOT ) then
              call hflip3_ ( work3d, imr,jnp,nl )
              if (imr/=nlon .or. jnp/=nlat ) then
                  if (which==ADMvect) then
                      work4d = zero
                      call interpack_terpv_ad ( imr,jnp,nl,work4d(:,:,:,1),work4d(:,:,:,1), nlon,nlat,nsig,work3d, ierr )
                  else if (which==TLMvect) then
                      work4d = zero
                      call interpack_terpv    ( nlon,nlat,nsig,work3d, imr,jnp,nl,work4d(:,:,:,1), ierr )
                  else
                      call die ( myname_,': invalid option' )
                  endif
              else
                  work4d(:,:,:,1) = work3d(:,:,:)
              endif
           call SwapV_ ( work4d(:,:,:,1) )
      endif

!     Scatter perturbations to GCM decomposition
!     ------------------------------------------
      call myscatter_ ( which, imr, jnp, nl, work4d(:,:,:,1), fld )

!     Swap work memory
!     ----------------
      deallocate ( work4d, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(work4d)'
            return
        end if

      deallocate ( work3d, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': delloc(work3d)'
            return
        end if

      end subroutine gsi2pert_

      end subroutine gsi2pgcmIFx_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pgcm2gsiIFx_:  Convert gcm-adm/tlm vector to gsi vector
!
! !INTERFACE:

      subroutine pgcm2gsiIFx_ ( spert, xx, which, stat )

! !USES:

      use ESMF, only: ESMF_State
      use ESMF, only: ESMF_SUCCESS
      use MAPL_Mod, only: ESMFL_StateGetPointerToData

      use constants, only: one,zero
      use gsi_bundlemod,only: gsi_bundle     ! GSI state vector
      use gsi_bundlemod,only: gsi_bundlegetpointer
      use gsi_bundlemod,only: assignment(=)
      use state_vectors, only: dot_product
      use kinds, only: r_kind, r_quad
 
      implicit none

! !INPUT PARAMETERS:

      type(ESMF_State), intent(inout)  :: spert  ! GCM perturbation vector 
      character(len=*), intent(in   )  :: which  ! adm or tlm

! !OUTPUT PARAMETERS:

      type(gsi_bundle), intent(inout) :: xx     ! GSI increment
      integer(i_kind) , intent(  out) :: stat

! !DESCRIPTION: Convert GEOS-5 perturbation vector to GSI increment vector
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  21Sep2007  Todling   Handles for O3 and CW.
!  22Feb2008  Todling   Handle for forecast (error) gradient vector.
!  19Nov2008  Todling   Update to use gsi-3d pressure instead of ps.
!  15May2010  Todling   Update to use GSI_Bundle
!  20Feb2011  Todling   Adapt to work with MAPL-SimpleBundle
!  30Mar2011  Guo       Reworked to apply to an ESMF_STATE:pert
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'::pgcm2gsiIFx_'

      real(r_kind),  allocatable, dimension(:,:,:) :: sub_u,sub_v,sub_dp,sub_q,sub_tv,sub_p3
      real(r_kind),  allocatable, dimension(:,:  ) :: sub_ps
      real(r_kind),  allocatable, dimension(:,:,:) :: sub_oz			! nc>=1
      real(r_kind),  allocatable, dimension(:,:,:) :: sub_ql,sub_qi,sub_cw	! nc>=2

      integer(i_kind) i,j,k,ij,ijk,id,ng_d,ng_s
      integer(i_kind) i_u,i_v,i_t,i_q,i_oz,i_cw,i_p

      real(r_ESkind),pointer,dimension(:,:,:):: p_u,p_v,p_tv,p_dp,p_qv,p_o3,p_ql,p_qi
      integer(i_kind) ierr,istatus,status
      integer(i_kind), parameter :: nc=2
      real(r_quad):: d

      _ENTRY_(myname_)

      stat = 0
      stat = ESMF_SUCCESS
      select case(which)
      case(ADMvect,TLMvect)
      case default
        call perr(myname_,'unknown mode, which =',which)
	stat=GENERIC_ERROR
	_EXIT_(myname_)
	return
      endselect

      allocate (sub_u (lat2,lon2,nsig), sub_v (lat2,lon2,nsig),&
                sub_tv(lat2,lon2,nsig), sub_q (lat2,lon2,nsig),&
                sub_ps(lat2,lon2     ), sub_dp(lat2,lon2,nsig),&
                sub_p3(lat2,lon2,nsig+1), stat=ierr )
	if ( ierr/=0 ) then
	    call perr(myname_,'allocate(sub_u, ...), stat =',ierr)
	    stat=ierr
	    _EXIT_(myname_)
	    return
	end if
      if ( nc>=1 ) then
          allocate (sub_oz(lat2,lon2,nsig), stat=ierr )
            if ( ierr/=0 ) then
	        call perr(myname_,'allocate(sub_oz), stat =',ierr)
	        stat=ierr
	        _EXIT_(myname_)
                return
            end if
      endif
      if ( nc>=2 ) then
          allocate (sub_cw(lat2,lon2,nsig), &
                    sub_ql(lat2,lon2,nsig), &
                    sub_qi(lat2,lon2,nsig), stat=ierr )
            if ( ierr/=0 ) then
	        call perr(myname_,'allocate(sub_cw, sub_ql, sub_qi), stat =',ierr)
	        stat=ierr
	        _EXIT_(myname_)
                return
            end if
      endif

!     Gather from GCM/Scatter to GSI subdomains
!     -----------------------------------------
	  p_u => null()
	  call ESMFL_StateGetPointerToData(spert, p_u ,  u_GCMpert)
	  p_v => null()
	  call ESMFL_StateGetPointerToData(spert, p_v ,  v_GCMpert)
	  p_tv => null()
	  call ESMFL_StateGetPointerToData(spert, p_tv, tv_GCMpert)
	  p_dp => null()
	  call ESMFL_StateGetPointerToData(spert, p_dp, dp_GCMpert)
	  p_qv => null()
	  call ESMFL_StateGetPointerToData(spert, p_qv, qv_GCMpert)

	  if( .not.associated(p_u ) .or. .not.associated(p_v ) .or. &
	      .not.associated(p_tv) .or. .not.associated(p_dp) .or. &
	      .not.associated(p_qv) ) then
	    if(.not.associated(p_u )) call perr(myname_,'.not.associated(p_u )')
	    if(.not.associated(p_v )) call perr(myname_,'.not.associated(p_v )')
	    if(.not.associated(p_tv)) call perr(myname_,'.not.associated(p_tv)')
	    if(.not.associated(p_dp)) call perr(myname_,'.not.associated(p_dp)')
	    if(.not.associated(p_qv)) call perr(myname_,'.not.associated(p_qv)')
	    call die(myname_)
	  endif

          call pert2gsi_ ( p_u , sub_u , ierr )
	  	if(ierr/=0) call die(myname_,'pert2gsi(p_u ), ierr',ierr)
          call pert2gsi_ ( p_v , sub_v , ierr )
	  	if(ierr/=0) call die(myname_,'pert2gsi(p_v ), ierr',ierr)
          call pert2gsi_ ( p_tv, sub_tv, ierr )
	  	if(ierr/=0) call die(myname_,'pert2gsi(p_tv), ierr',ierr)
          call pert2gsi_ ( p_dp, sub_dp, ierr )
	  	if(ierr/=0) call die(myname_,'pert2gsi(p_dp), ierr',ierr)
          call pert2gsi_ ( p_qv, sub_q , ierr )
	  	if(ierr/=0) call die(myname_,'pert2gsi(p_qv), ierr',ierr)

      if (nc>=1) then
	  p_o3 => null()
	  call ESMFL_StateGetPointerToData(spert, p_o3, o3_GCMpert)
	  if( .not.associated(p_o3) ) then
	    call perr(myname_,'.not.associated(p_o3)')
	    call die(myname_)
	  endif
          call pert2gsi_ ( p_o3, sub_oz, ierr )
		if ( ierr/=0 ) then
		  stat = 99
		  if(mype==ROOT) print*, trim(myname_), ': unfinished O3 convertion ...'
		  _EXIT_(myname_)
		  return
		end if
      endif
      if (nc>=2) then
	  p_ql => null()
	  call ESMFL_StateGetPointerToData(spert, p_ql, ql_GCMpert)
	  p_qi => null()
	  call ESMFL_StateGetPointerToData(spert, p_qi, qi_GCMpert)

	  if( .not.associated(p_ql) .or. .not.associated(p_qi) ) then
	    if(.not.associated(p_ql)) call perr(myname_,'.not.associated(p_ql)')
	    if(.not.associated(p_qi)) call perr(myname_,'.not.associated(p_qi)')
	    call die(myname_)
	  endif

          call pert2gsi_ ( p_ql, sub_ql, ierr )
		if ( ierr/=0 ) then
		  stat = 99
		  if(mype==ROOT) print*, trim(myname_), ': unfinished QL convertion ...'
		  _EXIT_(myname_)
		  return
		end if
          call pert2gsi_ ( p_qi, sub_qi, ierr )
		if ( ierr/=0 ) then
		  stat = 99
		  if(mype==ROOT) print*, trim(myname_), ': unfinished QI convertion ...'
		  _EXIT_(myname_)
		  return
		end if
      endif

!     Get poiners to GSI state-vector
!     -------------------------------
      ierr=0
      call gsi_bundlegetpointer(xx,   u_GSIpert, i_u   , istatus); if(istatus/=0) i_u   =-1
      call gsi_bundlegetpointer(xx,   v_GSIpert, i_v   , istatus); if(istatus/=0) i_v   =-1
      call gsi_bundlegetpointer(xx,  tv_GSIpert, i_t   , istatus); if(istatus/=0) i_t   =-1
      call gsi_bundlegetpointer(xx,   q_GSIpert, i_q   , istatus); if(istatus/=0) i_q   =-1
      call gsi_bundlegetpointer(xx,  ps_GSIpert, i_p   , istatus); if(istatus/=0) i_p   =-1

      if(min(i_u, i_v, i_t, i_q, i_p)<0) then
        call perr(myname_,'gsi_bundlegetpointer(xx), missing basic variable(s)')
	if(i_u   <0) call perr(myname_,'i_u   =',i_u   )
	if(i_v   <0) call perr(myname_,'i_v   =',i_v   )
	if(i_t   <0) call perr(myname_,'i_t   =',i_t   )
	if(i_q   <0) call perr(myname_,'i_q   =',i_q   )
	if(i_p   <0) call perr(myname_,'i_p   =',i_p   )
	stat=GENERIC_ERROR
	_EXIT_(myname_)
	return
      endif

      call gsi_bundlegetpointer(xx,  oz_GSIpert, i_oz  , istatus); if(istatus/=0) i_oz  =-1
      call gsi_bundlegetpointer(xx,  cw_GSIpert, i_cw  , istatus); if(istatus/=0) i_cw  =-1

!     Calculate all other perturbation for GSI
!     ----------------------------------------
      select case(which)
      case(ADMvect)
          sub_ps(:,:) = 0.

	  	! dp -> prse -> (ps,tv)
          sub_p3(:,:,:) = 0.
          sub_p3(:,:,1:nsig  ) = sub_p3(:,:,1:nsig  ) + sub_dp(:,:,1:nsig)
          sub_p3(:,:,2:nsig+1) = sub_p3(:,:,2:nsig+1) - sub_dp(:,:,1:nsig)
          sub_p3(:,:,:) = Pa_per_kPa * sub_p3(:,:,:)
          call getprs_ad(sub_ps,sub_tv,sub_p3)
	  sub_dp(:,:,:) = 0.

	  	! adjoint-convert ozone unit in GG to PPMV
	  if(nc>=1) sub_oz(:,:,:) = sub_oz(:,:,:) / GG_per_PPMV

		! The forward operation is
		!   (QL) = (1)(cw)
		!   (QI) = (0)
		!
		! Its adjoint operation is
		!    cw* = (1 0)(QL*) = QL*
		!		(QI*)

	  if(nc>=2) then
	    sub_cw(:,:,:) = sub_ql (:,:,:)
	    sub_ql(:,:,:) = 0.
	    sub_qi(:,:,:) = 0.
	  endif

      case(TLMvect)
      	! dp -> ps
	  sub_ps(:,:)=0.
	  do k=1,size(sub_dp,3)
	    sub_ps(:,:) = sub_ps(:,:) + sub_dp(:,:,k)
	  enddo
	  sub_ps(:,:) = kPa_per_Pa * sub_ps(:,:)

	  	! convert ozone unit in PPMV to GG
	  if(nc>=1) sub_oz(:,:,:) = sub_oz(:,:,:) * GG_per_PPMV

		! The forward operation is
		!    cw   = (1 1)(QL) = QL+QI
		!		 (QI)
		!
		! Its adjoint operation is
		!   (QL*) = (1)(cw*)  = (cw*)
		!   (QI*) = (1)         (cw*)

	  if(nc>=2) sub_cw(:,:,:) = sub_ql(:,:,:) + sub_qi(:,:,:)
      endselect

!     Copy pertbations to xx in gsiPert space

      xx%r2(i_p)%q(:,:  ) = sub_ps(:,:  )
      xx%r3(i_t)%q(:,:,:) = sub_tv(:,:,:)
      xx%r3(i_q)%q(:,:,:) = sub_q (:,:,:)
      xx%r3(i_v)%q(:,:,:) = sub_v (:,:,:)
      xx%r3(i_u)%q(:,:,:) = sub_u (:,:,:)

      if(i_oz>0) then
        xx%r3(i_oz)%q(:,:,:) = 0.
        if(nc>=1) xx%r3(i_oz)%q(:,:,:) = sub_oz(:,:,:)
      endif

      if(i_cw>0) then
        xx%r3(i_cw)%q(:,:,:) = 0.
        if(nc>=2) xx%r3(i_cw)%q(:,:,:) = sub_cw(:,:,:)
      endif

!     The following will be left untouched
!     ------------------------------------
!     xx%sst(:) = xx%sst(:)

      deallocate (sub_u, sub_v, sub_tv, sub_q, sub_dp, sub_ps, stat=ierr )
        if ( ierr/=0 ) then
            stat = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_)'
	    _EXIT_(myname_)
            return
        end if
      if ( nc>=1 ) then
          deallocate (sub_oz, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_oz)'
		_EXIT_(myname_)
                return
            end if
       endif
      if ( nc>=2 ) then
          deallocate (sub_cw, sub_ql, sub_qi, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_cw, sub_ql, sub_qi)'
		_EXIT_(myname_)
                return
            end if
       endif

      _EXIT_(myname_)
      CONTAINS

      subroutine pert2gsi_ ( fld, sub, stat_ )

      implicit none
      real(r_ESkind), intent(in)  :: fld(:,:,:)
      real(r_kind),   intent(out) :: sub(:,:,:)
      integer(i_kind),intent(out) :: stat_

      character(len=*), parameter :: myname_ = myname//'*pert2gsi_'

      real(r_kind), allocatable :: work4d(:,:,:,:)   ! auxliar 4d array
      real(r_kind), allocatable :: work3d(:,:,:)     ! auxliar 3d array
      real(r_kind), allocatable :: work2d(:,:)       ! auxliar 2d array
      real(r_kind), allocatable :: work(:)
      integer(i_kind) i,j,k,mm1,ierr
      integer(i_kind) imr,jnp,nl

      mm1 = mype+1 
      stat_ = 0
      call mygetdims_(imr,jnp,nl)

      allocate ( work3d(nlon,nlat,nsig), stat=ierr )
	if ( ierr/=0 ) then
	  call perr(myname_,'allocate(work3d), stat =', ierr)
	  stat=ierr
	  return
	end if
      allocate ( work4d(imr,jnp,nl,1), stat=ierr )
	if ( ierr/=0 ) then
	  call perr(myname_,'allocate(work4d), stat =', ierr)
	  stat=ierr
	  return
	end if
                                                                                                                           
!     Gather GCM perturbations to root processor
!     ------------------------------------------
      call mygather_( which, imr, jnp, nl, fld, work4d(:,:,:,1) )

!     Flip horizontal and vertical
!     ----------------------------
      if ( mype==ROOT ) then
              call hflip3_ ( work4d, imr,jnp,nl )
              if (imr/=nlon .or. jnp/=nlat ) then
                  if (which==ADMvect) then
                      work3d = zero
                      call interpack_terpv_ad ( nlon,nlat,nsig,work3d,work3d, imr,jnp,nl,work4d(:,:,:,1), ierr )
                  else if (which==TLMvect) then
                      work3d = zero
                      call interpack_terpv    ( imr,jnp,nl,work4d(:,:,:,1),  nlon,nlat,nsig,work3d, ierr )
                  else
                      call die ( myname_,': invalid option' )
                  endif
              else
                  work3d(:,:,:) = work4d(:,:,:,1)
              endif
           call SwapV_ ( work3d )
      endif

!     Swap work memory
!     ----------------
      deallocate ( work4d, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(work4d)'
            return
        end if
      allocate ( work(itotsub), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work)'
            return
        end if
      allocate ( work2d(lat2,lon2), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work4d)'
            return
        end if

!     Scatter to GSI subdomains
!     -------------------------
      do k=1,nsig
         if (mype==ROOT) then
             call reorder21(work3d(:,:,k),work)
         endif
         call mpi_scatterv(work,ijn_s,displs_s,mpi_rtype,&
              work2d,ijn_s(mm1),mpi_rtype,root,mpi_comm_world,ierr)
         do j=1,lon2
            do i=1,lat2
               sub(i,j,k) = work2d(i,j)
            end do
         end do
      end do

!     Release work memory
!     -------------------
      deallocate ( work2d, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(work4d)'
            return
        end if
      deallocate ( work, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': delloc(work)'
            return
        end if
      deallocate ( work3d, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': delloc(work3d)'
            return
        end if

      end subroutine pert2gsi_

   end subroutine pgcm2gsiIFx_
!------ BELOW THIS POINT: Routines of general (internal-only) use --------

!-------------------------------------------------------------------------
   subroutine SwapV_(fld)
!-------------------------------------------------------------------------
   implicit none
   real(r_kind),intent(inout) ::  fld(:,:,:)
   real(r_kind),allocatable   :: work(:,:,:)
   integer(i_kind) im, jm, km
   im   = size(fld,1)
   jm   = size(fld,2)
   km   = size(fld,3)
   allocate (work(im,jm,km))
   work = fld
   fld(:,:,km:1:-1) = work(:,:,1:km:+1)
   deallocate (work)
   end subroutine SwapV_
!-------------------------------------------------------------------------
   subroutine SwapIJK_(aij,aji)
!-------------------------------------------------------------------------
! transpose IJK-ordered array to JIK-ordered array
   implicit none
   real(r_kind),dimension(:,:,:),intent(in)    :: aij
   real(r_kind),dimension(:,:,:),intent(inout) :: aji
   integer(i_kind) :: i,k,isz,jsz,ksz,kk
!
   isz=size(aij,1)
   jsz=size(aij,2)
   ksz=size(aij,3)
   kk=1
   do k=1,ksz
      do i=1,isz
         aji(1:jsz,i,kk)=aij(i,1:jsz,k)
      end do
      kk=kk+1
   end do
   call SwapV_(aji)
   end subroutine SwapIJK_

   subroutine hflip3_ ( q,im,jm,km )
      implicit none
      integer(i_kind)  im,jm,km,i,j,k
      real(r_kind), intent(inout) :: q(im,jm,km)
      real(r_kind), allocatable   :: dum(:)
      allocate ( dum(im) )
      do k=1,km
      do j=1,jm
      do i=1,im/2
         dum(i) = q(i+im/2,j,k)
         dum(i+im/2) = q(i,j,k)
      enddo
         q(:,j,k) = dum(:)
      enddo
      enddo
      deallocate ( dum )
   end subroutine hflip3_
  
   subroutine hflip2_ ( q,im,jm )
      implicit none
      integer(i_kind)  im,jm,i,j
      real(r_kind), intent(inout) :: q(im,jm)
      real(r_kind), allocatable   :: dum(:)
      allocate ( dum(im) )
      do j=1,jm
      do i=1,im/2
         dum(i) = q(i+im/2,j)
         dum(i+im/2) = q(i,j)
      enddo
         q(:,j) = dum(:)
      enddo
      deallocate ( dum )
   end subroutine hflip2_

   subroutine mygetdims_ ( im_world, jm_world, km_world )
   use GSI_GridCompMod, only: GSI_AgcmPertGrid
   use ESMF, only: ESMF_MAXGRIDDIM
   use MAPL_Mod, only: MAPL_GridGet
   use MAPL_Mod, only: MAPL_VRFY
   implicit none
   integer(i_kind), intent(out) :: im_world
   integer(i_kind), intent(out) :: jm_world
   integer(i_kind), intent(out) :: km_world
  
   character(len=*),parameter :: Iam=myname//'mygetdims_'
   integer  DIMS(ESMF_MAXGRIDDIM)
   integer  rc,status

   call MAPL_GridGet(GSI_AgcmPertGrid, globalCellCountPerDim=DIMS, RC=STATUS)
        VERIFY_(STATUS)

   im_world = DIMS(1)
   jm_world = DIMS(2)
   km_world = DIMS(3)
   end subroutine mygetdims_

   subroutine myscatter_ ( what, im_world, jm_world, km_world, vari, varo )
!  NOTE: there is not need for Halo here ... this is going onto a "model" field
   use MAPL_Mod, only: ArrayScatter
   use MAPL_Mod, only: MAPL_VRFY
   use GSI_GridCompMod, only: GSI_AgcmPertGrid
   implicit none
   character(len=*),intent(in)    :: what
   integer,         intent(in)    :: im_world,jm_world,km_world
   real(r_kind),    intent(in)    :: vari(im_world,jm_world,km_world)
   real(r_ESkind),  intent(inout) :: varo(:,:,:)

   character(len=*),parameter :: Iam=myname//'myscatter_'
   real(r_ESkind),allocatable :: work2d(:,:)
   integer  k, rc,status

    allocate(work2d(im_world,jm_world))
    if (what==ADMvect) then
        do k=1,km_world  ! this would be a good place to swapV
           work2d=vari(:,:,k)
           call ArrayScatter(varo(:,:,k), work2d, GSI_AgcmPertGrid, rc=status)   ! _RT: ideally this should be the AD of gather
                VERIFY_(STATUS)
        enddo
    else
        do k=1,km_world  ! this would be a good place to swapV
           work2d=vari(:,:,k)
           call ArrayScatter(varo(:,:,k), work2d, GSI_AgcmPertGrid, rc=status)   ! _RT: ideally this should be the AD of gather
                VERIFY_(STATUS)
        enddo
    endif
    deallocate(work2d)
    
   end subroutine myscatter_

   subroutine mygather_ ( what, im_world, jm_world, km_world, vari, varo )
!  NOTE: there is not need for Halo here ... this is coming from a "model" field
   use MAPL_Mod, only: ArrayGather
   use MAPL_Mod, only: MAPL_VRFY
   use GSI_GridCompMod, only: GSI_AgcmPertGrid
   implicit none
   character(len=*),        intent(in)    :: what
   integer,                 intent(in)    :: im_world,jm_world,km_world
   real(r_ESkind),          intent(in)    :: vari(:,:,:)
   real(r_kind),            intent(inout) :: varo(im_world,jm_world,km_world)

   character(len=*),parameter :: Iam=myname//'mygather_'
   real(r_ESkind),allocatable :: work2d(:,:)
   integer  k, rc,status

    allocate(work2d(im_world,jm_world))
    if (what==ADMvect) then
        do k=1,km_world  ! this would be a good place to swapV
           work2d=zero
           call ArrayGather(vari(:,:,k), work2d, GSI_AgcmPertGrid, rc=status)   ! _RT: ideally this should be the AD of gather
                VERIFY_(STATUS)
           varo(:,:,k) = work2d
        enddo
    else
        do k=1,km_world  ! this would be a good place to swapV
           work2d=zero
           call ArrayGather(vari(:,:,k), work2d, GSI_AgcmPertGrid, rc=status)
                VERIFY_(STATUS)
           varo(:,:,k) = work2d
        enddo
    endif
    deallocate(work2d)
    
   end subroutine mygather_
end module geos_pertState
