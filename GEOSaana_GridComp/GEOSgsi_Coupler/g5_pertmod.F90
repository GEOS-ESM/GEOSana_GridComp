!#define DEBUG_VERBOSE
!#define DEBUG_TRACE
!#define DEBUG_CLOCK
!#define DEBUG_STATES
!#define DEBUG_MPOUT

   module g5_pertmod
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: g5_pertmod --- a g5pert based GSI perturbation model component 
!
! !DESCRIPTION: This module provides a wrapper around g5pert/ component
!
! !USES:

   use MAPL     , only: MAPL_AM_I_ROOT
   use kinds    , only: i_kind,r_kind
   use constants, only: R3600

!  ADM/TLM entries
!  ---------------
#ifdef GEOS_PERT
      use precision
      use mod_comm,    only : mp_init
      use mod_comm,    only : mp_gather4d
      use mod_comm,    only : mp_scatter4d
      use stepon,      only : ng_d, ng_s
      use stepon,      only : stepon_set
      use stepon,      only : nymd, nhms
      use stepon,      only : pdt         ! time step of AD/TL models
      use prognostics, only : imr         ! no. of grid points in the zonal direction
      use prognostics, only : jnp         ! no. of grid points in the meridional direction
      use prognostics, only : nl          ! no. of levels
      use prognostics, only : nc          ! no. of tracers
      use prognostics, only : jfirst      ! pointer for lat decomposition
      use prognostics, only : jlast       ! pointer for lat decomposition
      use prognostics, only : prognostics_initial
      use prognostics, only : prognostics_final
      use prognostics, only : prognostics_dotp
      use prognostics, only : prognostics_dup
      use prognostics, only : prognostics_cnst
      use m_iostate,   only : getstate_init
      use m_iostate,   only : getstate
      use m_iostate,   only : getstate_clean
      use m_trjphys,   only : physdrv1_get_init
      use m_trjphys,   only : physdrv1_get_all
      use m_trjphys,   only : physdrv1_get_clean
      use m_trajmng,   only : getpert
      use m_trajmng,   only : putpert
      use m_trajmng,   only : inqpert_dims
      use prognostics, only : dyn_prog    ! GCM perturbation vector

      use m_interpack,    only : interpack_terpv
      use m_interpack_ad, only : interpack_terpv_ad

      use stepon_tl,   only : stepon_g4tog5_tl

      use m_die,       only : die

      use GSI_GridCompMod, only: PPMV2GpG

!  the following will be cleared once I update the interface to putpert
!  ....................................................................
      use stepon,      only : ak          ! GEOS-5 ADM/TLM pressure levels
      use stepon,      only : bk          ! GEOS-5 ADM/TLM pressure levels
      use stepon,      only : ts
      use stepon,      only : oro
      use stepon,      only : job
      use stepon,      only : nstep
      use stepon,      only : fvpsasdt

#endif /* GEOS_PERT */

!  GSI entries
!  -----------
      use kinds,       only : r_kind,r_single,i_kind
      use mpimod,      only : mype,mpi_rtype,mpi_comm_world
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
      use gsi_bundlemod,only: gsi_bundle     ! GSI state vector
      use gsi_bundlemod,only: gsi_bundlegetpointer
      use gsi_bundlemod,only: assignment(=)

   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: pertmod_setServices ! interface to ESMF SetServices
   	!!! call cplr::pertmod_setServices_(rc)

   public :: pertmod_initialize	! interface to ESMF Initialize
   	!!! call cplr::pertmod_initialize_(idmodel,rc)

   public :: pertmod_TLinit
   public :: pertmod_TLrun	! interface to ESMF Run with phase=TL (tangent linear)
   public :: pertmod_TLfin
	!!!
   	!!! -- in model_tl() --
	!!! iymd=20000101; ihms=000000
	!!!	! z(0) = 0; y = 0, and this is not a pertmod method
   	!!! call cplr::pertmod_TLinit_(xini(0),y,iymd,ihms,ndtsec,rc)
	!!! y=y+xini(0)
	!!! xobs(0)=y
	!!! do i=0,n-1,+1
	!!!	! y := y(i) in, and y:=M(z(i)+xini(i)) out
   	!!!   call cplr::pertmod_TLrun_(xini(i),y,iymd,ihms,ndt,rc)
	!!!   y=y+xini(i+1)
	!!!   xobs(i+1)=y
	!!!   call tick(iymd,ihms,+010000)
	!!! enddo

   public :: pertmod_ADinit
   public :: pertmod_ADrun	! interface to ESMF Run with phase=AD (adjoint)
   public :: pertmod_ADfin
	!!!
   	!!! -- in model_ad() --
	!!! iymd=20000101; ihms=060000
	!!!	! z(n) = 0; x = 0, and this is not a pertmod method
   	!!! call cplr::pertmod_ADinit_(x,xobs(n),iymd,ihms,ndtsec,rc)
	!!! x=x+xobs(n)
	!!! xini(n)=x
	!!! do i=n-1,0,-1
	!!!   call tick(iymd,ihms,-010000)
	!!!	! (x=x(i)) in, and x:=M''(z(i)+xobs(i+1)) out
   	!!!   call cplr::pertmod_ADrun_(x,xobs(i+1),iymd,ihms,ndt,rc)
	!!!   x=x+xobs(i)
	!!!   xini(i)=x
	!!! enddo

   public :: pertmod_finalize	! interface to ESMF finalize
   	!!! call cplr::pertmod_finalize_(rc)

! !MODULE VARIABLES:
!
! !PUBLIC DATA:
!
	public:: pertmod_internal_tl_state_created
	public:: pertmod_internal_tl_state
	public:: pertmod_internal_ad_state_created
	public:: pertmod_internal_ad_state

	logical,save:: pertmod_internal_tl_state_created
	logical,save:: pertmod_internal_ad_state_created
#ifdef SAVE_DYN_PROG_STATE
	type(dyn_prog),target,save:: pertmod_internal_tl_state
	type(dyn_prog),target,save:: pertmod_internal_ad_state
#else
	type(gsi_bundle),target,save:: pertmod_internal_tl_state
	type(gsi_bundle),target,save:: pertmod_internal_ad_state
#endif

! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  00Jul2010  Guo       First implementation to meet ESMF in half-way.
!  01Apr2010  Treadon   move strip to gridmod
!  15May2010  Todling   Update to use GSI_Bundle
!  20May2011  Guo       This version support either stub_GEOSagcmPert module
!		        and the real GEOS AgcmPert gridded component, and it
!		        has passed GSI built-in adjoint test in GSI evaljgrad().
!  31Aug2010  Guo       Added internal_state variables.
!  28Oct2013  Todling   Rename p3d to prse
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
!EOP
!-------------------------------------------------------------------------

! Vocabulary:
!	pertmod		perturbation model known as g5pert

! module objects
! --------------

   character(len=*),parameter:: myname    ="g5_pertmod"
   character(len=*),parameter:: childname ="g5pert/"

   character(len=*),parameter:: mycfname="GSI_GridComp.rc"

   logical,save :: pertmod_created_ = .false.
   logical,save :: pertmod_inidone_ = .false.
   integer(i_kind),save :: pertmod_initcnt_ = HUGE(1_i_kind)
   logical,save :: pertmod_IDmodel_ = .false.
   integer,save :: pertmod_dtsec_   = R3600

   integer(i_kind),parameter:: GENERIC_ERROR =-9999	! a pertmod generic error flag
   integer(i_kind),parameter:: GSI_SUCCESS   =0		! translated ESMF_SUCCESS, used
   							! only by external interfaces.

!!! Additional module parameters/variables for g5pert/

      character(len=*), parameter :: fnpert = 'fsens.eta.nc4'
      character(len=*), parameter :: fnxgsi = 'xxgsi.eta'

      integer(i_kind),  parameter :: ROOT = 0 ! should really come from above

      integer(i_kind), save :: mycount = 0
      real(r_kind), parameter :: kPa_per_Pa = 0.001_r_kind
      real(r_kind), parameter :: Pa_per_kPa = 1000._r_kind

      logical, parameter :: memtraj = .true.
      logical, parameter :: verbose = .false.

      logical, save :: traj_initzd_ = .false.
      logical, save :: initialized_ = .false.
      logical, save :: skiptraj_    = .false.
      logical, save :: bks_checked_ = .false.
 
      integer(i_kind), save :: vectype_     = 5       ! default is GEOS-5 vector

#include "mytrace.H"
#ifdef DEBUG_MPOUT
#  define _MPOPEN_() call stdout_open(myname)
#  define _MPCLOSE_() call stdout_close()
#else
#  define _MPOPEN_()
#  define _MPCLOSE_()
#endif
   CONTAINS

!-------------------------------------------------------------------------
!BOPI

! !IROUTINE: pertmod_SetServices -- invoke MAPL/ESMF SetServices
! !DESCRIPTION:
! !INTERFACE:

   subroutine pertmod_SetServices ( RC )

   use kinds   , only: i_kind
   use timermod, only: timer_ini,timer_fnl
   use mpeu_util, only: tell,warn,perr,die
   use mpeu_util, only: stdout_open
   use mpeu_util, only: stdout_close
   implicit none

! !ARGUMENTS:

   integer(i_kind), optional, intent(out)    :: RC  ! return code

! !REVISION HISTORY:
!
!EOPI
!-------------------------------------------------------------------------

  character(len=*),parameter:: myname_=myname//"::setServices"
  integer(i_kind):: ier

_MPOPEN_()
_ENTRY_(myname_)
  call timer_ini(myname)
  if(present(rc)) rc=GSI_SUCCESS

  pertmod_created_ = .true.
  pertmod_inidone_ = .false.
  pertmod_initcnt_ = 0

  call timer_fnl(myname)
_EXIT_(myname_)
_MPCLOSE_()
   end subroutine pertmod_SetServices

!-------------------------------------------------------------------
!BOPI

! !IROUTINE: pertmod_initialize -- non-ESMF initialize()
! !DESCRIPTION:
! !INTERFACE:

   subroutine pertmod_initialize ( idmodel, RC )

! !ARGUMENTS:

   use GSI_GridCompMod, only: GSI_AgcmPertGrid
   use kinds    , only: i_kind
   use timermod , only: timer_ini,timer_fnl
   use mpeu_util, only: tell,warn,perr,die
   use mpeu_util, only: stdout_open
   use mpeu_util, only: stdout_close
   		! use the ESMF_Grid already defined by GSI_aGridCompMod
		! for pertmod.
   implicit none
   logical        , optional, intent(in ) :: idmodel
   integer(i_kind), optional, intent(out) :: RC      ! Error code

! !REVISION HISTORY:
!   25May2011  Todling  incorrectly, GSI must call this routine twice; allow 
!                       it to do nothing on 2nd pass (add pertmod_inidone_)
!   26May2011  Guo      Changed to let the lifecycles of some objects to
!			complete, when pertmod_initialize() is referenced
!			more than once (pertmod_inidone_==.true.).
!EOPI
!-------------------------------------------------------------------------
   character(len=*),parameter:: myname_=myname//'::initialize'
   integer(i_kind):: ier

_MPOPEN_()
_ENTRY_(myname_)
  call timer_ini(myname)
  if(present(rc)) rc = GSI_SUCCESS

! Begin...

  	! configuring identity model, if it is required.
  pertmod_IDmodel_= .false.
  if(present(idmodel)) pertmod_IDmodel_=idmodel

  	! Verify the state of pertmod_gc

  if(.not.pertmod_created_) then
    call perr(myname_,'not created, pertmod')
    if(.not.present(rc)) call die(myname_)
    rc=GENERIC_ERROR
    call timer_fnl(myname)
    _EXIT_(myname_)
    _MPCLOSE_()
    return
  endif

  	! This allows pertmod_initialize() being referenced more than
	! once by GSI, without upsetting the lifecycles of some code
	! objects, such as, timer, DEBUG_TRACE, and DEBUG_MPOUT.
  pertmod_initcnt_ = pertmod_initcnt_+1
  if(pertmod_inidone_) then
    if(MAPL_AM_I_ROOT()) &
      call tell(myname_,' an extra call, count =',pertmod_initcnt_)
    call timer_fnl(myname)
    _EXIT_(myname_)
    _MPCLOSE_()
    return
  endif
  	! Now code below will not be executed twice, unless one has
	! completed the lifecycle of g5_pertmod by a call to routine
	! pertmod_finalize().

  pertmod_inidone_=.true.

  if(.not.pertmod_IDmodel_) then
		_TRACE_(myname_,"entering g5_init_()")
    ier=0
    call g5_init_(ier, skiptraj=.not.pertmod_IDmodel_)
		if(ier/=0) then
		  call perr(myname_,'g5_init_(), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
		_TRACE_(myname_,"g5_init_() returned")
  endif

  call timer_fnl(myname)
_EXIT_(myname_)
_MPCLOSE_()
  end subroutine pertmod_initialize

!-------------------------------------------------------------------
!BOPI

! !IROUTINE: pertmod_TLinit -- initialize a TLM
! !DESCRIPTION:
! !INTERFACE:

  subroutine pertmod_TLinit(xi,yo,iymd,ihms,tstep,rc)
!   machine:
!
!$$$  end subprogram documentation block

! !ARGUMENTS:

  use kinds        , only: i_kind
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundleCreate
  use gsi_bundlemod, only: assignment(=)
  use constants    , only: ZERO
  use constants    , only: R3600
  use kinds        , only: i_kind
  use timermod     , only: timer_ini,timer_fnl
  use mpeu_util    , only: tell,warn,perr,die
  use mpeu_util, only: stdout_open
  use mpeu_util, only: stdout_close

  use m_model_tl, only: initial_tl
  implicit none

  type(gsi_bundle),intent(in ):: xi	! input state in x-space
  type(gsi_bundle),intent(out):: yo	! output state in y-space
  integer(i_kind ),intent(in ):: iymd	! initial date (YYYYMMDD) of the perturbation state
  integer(i_kind ),intent(in ):: ihms	! initial time (HHMMSS) of the perturbation state
  integer(i_kind ),intent(out):: tstep	! a single TLM time step in seconds
  integer(i_kind ),optional,intent(out):: rc	! return status code

! !REVISION HISTORY:
!
!EOPI
!-------------------------------------------------------------------------
  character(len=*),parameter :: myname_=myname//'::TLinit'
  character(len=*),parameter :: yoname_='TLpertState'
  integer(i_kind):: ier

_MPOPEN_()
_ENTRY_(myname_)
  call timer_ini(myname)
#ifdef DEBUG_VERBOSE
  call tell(myname_,'at (iymd,ihms) =',(/iymd,ihms/))
#endif
  if(present(rc)) rc=GSI_SUCCESS

  	! Verify the component flag
  if(.not.pertmod_created_) then
    call perr(myname_,'component yet to be created')
    if(.not.present(rc)) call die(myname_)
    rc=GENERIC_ERROR
    call timer_fnl(myname)
    _EXIT_(myname_)
    _MPCLOSE_()
    return
  endif

  tstep=R3600
  if(.not.pertmod_IDModel_) then
    call initial_tl()
    tstep=pertmod_dtsec_	! pertmod_dtsec_ = pdt ... in g5_init_()
  endif

  	! Instanciate xobs, according to xini
  call gsi_bundleCreate(yo,xi,yoname_,ier)	! make a yo as xi
  		if(ier/=0) then
		  call perr(myname_,'gsi_bundleCreate("'//yoname_//'"), istatus =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
  yo=ZERO

  call timer_fnl(myname)
_EXIT_(myname_)
_MPCLOSE_()
end subroutine pertmod_TLinit

!-------------------------------------------------------------------
!BOPI

! !IROUTINE: pertmod_TLrun -- invoke Run() method with phase=pertmod_TL
! !DESCRIPTION:
! !INTERFACE:

   subroutine pertmod_TLrun ( p_xi, yo, iymd,ihms, ndt, rc )

   use m_model_tl    , only: amodel_tl
   use gsi_bundlemod , only: gsi_bundle
   use kinds         , only: i_kind
   use constants     , only: ZERO
   use timermod      , only: timer_ini,timer_fnl
   use mpeu_util     , only: tell,warn,perr,die
   use mpeu_util, only: stdout_open
   use mpeu_util, only: stdout_close
   implicit none

! !ARGUMENTS:

   type(gsi_bundle),               pointer:: p_xi ! input gsi_bundle, maybe null()
   type(gsi_bundle), intent(inout),target :: yo	  ! output gsi_bundle
   integer(i_kind) , intent(in) :: iymd,ihms	! time in (yyyymmdd,hhmmss)
   integer(i_kind) , intent(in) :: ndt	! time step count
   integer(i_kind) , optional,intent(out) :: RC	! Error code

! !REVISION HISTORY:
!
!EOPI
!-------------------------------------------------------------------------

! Locals
   character(len=*),parameter:: myname_=myname//"::TLrun"
   integer(i_kind):: ier
   type(dyn_prog):: xpert

_MPOPEN_()
_ENTRY_(myname_)
  call timer_ini(myname)
  if(present(rc)) rc=GSI_SUCCESS

    ! No action is need when in the IDmodel mode, but excercising the
    ! clock, because of the in-and-out definition of yo.

  if(.not.pertmod_IDmodel_) then
    call prognostics_initial(xpert)
    call gsi2pgcm1_(yo, xpert, 'tlm', ier)
     		if(ier/=GSI_SUCCESS) then
		  call perr(myname_,'gsi2pgcm1_(TL), stat =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif

		_TRACE_(myname_,"entering amodel_tl()")
    call amodel_tl(xpert,nymdi=iymd,nhmsi=ihms,ntsteps=ndt,g5pert=.true.)
		_TRACE_(myname_,"amodel_tl() returned")

   	! put xpert values into yo.
    call pgcm2gsi1_(xpert, yo, 'tlm', ier)
		if(ier/=GSI_SUCCESS) then
		  call perr(myname_,'pgcm2gsi1_(TL), stat =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
    call prognostics_final(xpert)
  endif		! .not.pertmod_IDmodel_

  call timer_fnl(myname)
_EXIT_(myname_)
_MPCLOSE_()
  end subroutine pertmod_TLrun

!-------------------------------------------------------------------
!BOPI

! !IROUTINE: pertmod_TLfin -- finalize ADM
! !DESCRIPTION:
! !INTERFACE:

  subroutine pertmod_TLfin(xi,yo,iymd,ihms,rc)

  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundleDestroy
  use kinds        , only: i_kind
  use timermod     , only: timer_ini,timer_fnl
  use mpeu_util    , only: tell,warn,perr,die
  use mpeu_util    , only: stdout_open
  use mpeu_util    , only: stdout_close
  implicit none

  type(gsi_bundle),intent(in   ):: xi	! untouched perturbation increment
  type(gsi_bundle),intent(inout):: yo	! destroyed perturbation state
  integer(i_kind ),intent(in   ):: iymd	! final date (YYYYMMDD) of the perturbation state
  integer(i_kind ),intent(in   ):: ihms	! final time (HHMMSS) of the perturbation state
  integer(i_kind ),optional,intent(out):: rc	! return status code
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter:: myname_=myname//"::TLfin"
  integer(i_kind):: ier

_MPOPEN_()
_ENTRY_(myname_)
  call timer_ini(myname)
#ifdef DEBUG_VERBOSE
  call tell(myname_,'with (iymd,ihms) =',(/iymd,ihms/))
#endif

  if(present(rc)) rc=GSI_SUCCESS

  call gsi_bundleDestroy(yo,ier)
  		if(ier/=0) then
		  call perr(myname_,'gsi_bundleDestroy(), istatus =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=GENERIC_ERROR
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif

  call timer_fnl(myname)
_EXIT_(myname_)
_MPCLOSE_()
end subroutine pertmod_TLfin

!-------------------------------------------------------------------
!BOPI

! !IROUTINE: pertmod_ADinit -- initialize a ADM
! !DESCRIPTION:
! !INTERFACE:

  subroutine pertmod_ADinit(xo,yi,iymd,ihms,tstep,rc)
!   machine:
!
!$$$  end subprogram documentation block

! !ARGUMENTS:

  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundleCreate
  use gsi_bundlemod, only: assignment(=)
  use constants    , only: ZERO
  use constants    , only: R3600
  use kinds        , only: i_kind
  use timermod     , only: timer_ini,timer_fnl
  use mpeu_util    , only: tell,warn,perr,die
  use mpeu_util    , only: stdout_open
  use mpeu_util    , only: stdout_close
  use m_model_ad   , only: initial_ad
  implicit none

  type(gsi_bundle),intent(out):: xo	! output state to be instanciated in x-space
  type(gsi_bundle),intent(in ):: yi	! input state in y-space
  integer(i_kind ),intent(in ):: iymd	! ending date (YYYYMMDD) of the perturbation state
  integer(i_kind ),intent(in ):: ihms	! ending time (HHMMSS) of the perturbation state
  integer(i_kind ),intent(out):: tstep	! a single ADM time step in seconds
  integer(i_kind ),optional,intent(out):: rc	! return status code

! !REVISION HISTORY:
!
!EOPI
!-------------------------------------------------------------------------
  character(len=*),parameter :: myname_=myname//'::ADinit'
  character(len=*),parameter :: xoname_='ADpertState'
  integer(i_kind):: ier

_MPOPEN_()
_ENTRY_(myname_)
  call timer_ini(myname)
#ifdef DEBUG_VERBOSE
  call tell(myname_,'at (iymd,ihms) =',(/iymd,ihms/))
#endif
  if(present(rc)) rc=GSI_SUCCESS

  	! Verify the component flag
  if(.not.pertmod_created_) then
    call perr(myname_,'component yet to be created')
    if(.not.present(rc)) call die(myname_)
    rc=GENERIC_ERROR
    call timer_fnl(myname)
    _EXIT_(myname_)
    _MPCLOSE_()
    return
  endif

  tstep=R3600
  if(.not.pertmod_IDModel_) then
    call initial_ad()
    tstep=pertmod_dtsec_	! pertmod_dtsec_ = pdt ... in g5_init_()
  endif

  	! Instanciate xini, according to xobs
  call gsi_bundleCreate(xo,yi,xoname_,ier)	! make a xo from yi
  		if(ier/=0) then
		  call perr(myname_,'gsi_bundleCreate("'//xoname_//'"), istatus =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
  xo=ZERO

  call timer_fnl(myname)
_EXIT_(myname_)
_MPCLOSE_()
end subroutine pertmod_ADinit

!-------------------------------------------------------------------
!BOPI

! !IROUTINE: pertmod_ADrun -- invoke registered Run() method with phase=pertmod_AD
! !DESCRIPTION:
! !INTERFACE:

   subroutine pertmod_ADrun ( xo, p_yi, iymd,ihms, ndt, rc )

   use m_model_ad    , only: amodel_ad
   use gsi_bundlemod , only: gsi_bundle
   use kinds         , only: i_kind,r_kind
   use mpimod        , only: mype
   use constants     , only: ZERO
   use timermod      , only: timer_ini,timer_fnl
   use mpeu_util     , only: tell,warn,perr,die
   use mpeu_util, only: stdout_open
   use mpeu_util, only: stdout_close
   implicit none

! !ARGUMENTS:

   type(gsi_bundle), intent(inout),target ::   xo ! output gsi_bundle
   type(gsi_bundle),               pointer:: p_yi ! input gsi_bundle, maybe null()
   integer(i_kind) , intent(in) :: iymd,ihms	! time in (yyyymmdd,hhmmss)
   integer(i_kind) , intent(in) :: ndt	! time step count
   integer(i_kind) , optional,intent(out) :: RC	! Error code

! !REVISION HISTORY:
!
!EOPI
!-------------------------------------------------------------------------

! Locals
   character(len=*),parameter:: myname_=myname//"::ADrun"
   integer(i_kind):: ier
   type(dyn_prog):: xpert

_MPOPEN_()
_ENTRY_(myname_)
  call timer_ini(myname)
  if(present(rc)) rc=GSI_SUCCESS

    ! No action is need when in the IDmodel mode, but excercising the
    ! clock, because of the in-and-out definition of xo.

  if(.not.pertmod_IDmodel_) then
    call prognostics_initial(xpert)
    call gsi2pgcm1_(xo,xpert,'adm',ier)
     		if(ier/=GSI_SUCCESS) then
		  call perr(myname_,'gsi2pgcm1_(AD), stat =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif

		_TRACE_(myname_,"entering amodel_ad()")
    call amodel_ad(xpert,nymdi=iymd,nhmsi=ihms,ntsteps=ndt,g5pert=.true.)
		_TRACE_(myname_,"amodel_ad() returned")

    call pgcm2gsi1_(xpert,xo,'adm',ier)
		if(ier/=GSI_SUCCESS) then
		  call perr(myname_,'pgcm2gsi1_(AD), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
    call prognostics_final(xpert)
  endif		! .not.pertmod_IDmodel_

  call timer_fnl(myname)
_EXIT_(myname_)
_MPCLOSE_()
end subroutine pertmod_ADrun

!-------------------------------------------------------------------
!BOPI

! !IROUTINE: pertmod_ADfin -- finalize ADM
! !DESCRIPTION:
! !INTERFACE:

  subroutine pertmod_ADfin(xo,yi,iymd,ihms,rc)

  use kinds, only: i_kind
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundleDestroy
  use timermod     , only: timer_ini,timer_fnl
  use mpeu_util, only: tell,warn,perr,die
   use mpeu_util, only: stdout_open
   use mpeu_util, only: stdout_close
  implicit none

  type(gsi_bundle),intent(inout):: xo	! destroyed perturbation state
  type(gsi_bundle),intent(in   ):: yi	! untouched perturbation increment
  integer(i_kind ),intent(in   ):: iymd	! final date (YYYYMMDD) of the perturbation state
  integer(i_kind ),intent(in   ):: ihms	! final time (HHMMSS) of the perturbation state
  integer(i_kind ),optional,intent(out):: rc	! return status code
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter:: myname_=myname//"::ADfin"
  integer(i_kind):: ier

_MPOPEN_()
_ENTRY_(myname_)
  call timer_ini(myname)
#ifdef DEBUG_VERBOSE
  call tell(myname_,'with (iymd,ihms) =',(/iymd,ihms/))
#endif

  if(present(rc)) rc=GSI_SUCCESS

  call gsi_bundleDestroy(xo,ier)
  		if(ier/=0) then
		  call perr(myname_,'gsi_bundleDestroy(), istatus =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=GENERIC_ERROR
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif

  call timer_fnl(myname)
_EXIT_(myname_)
_MPCLOSE_()
end subroutine pertmod_ADfin

!-------------------------------------------------------------------
!BOPI

! !IROUTINE: pertmod_finalize -- Finalize method

! !INTERFACE:

   subroutine pertmod_finalize ( RC )

! !ARGUMENTS:

   use kinds, only: i_kind
   use timermod , only: timer_ini,timer_fnl
   use mpeu_util, only: tell,warn,perr,die
   use mpeu_util, only: stdout_open
   use mpeu_util, only: stdout_close
   implicit none
   integer, optional,   intent(  out) :: RC      ! Error code

! !DESCRIPTION: 
!   26May2011  Guo      reset _inidone_ and _initcnt_ here, to let the life-
!			cycle of pertmod to complete, and the lifecycle
!			of pertmod can start again.
!EOPI
!-------------------------------------------------------------------------
! Locals
  character(len=*),parameter:: myname_=myname//'::finalize'
  integer(i_kind):: ier

! Begin...

_MPOPEN_()
_ENTRY_(myname_)
  call timer_ini(myname)
  if(present(rc)) rc=GSI_SUCCESS
 
  pertmod_initcnt_=0
  pertmod_inidone_=.false.	! now one can start from the top, i.e.
  				! call pertmod_initialize().

  if(.not.pertmod_IDmodel_) then
		_TRACE_(myname_,"entering g5_clean_()")
    call g5_clean_()
		_TRACE_(myname_,"g5_clean_() returned")
  endif

  pertmod_IDModel_=.true.

  call timer_fnl(myname)
_EXIT_(myname_)
_MPCLOSE_()
end subroutine pertmod_finalize

#ifdef GEOS_PERT
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pgcm2gsi0_:  Convert gcm-adm/tlm vector to gsi vector
!
! !INTERFACE:

      subroutine pgcm2gsi0_ ( xx, which, stat, &
                              filename, skiptraj, nymd_in, nhms_in )

! !USES:

      implicit none

! !INPUT PARAMETERS:

      character(len=*),           intent(in)  :: which    ! adm or tlm
      character(len=*), optional, intent(in)  :: filename ! name of file w/ perturbation
      logical,          optional, intent(in)  :: skiptraj ! allows skip of trajectory/ics 
                                                                                                                           
! !OUTPUT PARAMETERS:

      type(gsi_bundle),         intent(inout) :: xx       ! GSI increment

      integer(i_kind),            intent(out) :: stat

      integer(i_kind),optional,   intent(out) :: nymd_in
      integer(i_kind),optional,   intent(out) :: nhms_in

! !DESCRIPTION: Convert GEOS-5 perturbation vector to GSI increment vector
!               (as pgcm2gsi1_, but reads GCM perturbation from file)
!
! !REMARKS:
!
!  31Oct2007  Todling   De-activated interpolation capability; not fully coded yet
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  17Jul2007  Todling   Add ability to read in perturbation at diff resolution
!  17Jul2008  Todling   Add filename as optional argument
!  15May2010  Todling   Update to use GSI_Bundle
!
!EOP
!-----------------------------------------------------------------------

     character(len=*), parameter :: myname_ = myname//'*pgcm2gsi0_'
     type(dyn_prog) xpert
     type(dyn_prog) ypert

     character(len=255) :: fname
     integer(i_kind) myimr,myjnp,mynl,mync
     integer(i_kind) ierr,nymdp,nhmsp
     real(r_kind)    dmodel,dgsi

     stat = 0
     xx   = zero

!    Initializes this package
!    ------------------------
     call g5_init_ ( ierr, skiptraj=skiptraj )
        if(ierr/=0) return

!    Set file to be read
!    -------------------
     fname = fnpert
     if(present(filename)) fname = trim(filename)

#ifdef GEOS_PERT 
!    Get dims of incoming vector from input file
!    -------------------------------------------
     call inqpert_dims ( trim(fname), myimr, myjnp, mynl, mync, stat=ierr )

!    Create GCM perturbation vector
!    ------------------------------
     if ( myimr/=imr .or. myjnp/=jnp .or. mynl/=nl .or. mync<nc ) then
          stat = 89
          if (mype==ROOT) then
             print*, 'myimr,myjnp,mynl,mync ', myimr,myjnp,mynl,mync
             print*, '  imr,  jnp,  nl,  nc ',   imr,  jnp,  nl,  nc
             print*, trim(myname_), ': Cannot handle resolution inconsistency '
          endif
          return
     else
         call prognostics_initial ( xpert )
     endif

!    Read in perturbation
!    --------------------
     nymdp = 0; nhmsp = 0
     call getpert ( trim(fname), nymdp, nhmsp, xpert, pick=.false., stat=ierr, vectype=vectype_, forceflip=.true. )
       if(ierr/=0)then
           stat = 90
           if(mype==ROOT) print*, trim(myname_), ': Error retrieving perturbation'
           return
       endif
       dmodel = prognostics_dotp(xpert,xpert)

       if (present(nymd_in)) then
           nymd_in = nymdp
       endif
       if (present(nhms_in)) then
           nhms_in = nhmsp
       endif

!    Convert to GSI perturbation vector
!    ----------------------------------
    if (nlon/=imr .or. nlat/=jnp ) then
        call die ( myname_,': this option is not fully implemented yet' )  ! RT: I am de-activating this for now
        if(mype==ROOT) print*, trim(myname_), ': Interpolating perturbation vector to GSI resolution: '

        call prognostics_initial ( ypert, nlon, nlat, nl, nc )

!       Interpolate input perturbation to internal resolution ...
!       ---------------------------------------------------------
        call pgcm2pgcm_ ( myimr,myjnp,mynl,xpert,  ypert, ierr )

!       ... convert GEOS-4 to GEOS-5 like perturbation and ...
!       ------------------------------------------------------
        if(vectype_==4) call stepon_g4tog5_tl ( nymdp, nhmsp, ypert )

!       ... then convert to GSI
!       -----------------------
        call pgcm2gsi1_ ( ypert, xx, which, stat, jgradf=.true. )

        call prognostics_final ( ypert)

    else

!       Simply convert
!       --------------
        call pgcm2gsi1_ ( xpert, xx, which, stat, jgradf=.true. )

    endif
    dgsi = dot_product(xx,xx)
    if(mype==ROOT) write(6,'(2a,1p,e25.18)') trim(myname_), ': magnitude of input vector in model    space ', dmodel
    if(mype==ROOT) write(6,'(2a,1p,e25.18)') trim(myname_), ': magnitude of input vector in analysis space ', dgsi

!    Release GCM perturbation vector
!    -------------------------------
     call prognostics_final ( xpert )
#endif /* GEOS_PERT */ 

     end subroutine pgcm2gsi0_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pgcm2gsi1_:  Convert gcm-adm/tlm vector to gsi vector
!
! !INTERFACE:

      subroutine pgcm2gsi1_ ( xpert, xx, which, stat, jgradf )

! !USES:

      implicit none

! !INPUT PARAMETERS:

      type(dyn_prog),             intent(in)  :: xpert  ! GCM perturbation vector 
      character(len=*),           intent(in)  :: which  ! adm or tlm
      logical, optional,          intent(in)  :: jgradf ! specify when input is forecast (error) gradient

! !OUTPUT PARAMETERS:

      type(gsi_bundle),         intent(inout) :: xx     ! GSI increment

      integer(i_kind),            intent(out) :: stat

! !DESCRIPTION: Convert GEOS-5 perturbation vector to GSI increment vector
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  21Sep2007  Todling   Handles for O3 and CW.
!  22Feb2008  Todling   Handle for forecast (error) gradient vector.
!  19Nov2008  Todling   Update to use gsi-3d pressure instead of ps.
!  15May2010  Todling   Update to use GSI_Bundle
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*pgcm2gsi1_'

      real(r_kind),  allocatable, dimension(:,:,:) :: sub_u,sub_v,sub_delp,sub_q,sub_tv
      real(r_kind),  allocatable, dimension(:,:,:) :: sub_oz,sub_cw
      real(r_kind),  allocatable, dimension(:,:)   :: sub_ps

      integer(i_kind) i,j,k,ij,ijk
      integer(i_kind) i_u,i_v,i_t,i_q,i_oz,i_cw,i_prse,i_p,i_tsen
      integer(i_kind) ierr,istatus
      logical         scaleit

      scaleit = .true.  ! default: scale input vector as original var from G-5 GCM
      stat = 0
      if ( present(jgradf) ) then
           if(jgradf) scaleit = .false. ! input vector is a gradient, don''t scale vars
      endif

!     Initializes this package
!     ------------------------
      call g5_init_ ( ierr )
        if(ierr/=0) return
      call check_bks()

      allocate (sub_tv  (lat2,lon2,nsig), sub_u (lat2,lon2,nsig),&
                sub_v   (lat2,lon2,nsig), sub_q (lat2,lon2,nsig),&
                sub_delp(lat2,lon2,nsig), sub_ps(lat2,lon2), stat=ierr )
        if ( ierr/=0 ) then
            stat = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_)'
            return
        end if
      if ( nc>1 ) then
          allocate (sub_oz  (lat2,lon2,nsig), stat=ierr )
            if ( ierr/=0 ) then
                stat = 91
                if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_oz)'
                return
            end if
      endif
      if ( nc>2 ) then
          allocate (sub_cw  (lat2,lon2,nsig), stat=ierr )
            if ( ierr/=0 ) then
                stat = 91
                if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_cw)'
                return
            end if
      endif

!     Gather from GCM/Scatter to GSI subdomains
!     -----------------------------------------
                call pert2gsi_ ( xpert%u,          sub_u   , ng_d, ng_s, ierr )
                call pert2gsi_ ( xpert%v,          sub_v   , ng_s, ng_d, ierr )
                call pert2gsi_ ( xpert%pt,         sub_tv  , ng_d, ng_d, ierr )
                call pert2gsi_ ( xpert%delp,       sub_delp,    0,    0, ierr )
                call pert2gsi_ ( xpert%q(:,:,:,1), sub_q   , ng_d, ng_d, ierr )
      if (nc>1) call pert2gsi_ ( xpert%q(:,:,:,2), sub_oz  , ng_d, ng_d, ierr )
      if (nc>2) call pert2gsi_ ( xpert%q(:,:,:,3), sub_cw  , ng_d, ng_d, ierr )
        if ( ierr/=0 ) then
            stat = 99
            if(mype==ROOT) print*, trim(myname_), ': unfinished convertion ...'
            return
        end if

!     Get poiners to GSI state-vector
!     -------------------------------
      ierr=0
      call gsi_bundlegetpointer(xx,'u' ,  i_u,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'v' ,  i_v,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'tv',  i_t,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'q' ,  i_q,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'oz',  i_oz , istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'cw',  i_cw , istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'prse',i_prse,istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'ps',  i_p,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'tsen',i_tsen,istatus);ierr=ierr+istatus
      if(ierr/=0) then
         write(6,*) myname_, ': trouble getting pointers'
         call stop2(999)
      endif

!     Calculate perturbation ps for GSI
!     ---------------------------------
      if (which == 'adm') then
          if ( scaleit ) then
               call ps2delp_ad_ ( Pa_per_kPa ) 
          else
               call delp2ps_    ( one ) 
          endif
      else if (which == 'tlm') then
          call delp2ps_    ( kPa_per_Pa ) 
      else
          call die ( myname_,': invalid option' )
      endif

!     Calculate all other perturbation for GSI
!     ----------------------------------------
      do k=1,nsig
         do j=1,lon2
            do i=1,lat2
               xx%r3(i_u)%q(i,j,k) = sub_u (i,j,k)
               xx%r3(i_v)%q(i,j,k) = sub_v (i,j,k)
               xx%r3(i_t)%q(i,j,k) = sub_tv(i,j,k)
               xx%r3(i_q)%q(i,j,k) = sub_q (i,j,k)
            enddo
         enddo
      enddo
      if ( nc>1 ) then
           do k=1,nsig
              do j=1,lon2
                 do i=1,lat2
                    xx%r3(i_oz)%q(i,j,k) = sub_oz (i,j,k)
                 enddo
              enddo
           enddo
           if(scaleit) xx%r3(i_oz)%q  = xx%r3(i_oz)%q * PPMV2GpG
      endif
      if ( nc>2 ) then
           do k=1,nsig
              do j=1,lon2
                 do i=1,lat2
                    xx%r3(i_cw)%q(i,j,k) = sub_cw (i,j,k)
                 enddo
              enddo
           enddo
      endif
      if ( which == 'tlm' ) then
          call getprs_tl   (xx%r2(i_p)%q,xx%r3(i_t)%q,xx%r3(i_prse)%q)
          call tv_to_tsen  (xx%r3(i_t)%q,xx%r3(i_q)%q,xx%r3(i_tsen)%q)
      endif

!     The following will be left untouched
!     ------------------------------------
!     xx%sst(:) = xx%sst(:)

      deallocate (sub_tv, sub_u, sub_v, sub_q, sub_delp, sub_ps, stat=ierr )
        if ( ierr/=0 ) then
            stat = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_)'
            return
        end if
      if ( nc>1 ) then
          deallocate (sub_oz, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_oz)'
                return
            end if
       endif
      if ( nc>2 ) then
          deallocate (sub_cw, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_cw)'
                return
            end if
       endif

      CONTAINS

      subroutine delp2ps_ ( alpha )  ! delp2p
! inverse
      implicit none
      real(r_kind), intent(in) :: alpha
      real(r_kind) :: bkweight
      real(r_kind),  allocatable, dimension(:,:,:) :: sub_aux
      integer(i_kind) i,j,k
        xx%r2(i_p)%q = zero
        do k=1,nsig
           do j=1,lon2
              do i=1,lat2
                 xx%r2(i_p)%q(i,j) = xx%r2(i_p)%q(i,j) + alpha * sub_delp(i,j,k)
              end do
           end do
        end do
        allocate ( sub_aux(lat2,lon2,nsig+1) )
           k=nsig+1
           do j=1,lon2
              do i=1,lat2
                 sub_aux(i,j,k) = zero
              end do
           end do
        do k=nsig,1,-1
           do j=1,lon2
              do i=1,lat2
                 sub_aux(i,j,k) = sub_aux(i,j,k+1) + alpha * sub_delp(i,j,k)
              end do
           end do
        end do
        do k=1,nsig+1
           do j=1,lon2
              do i=1,lat2
                 xx%r3(i_prse)%q(i,j,k) = sub_aux(i,j,k)
              end do
           end do
        end do
        deallocate ( sub_aux )
      end subroutine delp2ps_

      subroutine ps2delp_ad_ ( alpha )
! adm-only
      implicit none
      real(r_kind), intent(in) :: alpha
      real(r_kind), allocatable :: sub_aux(:,:,:) 
      integer(i_kind) i,j,k
        allocate ( sub_aux(lat2,lon2,nsig+1) )
        sub_aux = zero
        do k=1,nsig
           do j=1,lon2
              do i=1,lat2
                 sub_aux(i,j,k+1) = sub_aux(i,j,k+1) - sub_delp(i,j,k)
                 sub_aux(i,j,k)   = sub_aux(i,j,k)   + sub_delp(i,j,k)
                 sub_delp(i,j,k)  = zero
              enddo
           enddo
        enddo
        do k=1,nsig+1
           do j=1,lon2
              do i=1,lat2
                 xx%r3(i_prse)%q(i,j,k) = xx%r3(i_prse)%q(i,j,k) + alpha * sub_aux(i,j,k)
              enddo
           enddo
        enddo
!??     xx%p = zero
        deallocate ( sub_aux )
      end subroutine ps2delp_ad_

      subroutine pert2gsi_ ( fld, sub, ngd, ngs, stat_ )

      integer(i_kind),intent(in) :: ngd, ngs
      real(r8),       intent(in)  :: fld(:,:,:)
      real(r_kind),   intent(out) :: sub(:,:,:)
      integer(i_kind),intent(out) :: stat_

      character(len=*), parameter :: myname_ = myname//'*pert2gsi_'

      real(r_kind), allocatable :: work4d(:,:,:,:)   ! auxliar 4d array
      real(r_kind), allocatable :: work3d(:,:,:)     ! auxliar 3d array
      real(r_kind), allocatable :: work2d(:,:)       ! auxliar 2d array
      real(r_kind), allocatable :: work(:)
      integer(i_kind) i,j,k,mm1

      mm1 = mype+1 
      stat_ = 0

      allocate ( work3d(nlon,nlat,nsig), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work3d)'
            return
        end if
      allocate ( work4d(imr,jnp,nl,1), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work4d)'
            return
        end if
                                                                                                                           
!     Gather GCM perturbations to root processor
!     ------------------------------------------
      call mp_gather4d(fld, work4d, imr, jnp, nl, 1, jfirst, jlast, 1, nl, ngd, ngs, root)

!     Flip horizontal and vertical
!     ----------------------------
      if ( mype==ROOT ) then
              if (imr/=nlon .or. jnp/=nlat ) then
                  if (which=='adm') then
                      work3d = zero
                      call interpack_terpv_ad ( nlon,nlat,nsig,work3d,work3d, imr,jnp,nl,work4d(:,:,:,1), ierr )
                  else if (which=='tlm') then
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

      end subroutine pgcm2gsi1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: gsi2pgcm0_:  Convert GSI increments to GCM perturbations
!
! !INTERFACE:

      subroutine gsi2pgcm0_ ( nymd, nhms, xx, which, stat, &
                              xp, filename )  ! optionals

! !USES:

      implicit none

! !INPUT PARAMETERS:
                                                                                                                           
      integer(i_kind),    intent(in)    :: nymd   ! date as in YYYYMMDD
      integer(i_kind),    intent(in)    :: nhms   ! time as in HHMMSS
      type(gsi_bundle),   intent(inout) :: xx     ! GSI increment
      character(len=*),   intent(in)    :: which  ! adm or tlm

      character(len=*),optional, intent(in) :: filename ! output filename

! !OUTPUT PARAMETERS:

      type(dyn_prog),optional,intent(out) :: xp
      integer(i_kind),        intent(out) :: stat

! !DESCRIPTION: Convert GSI increment vector to GEOS-5 perturbation vector
!               (as gsi2pgcm1_, but output GCM perturbation to file)
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  30Sep2007  Todling   Updated interface to putpert.
!  17Jan2010  Todling   Update interface to putpert (add forceflip).
!  15May2010  Todling   Update to use GSI_Bundle
!
!EOP
!-----------------------------------------------------------------------

     character(len=*), parameter :: myname_ = myname//'*gsi2pgcm0_'
     type(dyn_prog) xpert

     character(len=255) fname
     integer(i_kind) ierr
     integer(i_kind) idim,jdim,kdim,i,j,k
     real(r_kind), allocatable :: ps(:,:)

     stat = 0

!    Initializes this package
!    ------------------------
     call g5_init_ ( ierr )
        if(ierr/=0) return
     call check_bks()

!    Create GCM perturbation vector
!    ------------------------------
     call prognostics_initial ( xpert )

!    Convert to GSI perturbation vector
!    ----------------------------------
     call gsi2pgcm1_ ( xx, xpert, which, stat )

!    Build surface pressure perturbation - output purposes
!    -----------------------------------------------------
     idim = size(xpert%delp,1)
     jdim = size(xpert%delp,2)
     kdim = size(xpert%delp,3)
     allocate ( ps(idim,jdim) )
     ps = zero
     do k=1,kdim
        do j=1,jdim
           do i=1,idim
              ps(i,j) = ps(i,j) + xpert%delp(i,j,k)
           end do
       end do
     end do

!    Write out perturbation
!    ----------------------
     mycount = mycount + 1
     if (present(filename)) then
       write(fname,'(3a)')      trim(job), '.', trim(filename)
     else
       write(fname,'(4a,i3.3)') trim(job), '.', trim(fnxgsi), '_', mycount
     endif
     call putpert ( job, nymd, nhms, xpert, fvpsasdt, nstep, &
                    ak, bk, Ts, oro, ps, fname, vectype=vectype_, forceflip=.true. )
       if(ierr/=0)then
           stat = 90
           if(mype==ROOT) print*, trim(myname_), ': Error retrieving perturbation'
           return
       endif

     if (present(xp)) then
         call prognostics_dup ( xpert, xp )
     endif

!    Release GCM perturbation vector
!    -------------------------------
     deallocate ( ps )
     call prognostics_final ( xpert )

     end subroutine gsi2pgcm0_

!------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1    !
!------------------------------------------------------------------------
!BOP
!
! !ROUTINE: gsi2pgcm1_:  Convert GSI increments to GCM perturbations
!
! !INTERFACE:

      subroutine gsi2pgcm1_ ( xx, xpert, which, stat )

! !USES:

      implicit none

! !INPUT PARAMETERS:

      type(gsi_bundle),   intent(inout) :: xx    ! GSI increment vector
      character(len=*),   intent(in)    :: which ! adm or tlm

! !OUTPUT PARAMETERS:

      type(dyn_prog),     intent(out)   :: xpert ! GCM perturbation vector
      integer(i_kind),    intent(out)   :: stat  ! return error code

! !DESCRIPTION: Converts GSI increments in to ADM/TLM perturbations.
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  21Sep2007  Todling   Handles for O3 and CW.
!  19Nov2008  Todling   Update to use gsi-3d pressure instead of ps.
!  15May2010  Todling   Update to use GSI_Bundle
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*gsi2pgcm_'

      real(r_kind), allocatable, dimension(:,:,:) :: sub_tv,sub_u,sub_v,sub_q,sub_delp
      real(r_kind), allocatable, dimension(:,:,:) :: sub_oz,sub_cw
      real(r_kind), allocatable, dimension(:,:)   :: sub_ps

      integer(i_kind)  i,j,k,ijk,ij
      integer(i_kind)  i_u,i_v,i_t,i_q,i_oz,i_cw,i_prse,i_p,i_tsen
      integer(i_kind)  ierr,istatus
      character(len=255) :: whatin

      stat = 0

!     Initializes this package
!     ------------------------
      call g5_init_ ( ierr )
        if(ierr/=0) return

      allocate (sub_tv  (lat2,lon2,nsig), sub_u (lat2,lon2,nsig),&
                sub_v   (lat2,lon2,nsig), sub_q (lat2,lon2,nsig),&
                sub_delp(lat2,lon2,nsig), sub_ps(lat2,lon2), stat=ierr )
        if ( ierr/=0 ) then
            stat = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_)'
            return
        end if
      if ( nc>1 ) then
          allocate (sub_oz  (lat2,lon2,nsig), stat=ierr )
            if ( ierr/=0 ) then
                stat = 91
                if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_oz)'
                return
            end if
      endif
      if ( nc>2 ) then
          allocate (sub_cw  (lat2,lon2,nsig), stat=ierr )
            if ( ierr/=0 ) then
                stat = 91
                if(mype==ROOT) print*, trim(myname_), ': Alloc(sub_cw)'
                return
            end if
      endif

!     Get poiners to GSI state-vector
!     -------------------------------
      ierr=0
      call gsi_bundlegetpointer(xx,'u' ,  i_u,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'v' ,  i_v,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'tv',  i_t,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'q' ,  i_q,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'oz',  i_oz , istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'cw',  i_cw , istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'prse',i_prse,istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'ps',  i_p,   istatus);ierr=ierr+istatus
      call gsi_bundlegetpointer(xx,'tsen',i_tsen,istatus);ierr=ierr+istatus
      if(ierr/=0) then
         write(6,*) myname_, ': trouble getting pointers'
         call stop2(999)
      endif

      if (which == 'adm') then
          call tv_to_tsen_ad(xx%r3(i_t)%q,xx%r3(i_q)%q,xx%r3(i_tsen)%q)
          call getprs_ad    (xx%r2(i_p)%q,xx%r3(i_t)%q,xx%r3(i_prse)%q)
      endif

!     Fill in subdomain arrays
!     ------------------------
      do k=1,nsig
         do j=1,lon2
            do i=1,lat2
               sub_u (i,j,k) = xx%r3(i_u)%q(i,j,k)
               sub_v (i,j,k) = xx%r3(i_v)%q(i,j,k)
               sub_tv(i,j,k) = xx%r3(i_t)%q(i,j,k)
               sub_q (i,j,k) = xx%r3(i_q)%q(i,j,k)
            enddo
         enddo
      enddo
      if ( nc>1 ) then
           do k=1,nsig
              do j=1,lon2
                 do i=1,lat2
                    sub_oz(i,j,k) = xx%r3(i_oz)%q(i,j,k)
                 enddo
              enddo
           enddo
           sub_oz = sub_oz / PPMV2GpG
      endif
      if ( nc>2 ) then
           do k=1,nsig
              do j=1,lon2
                 do i=1,lat2
                    sub_cw(i,j,k) = xx%r3(i_cw)%q(i,j,k)
                 enddo
              enddo
           enddo
      endif
      do j=1,lon2
         do i=1,lat2
            sub_ps(i,j) = zero
         enddo
      enddo

!     Calculate perturbation delp
!     ---------------------------
      if (which == 'adm') then
          call delp2ps_ad_ ( kPa_per_Pa ) 
      else if (which == 'tlm') then
          call ps2delp_    ( Pa_per_kPa ) 
      else
          call die ( myname_, ': invalid option' )
      endif

!     Gather from GSI subdomains/Scatter to GCM
!     -----------------------------------------
               call gsi2pert_ ( sub_u,    xpert%u,          ng_d, ng_s, ierr )
               call gsi2pert_ ( sub_v,    xpert%v,          ng_s, ng_d, ierr )
               call gsi2pert_ ( sub_tv,   xpert%pt,         ng_d, ng_d, ierr )
               call gsi2pert_ ( sub_delp, xpert%delp,          0,    0, ierr )
               call gsi2pert_ ( sub_q ,   xpert%q(:,:,:,1), ng_d, ng_d, ierr )
      if(nc>1) call gsi2pert_ ( sub_oz,   xpert%q(:,:,:,2), ng_d, ng_d, ierr )
      if(nc>2) call gsi2pert_ ( sub_cw,   xpert%q(:,:,:,3), ng_d, ng_d, ierr )
        if ( ierr/=0 ) then
            stat = 98
            if(mype==ROOT) print*, trim(myname_), ': unfinished convertion ...'
            return
        end if

      deallocate (sub_tv, sub_u, sub_v, sub_q, sub_delp, sub_ps, stat=ierr )
        if ( ierr/=0 ) then
            stat = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_)'
            return
        end if
      if ( nc>1 ) then
          deallocate (sub_oz, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_oz)'
                return
            end if
       endif
      if ( nc>2 ) then
          deallocate (sub_cw, stat=ierr )
            if ( ierr/=0 ) then
                stat = 99
                if(mype==ROOT) print*, trim(myname_), ': Dealloc(sub_cw)'
                return
            end if
       endif


      CONTAINS

      subroutine ps2delp_ ( alpha )  ! p2delp
! tlm-only
      implicit none
      real(r_kind), intent(in) :: alpha
      real(r_kind) bkweight
      real(r_kind), allocatable, dimension(:,:,:) :: sub_aux
      integer(i_kind) i,j,k,ij,ijk
        allocate ( sub_aux(lat2,lon2,nsig+1) )
        do k=1,nsig+1
           do j=1,lon2
              do i=1,lat2
                 sub_aux(i,j,k) = alpha * xx%r3(i_prse)%q(i,j,k)
              enddo
           enddo
        enddo
        do k=1,nsig
           do j=1,lon2
              do i=1,lat2
                 sub_delp(i,j,k) = sub_aux(i,j,k) - sub_aux(i,j,k+1)
              enddo
           enddo
        enddo
        deallocate ( sub_aux )
      end subroutine ps2delp_

      subroutine delp2ps_ad_ ( alpha )
! inverse-adm
      implicit none
      real(r_kind), intent(in) :: alpha
      real(r_kind), allocatable :: sub_aux(:,:,:)
      integer(i_kind) i,j,k
        allocate ( sub_aux(lat2,lon2,nsig+1) )
        sub_aux=zero
        do k=1,nsig+1
           do j=1,lon2
              do i=1,lat2
                 sub_aux(i,j,k) = sub_aux(i,j,k) + xx%r3(i_prse)%q(i,j,k)
                 xx%r3(i_prse)%q(i,j,k) = zero
              end do
           end do
        end do
        do k=1,nsig
           do j=1,lon2
              do i=1,lat2
                 sub_aux(i,j,k+1) = sub_aux(i,j,k+1) + sub_aux(i,j,k)
                 sub_delp(i,j,k) = sub_delp(i,j,k) + alpha * sub_aux(i,j,k)
                 sub_aux(i,j,k) = zero
              end do
           end do
        end do
        deallocate ( sub_aux )
        do k=nsig,1,-1
           do j=1,lon2
              do i=1,lat2
                 sub_delp(i,j,k) = sub_delp(i,j,k) + alpha * xx%r2(i_p)%q(i,j)
              end do
           end do
        end do
        xx%r2(i_p)%q=zero
      end subroutine delp2ps_ad_

      subroutine gsi2pert_ ( sub, fld, ngd, ngs, stat_ )

      integer(i_kind),intent(in)  :: ngd, ngs
      real(r_kind),   intent(in)  :: sub(:,:,:)
      real(r8),       intent(out) :: fld(:,:,:)
      integer(i_kind),intent(out) :: stat_

      character(len=*), parameter :: myname_ = myname//'*gsi2pert_'

      real(r_kind), allocatable :: fldsm(:,:)
      real(r_kind), allocatable :: work4d(:,:,:,:)   ! auxliar 4d array
      real(r_kind), allocatable :: work3d(:,:,:)     ! auxliar 3d array
      real(r_kind), allocatable :: work(:)

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

      allocate ( work4d(imr,jnp,nl,1), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work4d)'
            return
        end if

!     Flip horizontal and vertical
!     ----------------------------
      if ( mype==ROOT ) then
              if (imr/=nlon .or. jnp/=nlat ) then
                  if (which=='adm') then
                      work4d = zero
                      call interpack_terpv_ad ( imr,jnp,nl,work4d(:,:,:,1),work4d(:,:,:,1), nlon,nlat,nsig,work3d, ierr )
                  else if (which=='tlm') then
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
      call mp_scatter4d ( work4d, fld, imr, jnp, nl, 1, jfirst, jlast, 1, nl, ngd, ngs, root )

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

      end subroutine gsi2pgcm1_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: pgcm2pgcm_:  Convert gcm-adm/tlm vectors between diff resolutions
!
! !INTERFACE:

      subroutine pgcm2pgcm_ ( myimr,myjnp,mynl,ypert,  xpert, stat )

! !USES:

      implicit none

! !INPUT PARAMETERS:

      integer(i_kind), intent(in) :: myimr,myjnp,mynl
      type(dyn_prog),  intent(in) :: ypert  ! incoming GCM perturbation vector

! !OUTPUT PARAMETERS:

      type(dyn_prog), intent(out) :: xpert  ! interpolated GCM perturbation vector

      integer(i_kind),intent(out) :: stat

! !DESCRIPTION: Interpolate and convert GEOS-5 perturbation vector 
!               into internal GEOS-5 perturbation vector.
!
! !REVISION HISTORY:
!
!  17Jul2007  Todling   Initial code.
!
!EOP
!-----------------------------------------------------------------------

      integer(i_kind) n,ierr

      stat = 0

      call pert2pert_ ( ng_d,ng_s, myimr,myjnp,mynl,ypert%u,           xpert%u,          ierr )
      call pert2pert_ ( ng_s,ng_d, myimr,myjnp,mynl,ypert%v,           xpert%v,          ierr )
      call pert2pert_ ( ng_d,ng_d, myimr,myjnp,mynl,ypert%pt,          xpert%pt,         ierr )
      call pert2pert_ (    0,   0, myimr,myjnp,mynl,ypert%delp,        xpert%delp,       ierr )
      do n = 1, nc
      call pert2pert_ ( ng_d,ng_d, myimr,myjnp,mynl,ypert%q(:,:,:,n),  xpert%q(:,:,:,n), ierr )
      enddo

      stat = ierr

      CONTAINS

      subroutine pert2pert_ ( ngd,ngs, myimr,myjnp,mynl,fldi,  fldo, stat_ )

      integer(i_kind),  intent(in)  :: ngd,ngs
      integer(i_kind),  intent(in)  :: myimr,myjnp,mynl
      real(r8), intent(in)  :: fldi(:,:,:)
      real(r8), intent(out) :: fldo(:,:,:)
      integer(i_kind),  intent(out) :: stat_

      character(len=*), parameter :: myname_ = myname//'*pert2pert_'

      real(r_kind), allocatable :: work4di(:,:,:,:)   ! auxliar 4d array
      real(r_kind), allocatable :: work4do(:,:,:,:)   ! auxliar 4d array

      stat_ = 0

      allocate ( work4do(nlon,nlat,nl,1), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work3d)'
            return
        end if
      allocate ( work4di(myimr,myjnp,mynl,1), stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 91
            if(mype==ROOT) print*, trim(myname_), ': Alloc(work4d)'
            return
        end if
                                                                                                                           
!     Gather GCM perturbations to root processor
!     ------------------------------------------
      call mp_gather4d(fldi, work4di, myimr, myjnp, nl, 1, jfirst, jlast, 1, nl, ngd, ngs, root)

!     Interpolate perturbation to internal resolution
!     -----------------------------------------------
      if ( mype==ROOT ) then
           work4do = zero
           call interpack_terpv ( myimr,myjnp,mynl,work4di(:,:,:,1),  nlon,nlat,nl,work4do(:,:,:,1), ierr )
      endif

!     Scatter interpolated perturbations to internal vector
!     -----------------------------------------------------
      call mp_scatter4d ( work4do, fldo, nlon, nlat, nl, nc, jfirst, jlast, 1, nl, ngd, ngs, root )

!     Release work memory
!     -------------------
      deallocate ( work4di, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': Dealloc(work4di)'
            return
        end if
      deallocate ( work4do, stat=ierr )
        if ( ierr/=0 ) then
            stat_ = 99
            if(mype==ROOT) print*, trim(myname_), ': delloc(work4do)'
            return
        end if

      end subroutine pert2pert_

      end subroutine pgcm2pgcm_

#endif /* GEOS_PERT */

!------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1    !
!------------------------------------------------------------------------
!BOP
!
! !ROUTINE: g5_init_:  Initializes GEOS-5 TLM/ADM
!
! !INTERFACE:

      subroutine g5_init_ ( stat, skiptraj )

      implicit none

! !INPUT PARAMETERS:

!     integer(i_kind),   intent(in) :: nymd
!     integer(i_kind),   intent(in) :: nhms
      logical, optional, intent(in) :: skiptraj     ! when .t., trajectory not read in
      
! !OUTPUT PARAMETERS:
  
      integer(i_kind), intent(out) :: stat

! !DESCRIPTION: Initializes GEOS-5 TLM and ADM
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!  19Nov2008  Todling   Allow upto 3 tracers.
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*g5_init_'
#ifdef GEOS_PERT
      type(dyn_prog) :: prog
      integer(i_kind) m,n
#endif /* GEOS_PERT */

      stat = 0

#ifdef GEOS_PERT

!     If already initialized, there is nothing to do
!     ----------------------------------------------
      if(initialized_) return
      if(present(skiptraj)) skiptraj_ = skiptraj

!     Consistency checking between ADM/TLM and GSI dims
!     -------------------------------------------------
      if ( nc>3 ) then
           stat = 90
           if(mype==ROOT) print*, trim(myname_), ': unacceptable number of tracers'
           return
      endif
      if ( nsig/=nl ) then
           stat = 91
           if(mype==ROOT) print*, trim(myname_), ': inconsistent number of levels, nsig,nl: ', nsig,nl
           return
      endif

      call setfunc ( n, m, prog )

      call stepon_set ( prog, rstskip=skiptraj )

      call prognostics_final ( prog )

!     Set public time step of TL and AD models
!     ----------------------------------------
      pertmod_dtsec_ = pdt

      if ( .not. skiptraj_ ) then

!         Initialize dynamics trajectory handle
!         -------------------------------------
          call getstate_init ( nymd, nhms, memtrj=memtraj, verbose=verbose )
          call getstate      ( nymd, nhms )

!         Initialize physics trajectory handle
!         -------------------------------------
          call physdrv1_get_init ( nymd, nhms, memphys=memtraj, verbose=verbose )
          call physdrv1_get_all  ( nymd, nhms )

          traj_initzd_ = .true.

      endif
#else /* GEOS_PERT */

!     Set public DUMMY time step of TL and AD models
!     ----------------------------------------------
      pertmod_dtsec_ = R3600

#endif /* GEOS_PERT */

      initialized_ = .true.

      end subroutine g5_init_

!------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 601.1    !
!------------------------------------------------------------------------
!BOP
!
! !ROUTINE: g5_clean_:  Clean GEOS-5 TLM/ADM
!
! !INTERFACE:

      subroutine g5_clean_ ( )

      implicit none

! !INPUT PARAMETERS:
      
! !OUTPUT PARAMETERS:
  

! !DESCRIPTION: Initializes GEOS-5 TLM and ADM
!
! !REVISION HISTORY:
!
!  08May2007  Todling   Initial code.
!
!EOP
!-----------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'*g5_clean_'
#ifdef GEOS_PERT
      type(dyn_prog) :: prog
#endif /* GEOS_PERT */

!     If not initialized, there is nothing to do
!     ------------------------------------------
      if(.not.initialized_) return

#ifdef GEOS_PERT

!_RT  call postfunc ( prog )

!_RT  call prognostics_final ( prog )

      if ( traj_initzd_ ) then

!         Initialize dynamics trajectory handle
!         -------------------------------------
          call getstate_clean ( )

!         Initialize physics trajectory handle
!         -------------------------------------
          call physdrv1_get_clean ( )
          traj_initzd_ = .false.

      endif

#endif /* GEOS_PERT */

      initialized_ = .false.

      end subroutine g5_clean_

    subroutine check_bks
    implicit none
#ifdef GEOS_PERT
    integer(i_kind) k,kk 
    if(bks_checked_) return
    do k = 1, nsig+1
       kk = nsig-k+2
       if(abs(bk(k)-bk5(kk))>0.00001_r_kind)then
          if(mype==ROOT)then
             print*,'bk5',bk5
             print*,'bk',bk
          endif
          write(6,*)'check_bks: troubled vertical coord system'
          call stop2(126)
       endif
    enddo
#endif /* GEOS_PERT */
    bks_checked_ = .true.
    end subroutine check_bks


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
  
end module g5_pertmod
