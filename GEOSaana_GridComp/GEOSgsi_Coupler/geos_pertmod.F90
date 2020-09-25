#include "MAPL_ErrLog.h"
!#define DEBUG_VERBOSE
!#define DEBUG_TRACE
!#define DEBUG_CLOCK
!#define DEBUG_STATES
!#define DEBUG_MPOUT
!#define DEBUG_METACOMP

#define FORCING
	! If defined(FORCING), imports of RUN() are increments to the
	! perturbation states. Otherwise, imports of RUN() are perturbation
	! states themselves.

!#define STUB_AGCMPERT
	! If define(STUB_AGCMPERT), a local stub_GEOSagcmPert_GridCompMod
	! is used to simulate an identity perturbation model.
#ifdef STUB_AGCMPERT
#  ifndef FORCING
#     define FORCING
#  endif
#endif

   module geos_pertmod
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: geos_pertmod --- GEOS/GSI perturbation model component 
!
! !DESCRIPTION: This gridded component provides a wrapper around the
!               GEOSagcmPert_GridCompMod, as well as a placeholder for
!		default user geos_pertmod data elements.
!
! !USES:

   use ESMF, only: ESMF_GridComp
   use ESMF, only: ESMF_State
   use ESMF, only: ESMF_Clock
   use ESMF, only: ESMF_Config
   use kinds, only: i_kind,r_kind

#if defined(STUB_AGCMPERT) || defined(G5_PERTMOD)
   use stub_GEOSagcmPert_GridCompMod, only: agcmPert_setServices => setServices
   use stub_GEOSagcmPert_GridCompMod, only: agcmPert_TLMphase    => TLMphase
   use stub_GEOSagcmPert_GridCompMod, only: agcmPert_ADMphase    => ADMphase
   use stub_GEOSagcmPert_GridCompMod, only: phase
   use stub_GEOSagcmPert_GridCompMod, only: phase_opAdjoint
   use stub_GEOSagcmPert_GridCompMod, only: phase_addImport
   use stub_GEOSagcmPert_GridCompMod, only: phase_getExport

#else
   use GEOS_AgcmPertGridCompMod, only: agcmPert_setServices => setServices
   use GEOS_PertSharedMod, only: agcmPert_TLMphase => TLMphase
   use GEOS_PertSharedMod, only: agcmPert_ADMphase => ADJphase
   use GEOS_PertSharedMod, only: phase
   use GEOS_PertSharedMod, only: phase_opAdjoint
   use GEOS_PertSharedMod, only: phase_addImport
   use GEOS_PertSharedMod, only: phase_getExport

   use GEOS_PertSharedMod, only: pert_ak,pert_bk
   use gridmod,            only: ak5,bk5,nsig

#endif
#ifdef DEBUG_METACOMP
   use MAPL, only: MAPL_MetaCompPrint
#endif

   use geos_pertState, only: pertState_TLMvect => TLMvect
   use geos_pertState, only: pertState_ADMvect => ADMvect

!   use mpeu_util, only: tell,warn,perr,die

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

!!   public :: pertmod_get	! inquire for certain pertmod_gc parameters.
!!   public :: pertmod_set	! configurate certain pertmod_gc parameters.
!!   public :: pertmod_create	! create the pertmod_gc
!!   public :: pertmod_destroy	! destroy the pertmod_gc

	!!! -- in GSI_GridComp_finalize()//cplr::finalize(rc=..) --
   	!!! use pertmod, only:: pertmod_finalize

!   call pertmod_setServices(rc)
!
!   call pertmod_initialize(idmodel,rc)
!
!   call clock_set(iymd,ihms,begin_iymd,begin_ihms)
!   call pertmod_TLinit(xini(0),yo,iymd,ihms,tstep,rc)
!   do i=1,n,+1
!     call clock_tick(iymd,ihms,+ntstep*tstep)
!     call pertmod_TLrun (p_xini,yo,iymd,ihms,ntstep,rc)
!   enddo
!   call pertmod_TLfin (p_xini,xobs,iymd,ihms,rc)
!
!   call pertmod_ADinit(p_xini,p_xobs,iymd,ihms,tstep,rc)
!   do i=n,1,-1
!     call pertmod_ADrun(xini,p_xobs,iymd,ihms,ntstep,rc)
!     call clock_tick(iymd,ihms,-ntstep*tstep)
!   enddo
!   call pertmod_ADfin(xini,p_xobs,iymd,ihms,rc)
!   call clock_end(iymd,ihms,begin_iymd,begin_ihms)
!
!   call pertmod_finalize(rc)

! !MODULE VARIABLES:

!!   public :: pertmod_gc		! type( ESMF_GridComp )
!!   public :: pertmod_im		! type( ESMF_State    )
!!   public :: pertmod_ex		! type( ESMF_State    )
!!   public :: pertmod_clock	! type( ESMF_Clock    )

! !REVISION HISTORY:
!
!  00Jul2010  Guo  First implementation to meet ESMF in half-way.
!  20May2011  Guo  This version support either stub_GEOSagcmPert module
!		   and the real GEOS AgcmPert gridded component, and it
!		   has passed GSI built-in adjoint test in GSI evaljgrad().
!
!EOP
!-------------------------------------------------------------------------

! Vocabulary:
!	pertmod		perturbation model known as GEOSagcmPert

! module objects
! --------------
   type(ESMF_GridComp), target, save :: pertmod_gc	! myself as a ESMF GridComp
   type(ESMF_State)   , target, save :: pertmod_im	! my ESMF import state
   type(ESMF_State)   , target, save :: pertmod_ex	! my ESMF export state
   type(ESMF_Clock)   , target, save :: pertmod_clock	! my ESMF clock

!!   type(ESMF_Clock),allocatable, save :: initial_pertmod_clock_ ! the initial clock
   type(ESMF_Clock),allocatable, save :: pertmod_clock_forward_ ! a clock_forward
   type(ESMF_Clock),allocatable, save :: pertmod_clock_reverse_ ! a clock_reverse

   type(ESMF_Config)  , target, save :: pertmod_config

   character(len=*),parameter:: myname    ="geos_pertmod"
#ifdef STUB_AGCMPERT
   character(len=*),parameter:: childname ="stub_AGCMPERT/"
#else
   character(len=*),parameter:: childname ="GEOS_AgcmPert/"
#endif
!  character(len=*),parameter:: mycfname="agcmPert_GridComp.rc"
   character(len=*),parameter:: mycfname="GSI_GridComp.rc"

   logical,parameter :: pertState_testson_ = .false.

   logical,save :: pertmod_created_ = .false.
   logical,save :: pertmod_inidone_ = .false.
   integer(i_kind),save :: pertmod_initcnt_ = HUGE(1_i_kind)
   logical,save :: pertmod_IDmodel_ = .false.

   integer(i_kind),parameter:: GENERIC_ERROR =-9999	! a pertmod generic error flag
   integer(i_kind),parameter:: GSI_SUCCESS   =0		! translated ESMF_SUCCESS, used
   							! only by external interfaces.

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

   use ESMF, only: ESMF_SUCCESS
   use ESMF, only: ESMF_GridCompCreate
!   use ESMF, only: ESMF_ATM
   use ESMF, only: ESMF_Config
   use ESMF, only: ESMF_ConfigCreate
   use ESMF, only: ESMF_ConfigLoadFile
   use ESMF, only: ESMF_ConfigGetAttribute

   use ESMF, only: ESMF_GridCompSetServices
!!   use MAPL, only: MAPL_GenericSetServices
   use MAPL, only: MAPL_am_I_root
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
  character(len=80):: expid

_MPOPEN_()
_ENTRY_(myname_)
  call timer_ini(myname)
  if(present(rc)) rc=GSI_SUCCESS
  if(GSI_SUCCESS/=ESMF_SUCCESS) then
    call warn(myname_,'mismatched return code')
    call warn(myname_,'         GSI_SUCCESS =', GSI_SUCCESS)
    call warn(myname_,'        ESMF_SUCCESS =',ESMF_SUCCESS)
  endif

! Begin...

  	! MAPL-Rule 4: a MAPL GridComp will expect a configuration to be
	! open when setServices() is called.  Do I need to open(configfile)?

!  Create and load the configuration file
!  ----------------------------------
  pertmod_config = ESMF_ConfigCreate (rc=ier)
		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,"ESMF_ConfigCreate(), rc =",ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif

  call ESMF_ConfigLoadFile   ( pertmod_config, mycfname,rc=ier )
   		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_ConfigLoadFile("'//trim(mycfname)//'"), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif

  call ESMF_ConfigGetAttribute( pertmod_config, expid, label ='expid:', default='expid:', rc=ier)
   		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_ConfigGetAttribute(expid), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif

  pertmod_gc = ESMF_GridCompCreate(	name = childname     , &
!  				GridCompType = ESMF_ATM	     , &
				      Config = pertmod_config, rc=ier)
   		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,"ESMF_GridCompCreate(), rc =",ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
   	
  pertmod_created_ = .true.
  pertmod_inidone_ = .false.
  pertmod_initcnt_ = 0

! SetServices() -- setservices() for pertmod_gc itself only.
! -------------------------------------------
	! It is expected that, GEOSagcmPert_setServices() will add both
	! import-specs and export-specs into the pertmod_gc internal state.

		_TRACE_(myname_,"entering ESMF_GridCompSetServices(agcmPert_SetServices)")
   call ESMF_GridCompSetServices ( pertmod_gc, agcmPert_SetServices, rc=ier)
   		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,"ESMF_GridCompSetServices(), rc =",ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
		_TRACE_(myname_,"ESMF_GridCompSetServices(agcmPert_SetServices) returned")

#ifdef DEBUG_METACOMP
  call MAPL_MetaCompPrint(pertmod_gc,header='after ESMF_GridCompSetServices(agcmPert_SetServices)',deep=.true.)
#endif
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

   use ESMF, only: ESMF_SUCCESS
   use ESMF, only: ESMF_GridCompSet, ESMF_GridValidate
   use ESMF, only: ESMF_GridCompInitialize
   use ESMF, only: ESMF_StateCreate, ESMF_STATEINTENT_IMPORT, ESMF_STATEINTENT_EXPORT
   use ESMF, only: ESMF_StatePrint
   use ESMF, only: ESMF_ClockCreate
!!   use MAPL, only: MAPL_GenericInitialize
   use MAPL, only: WRITE_PARALLEL
   use MAPL, only: MAPL_am_I_root

   use geos_pertState , only: pertState_setNeeded
   use GSI_GridCompMod, only: GSI_AgcmPertGrid
   use kinds    , only: i_kind,r_kind
   use timermod , only: timer_ini,timer_fnl
   use mpeu_util, only: tell,warn,perr,die
   use mpeu_util, only: stdout_open
   use mpeu_util, only: stdout_close
   		! use the ESMF_Grid already defined by GSI_GridCompMod
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
#ifdef DEBUG_METACOMP
  call MAPL_MetaCompPrint(pertmod_gc,deep=.true.)
#endif
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
	! completed the lifecycle of geos_pertmod by a call to routine
	! pertmod_finalize().

  pertmod_inidone_=.true.

  	! Initialize the GridComp (pertmod_gc) with the same grid
	! defined for GSI_GridComp.

  call ESMF_GridCompSet(pertmod_gc, grid=GSI_AgcmPertGrid, rc=ier)
  		if(ier/=0) then
		  call perr(myname_,'ESMF_GridCompSet(), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif

   	! Create its import state
   pertmod_im = ESMF_StateCreate (name="agcmPert_Import",  &
        stateintent = ESMF_STATEINTENT_IMPORT, rc = ier)
		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_StateCreate(ESMF_STATEINTENT_IMPORT), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif

   	! Create its export state
   pertmod_ex = ESMF_StateCreate (name="agcmPert_Export",  &
        stateintent = ESMF_STATEINTENT_EXPORT, rc = ier)
		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_StateCreate(ESMF_STATEINTENT_EXPORT), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif

  	! Create a clock for intialization, pertmod_clock_forward_.

  call myClock_create_(pertmod_clock)

!  	! Envoke Initialize_(), which has been registered previously by
!	! SetServicies_()
!
!  call MAPL_GenericInitialize( gc=pertmod_gc,	&
!   		import=pertmod_im,	&
!		export=pertmod_ex,	&
!		 clock=pertmod_clock, rc=ier)
!
!	  	if(ier/=0) then
!		  call perr(myname_,'MAPL_GenericInitialize(), rc =',ier)
!		  if(.not.present(rc)) call die(myname_)
!		  rc=ier
!		  call timer_fnl(myname)
!		  _EXIT_(myname_)
!		  _MPCLOSE_()
!		  return
!		endif

  if(.not.pertmod_IDmodel_) then
		_TRACE_(myname_,"entering ESMF_GridCompInitialize()")
    call ESMF_GridCompInitialize(pertmod_gc, importState=pertmod_im, &
         exportState=pertmod_ex, clock=pertmod_clock, rc=ier)
		if(ier/=0) then
		  call perr(myname_,'ESMF_GridCompInitialize(), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
		_TRACE_(myname_,"ESMF_GridCompInitialize() returned")
    call pertState_setNeeded(pertmod_ex)
  endif

#ifdef DEBUG_VERBOSE
  call WRITE_PARALLEL ( myname_//": IMPORT State" )
	if ( MAPL_am_I_root() ) call ESMF_StatePrint ( pertmod_im )
  call WRITE_PARALLEL ( myname_//": EXPORT State" )
	if ( MAPL_am_I_root() ) call ESMF_StatePrint ( pertmod_ex )
#endif

! Initialize ak/bk - not sure this is the best place for this!
  allocate(pert_ak(0:nsig),pert_bk(0:nsig)) 
  pert_ak(0:nsig)=ak5(nsig+1:1:-1)*1000. ! back to hPa
  pert_bk(0:nsig)=bk5(nsig+1:1:-1)

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

  use ESMF     , only: ESMF_SUCCESS
  use ESMF     , only: ESMF_GridCompSet
  use MAPL     , only: MAPL_am_I_root
  use kinds        , only: i_kind
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundleCreate
  use gsi_bundlemod, only: assignment(=)
  use constants    , only: ZERO
  use kinds        , only: i_kind
  use timermod     , only: timer_ini,timer_fnl
  use mpeu_util    , only: tell,warn,perr,die
  use mpeu_util, only: stdout_open
  use mpeu_util, only: stdout_close
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

  	! Config forward pertmod_clock, where
	!   1) the begin time equals the current time (iymd,ihms);
	!   2) the model time step is forward (tstep>0);
	!   3) return tstep.

  call myClock_TLconfig_(pertmod_clock,iymd,ihms,tstep,rc=ier)
  		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'myClock_TLconfig_(), rc =',ier)
		  call perr(myname_,'myClock_TLconfig_(), iymd =',iymd)
		  call perr(myname_,'myClock_TLconfig_(), ihms =',ihms)
		  call perr(myname_,'myClock_TLconfig_(), tstep =',tstep)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
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

  call myClock_show_(pertmod_clock,myname_,'pertmod_clock',show=.false.)

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

   use ESMF, only: ESMF_SUCCESS
   use ESMF, only: ESMF_GridCompRun
   use ESMF, only: ESMF_ClockAdvance
   use MAPL, only: MAPL_am_I_root

   use geos_pertState, only: pertState_set	! set a pertState to a value or another pertState
   use geos_pertState, only: pertState_show	! show a pertState
   use geos_pertState, only: pertState_getIFx	! get a pertState from xi
   use geos_pertState, only: pertState_putIFx	! put a pertState into yo

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
   type(gsi_bundle),pointer:: p_yo
   integer(i_kind):: ier,istep
   integer(i_kind):: this_phase
   logical:: yesImport,yesExport

_MPOPEN_()
_ENTRY_(myname_)
  call timer_ini(myname)
  if(present(rc)) rc=GSI_SUCCESS

    ! No action is need when in the IDmodel mode, but excercising the
    ! clock, because of the in-and-out definition of yo.

  yesImport=associated(p_xi)
  yesExport=.true.

  if(.not.pertmod_IDmodel_) then
#ifdef FORCING
    call pertState_getIFx(pertmod_im,p_xi, pertState_TLMvect, rc=ier)
#else /* FORCING */
	.. udef(FORCING) is not supported ..
    p_yo => yo
    call pertState_getIFx(pertmod_im,p_yo, pertState_TLMvect, rc=ier)
    p_yo => NULL()
#endif /* FORCING */
     		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'pertState_getIFx(TL), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif

  endif		! .not.pertmod_IDmodel_

      	! Invoke my registered Run() method
  do istep =1,ndt,+1
    this_phase = phase(adjoint=.false.,	import=yesImport.and.istep==1, &
					export=yesExport.and.istep==ndt)
	_TRACEV_(myname_,'istep =',istep)
	_TRACEV_(myname_,'phase =',this_phase)

    if(.not.pertmod_IDmodel_) then
#ifdef DEBUG_STATES
      if(phase_addImport(this_phase)) &
      	call pertState_show(pertmod_im, myname_,'pertmod_im')
#endif
		_TRACE_(myname_,"entering ESMF_GridCompRun(TL)")
      call ESMF_GridCompRun(pertmod_gc, importState=pertmod_im, &
                          exportState=pertmod_ex, &
    			  clock=pertmod_clock, &
			  phase=this_phase   , rc=ier)

		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_GridCompRun(TL), rc =',ier)
		  call perr(myname_,'ESMF_GridCompRun(TL), iymd =',iymd)
		  call perr(myname_,'ESMF_GridCompRun(TL), ihms =',ihms)
		  call perr(myname_,'ESMF_GridCompRun(TL), istep =',istep)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
		_TRACE_(myname_,"ESMF_GridCompRun(TL) returned")
#ifdef DEBUG_STATES
      if(phase_getExport(this_phase)) &
      	call pertState_show(pertmod_ex, myname_,'pertmod_ex')
#endif

#ifdef FORCING
      call pertState_set(pertmod_im,ZERO, rc=ier)
#else /* FORCING */
	.. udef(FORCING) is not supported ..
      call pertState_set(pertmod_im,pertmod_ex, rc=ier)
#endif /* FORCING */
		! set pertmod_im=0, so it has no effect on the next time step.
     		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'pertState_set(pertmod_im,0.), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
    endif		! .not.pertmod_IDmodel_

    call ESMF_ClockAdvance(pertmod_clock)
    call myClock_show_(pertmod_clock,myname_,'pertmod_clock',show=.false.)
  enddo

  if(.not.pertmod_IDmodel_) then

   	! put pertmod_ex values into yo, unless p_yo is null().
    p_yo => yo
    call pertState_putIFx(pertmod_ex,p_yo, pertState_TLMvect, rc=ier)
		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'pertState_putIFx(TL), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
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

  use ESMF     , only: ESMF_SUCCESS
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundleDestroy
  use kinds        , only: i_kind
  use timermod     , only: timer_ini,timer_fnl
  use mpeu_util    , only: tell,warn,perr,die
   use mpeu_util, only: stdout_open
   use mpeu_util, only: stdout_close
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

  use ESMF     , only: ESMF_SUCCESS
  use ESMF     , only: ESMF_GridCompSet
  use gsi_bundlemod, only: gsi_bundle
  use gsi_bundlemod, only: gsi_bundleCreate
  use gsi_bundlemod, only: assignment(=)
  use constants    , only: ZERO
  use kinds        , only: i_kind
  use timermod     , only: timer_ini,timer_fnl
  use mpeu_util    , only: tell,warn,perr,die
  use mpeu_util    , only: stdout_open
  use mpeu_util    , only: stdout_close
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

  	! Config backward pertmod_clock, where
	!   1) the end time equals the current time (iymd,ihms);
	!   2) the model time step is set backward (tstep<0);
	!   3) return tstep.

  call myClock_ADconfig_(pertmod_clock,iymd,ihms,tstep,rc=ier)
  		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'myClock_ADconfig_(), rc =',ier)
		  call perr(myname_,'myClock_ADconfig_(), iymd =',iymd)
		  call perr(myname_,'myClock_ADconfig_(), ihms =',ihms)
		  call perr(myname_,'myClock_ADconfig_(), tstep =',tstep)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
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

   use ESMF, only: ESMF_SUCCESS
   use ESMF, only: ESMF_GridCompRun
   use ESMF, only: ESMF_ClockAdvance
   use MAPL, only: MAPL_am_I_root
   use geos_pertState, only: pertState_set	! set a pertState to a value or another pertState
   use geos_pertState, only: pertState_show	! show a pertState
   use geos_pertState, only: pertState_getIFx	! get a pertState from yi
   use geos_pertState, only: pertState_putIFx	! put a pertState into xo
   use geos_pertState, only: pertState_dot_product
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
   type(gsi_bundle),pointer:: p_xo
   integer(i_kind):: ier,istep
   integer(i_kind):: this_phase
   real(r_kind) dp
   logical:: yesImport,yesExport

_MPOPEN_()
_ENTRY_(myname_)
  call timer_ini(myname)
  if(present(rc)) rc=GSI_SUCCESS

    ! No action is need when in the IDmodel mode, but excercising the
    ! clock, because of the in-and-out definition of xo.

  yesImport=associated(p_yi)
  yesExport=.true.

  if(.not.pertmod_IDmodel_) then
   	! get pertmod_im values from yi, or set to zero if null()
#ifdef FORCING
    call pertState_getIFx(pertmod_im,p_yi, pertState_ADMvect, rc=ier)
#else /* FORCING */
	.. udef(FORCING) is not supported ..
    p_xo => xo
    call pertState_getIFx(pertmod_im,p_xo, pertState_ADMvect, rc=ier)
    p_xo => NULL()
    dp=pertState_dot_product(pertmod_im)
    if(mype==0) print*, '_RT debug inside ADrun in  dp(xo) = ', dp
#endif /* FORCING */
     		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'pertState_getIFx(AD), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif

  endif		! .not.pertmod_IDmodel_

      	! Invoke Run() method of the PGCM Grid Component.
  do istep=ndt,1,-1
    this_phase = phase(adjoint=.true.,	import=yesImport.and.istep==ndt, &
					export=yesExport.and.istep==1    )
	_TRACEV_(myname_,'istep =',istep)
	_TRACEV_(myname_,'phase =',this_phase)

    call myClock_show_(pertmod_clock,myname_,'pertmod_clock',show=.false.)
    call ESMF_ClockAdvance(pertmod_clock)

    if(.not.pertmod_IDmodel_) then
#ifdef DEBUG_STATES
      if(phase_addImport(this_phase)) &
      	call pertState_show(pertmod_im, myname_,'pertmod_im')
#endif
		_TRACE_(myname_,"entering ESMF_GridCompRun(AD)")
      call ESMF_GridCompRun(pertmod_gc, importState=pertmod_im, &
                        exportState=pertmod_ex, &
			clock=pertmod_clock, &
			phase=this_phase   , rc=ier)
	     	if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_GridCompRun(AD), rc =',ier)
		  call perr(myname_,'ESMF_GridCompRun(AD), iymd =',iymd)
		  call perr(myname_,'ESMF_GridCompRun(AD), ihms =',ihms)
		  call perr(myname_,'ESMF_GridCompRun(AD), istep =',-istep)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
		_TRACE_(myname_,"ESMF_GridCompRun(AD) returned")
#ifdef DEBUG_STATES
      if(phase_getExport(this_phase)) &
        call pertState_show(pertmod_ex, myname_,'pertmod_ex')
#endif

#ifdef FORCING
    call pertState_set(pertmod_im,ZERO, rc=ier)
#else /* FORCING */
	.. udef(FORCING) is not supported ..
    dp=pertState_dot_product(pertmod_ex)
    if(mype==0) print*, '_RT debug inside ADrun in before  dp(ex) = ', dp
    call pertState_set(pertmod_im,pertmod_ex, rc=ier)
    dp=pertState_dot_product(pertmod_im)
    if(mype==0) print*, '_RT debug inside ADrun in after   dp(im) = ', dp
#endif /* FORCING */
     		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'pertState_set(pertmod_im,ZERO), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
    endif		! .not.pertmod_IDmodel_
  enddo

   	! put pertmod_ex values into xo, unless p_xo is null()
  if(.not.pertmod_IDmodel_) then
    p_xo => xo
    call pertState_putIFx(pertmod_ex,p_xo, pertState_ADMvect, rc=ier)
		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'pertState_putIFx(AD), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
!   dp=pertState_dot_product(pertmod_ex)
!   if(mype==0) print*, '_RT debug inside ADrun out dp(xo) = ', dp
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

  use ESMF, only: ESMF_SUCCESS
  use MAPL, only: MAPL_am_I_root
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

  call myClock_show_(pertmod_clock,myname_,'pertmod_clock',show=.false.)

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

   use ESMF, only: ESMF_SUCCESS
   use ESMF, only: ESMF_GridCompDestroy
   use ESMF, only: ESMF_GridCompFinalize
   use ESMF, only: ESMF_StateDestroy
   use ESMF, only: ESMF_ConfigDestroy
   use ESMF, only: myClock_destroy_ => ESMF_ClockDestroy

!!   use MAPL, only: MAPL_GenericFinalize
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
		_TRACE_(myname_,"entering ESMF_GridCompFinalize()")
    call ESMF_GridCompFinalize(pertmod_gc, importState=pertmod_im, &
         exportState=pertmod_ex, clock=pertmod_clock, rc=ier)
  		if(ier/=0) then
		  call perr(myname_,'ESMF_GridCompFinalize(), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
		_TRACE_(myname_,"ESMF_GridCompFinalize() returned")
  endif

!  -------------------------------------
!  call MAPL_GenericFinalize ( gc=pertmod_gc,	&
!   		import=pertmod_im,	&
!		export=pertmod_ex,	&
!		 clock=pertmod_clock, rc=ier)
!
!		if(ier/=ESMF_SUCCESS) then
!		  call perr(myname_,'MAPL_GenericFinalize(), rc =',ier)
!		  if(.not.present(rc)) call die(myname_)
!		  rc=ier
!		  call timer_fnl(myname)
!		  _EXIT_(myname_)
!		  _MPCLOSE_()
!		  return
!		endif

  call ESMF_GridCompDestroy(pertmod_gc, rc=ier)
#ifdef _EXITONERROR_
  		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_GridCompDestroy(pertmod_gc), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
#endif /* _EXITONERROR_ */

  call ESMF_StateDestroy(pertmod_im, rc=ier)
#ifdef _EXITONERROR_
  		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_StateDestroy(pertmod_im), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
#endif /* _EXITONERROR_ */

  call ESMF_StateDestroy(pertmod_ex, rc=ier)
#ifdef _EXITONERROR_
  		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_StateDestroy(pertmod_ex), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
#endif /* _EXITONERROR_ */

  call myClock_destroy_(pertmod_clock, rc=ier)
  		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'myClock_destroy_(pertmod_clock), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif

  call myClock_destroy_(pertmod_clock_forward_, rc=ier)
  		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'myClock_destroy_(pertmod_clock_forward_), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
  deallocate(pertmod_clock_forward_)

  if(allocated(pertmod_clock_forward_)) then
  	! This checking of allocated() is necessary, since the code may never
	! entered the AD mode.
    call myClock_destroy_(pertmod_clock_reverse_, rc=ier)
  		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'myClock_destroy_(pertmod_clock_reverse_), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif
    deallocate(pertmod_clock_reverse_)
  endif

  call ESMF_ConfigDestroy(pertmod_config, rc=ier)
   		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,"ESMF_ConfigDestroy(pertmod_config), rc =",ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  call timer_fnl(myname)
		  _EXIT_(myname_)
		  _MPCLOSE_()
		  return
		endif

! release ak/bk from pert model
  deallocate(pert_ak,pert_bk)

  call timer_fnl(myname)
_EXIT_(myname_)
_MPCLOSE_()
end subroutine pertmod_finalize

subroutine myClock_create_(myClock)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 subroutine myClock_create_
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:	 2010-08-13
!
! abstract: - create an ESMF_Clock
!
! program history log:
!   2010-08-13  j guo   - added this document block
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

! function interface:

  use ESMF, only: ESMF_SUCCESS
  use ESMF, only: ESMF_Clock
  use ESMF, only: ESMF_Time
  use ESMF, only: ESMF_TimeInterval
  use ESMF, only: ESMF_ClockGet
  use ESMF, only: ESMF_TimeIntervalSet
  use ESMF, only: ESMF_ClockCreate
  use ESMF, only: ESMF_ClockSet
  use ESMF, only: ESMF_DIRECTION_FORWARD
  use ESMF, only: ESMF_DIRECTION_REVERSE
  use ESMF, only: ESMF_ClockAdvance
  use ESMF, only: ESMF_Direction_Flag
  use MAPL, only: MAPL_am_I_root

  use kinds    , only: i_kind
  use mpeu_util, only: tell,warn,perr,die
  use mpeu_util, only: stdout_open
  use mpeu_util, only: stdout_close
  implicit none
  type(ESMF_clock)    ,intent(out):: myClock

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname_=myname//'::myClock_create_'
  type(ESMF_Time)            :: currTime, startTime, stopTime
  type(ESMF_TimeInterval)    :: timeStep
  integer(kind=i_kind):: rc
_ENTRY_(myname_)

  rc =ESMF_SUCCESS

  call myConfig_getTime_(pertmod_config,startTime,'BEG_DATE:')
  call myConfig_getTime_(pertmod_config, stopTime,'END_DATE:')
  call myConfig_getTimeInterval_(pertmod_config,timeStep,'RUN_DT:')

		if(allocated(pertmod_clock_forward_)) &
			call die(myname_,'already created_, pertmod_clock_forward_')

  	! Creating a clock_forward_ from scratch
  allocate(pertmod_clock_forward_)
  pertmod_clock_forward_ = ESMF_ClockCreate(name=myname//"_clock_forward_", &
	timeStep=timeStep, startTime=startTime, stopTime=stopTime, rc=rc)

		if(rc/=ESMF_SUCCESS) call die(myname_, &
			'pertmod_clock_forward_ = ESMF_ClockCreate(), rc =',rc)

  call ESMF_ClockSet(pertmod_clock_forward_, currTime=startTime, rc=rc)
		if(rc/=ESMF_SUCCESS) call die(myname_,'ESMF_ClockSet(rurrTime), rc =',rc)

  call ESMF_ClockSet(pertmod_clock_forward_, direction=ESMF_DIRECTION_FORWARD, rc=rc)
		if(rc/=ESMF_SUCCESS) call die(myname_,'ESMF_ClockSet(ESMF_DIRECTION_FORWARD), rc =',rc)

  call myClock_show_(pertmod_clock_forward_,myname_,'initial-pertmod_clock',show=.true.)

  myClock = ESMF_ClockCreate(pertmod_clock_forward_, rc=rc)
	  	if(rc/=ESMF_SUCCESS) call perr(myname_, &
			'myClock = ESMF_ClockCreate(pertmod_clock_forward_), rc =',rc)

_EXIT_(myname_)
end subroutine myClock_create_

subroutine myClock_show_(aClock,who,what,show)
  use ESMF, only: ESMF_SUCCESS
  use ESMF, only: ESMF_Clock
  use ESMF, only: ESMF_ClockGet
  use ESMF, only: ESMF_Time
  use ESMF, only: ESMF_TimeGet
  use ESMF, only: ESMF_TimeInterval
  use ESMF, only: ESMF_TimeIntervalGet
  use ESMF, only: ESMF_Direction_Flag
  use ESMF, only: ESMF_DIRECTION_FORWARD
  use ESMF, only: ESMF_DIRECTION_REVERSE
  use ESMF, only: ESMF_CalKind_Flag
  use ESMF, only: ESMF_CALKind_GREGORIAN
  use ESMF, only: operator(==)
  use MAPL, only: MAPL_AM_I_ROOT
  use mpeu_util, only: tell,warn,perr,die
  implicit none
  type(ESMF_Clock),intent(in):: aClock
  character(len=*),intent(in):: who
  character(len=*),intent(in):: what
  logical         ,intent(in):: show

  character(len=*),parameter:: myname_=myname//"::myClock_show_"
  character(len=20):: timestr
  type(ESMF_Time):: bTime,cTime,eTime
  type(ESMF_TimeInterval):: dTime
  type(ESMF_Direction_Flag):: direction
  type(ESMF_CalKind_Flag):: calKindFlag
  integer(i_kind):: d,h,m,s
  logical:: show_
  integer(i_kind):: ier

  if(.not.MAPL_AM_I_ROOT()) return
  	! Note: this "if(I_am_not_root()) return" should only apply to
	! local algorithms.  Otherwise, it may become a dead-lock.

  show_=show
#ifdef DEBUG_CLOCK
  show_=.true.
#endif
  if(.not.show_) return

  call ESMF_ClockGet(aClock, startTime=btime, stopTime=etime, &
  			      currTime=ctime, timeStep=dTime, &
		          calKindFlag=calKindFlag, direction=direction, rc=ier )
		if(ier/=ESMF_SUCCESS) &
			call die(myname_,'ESMF_ClockGet(), rc =',ier)

  call ESMF_TimeGet(btime,timeStringISOFrac=timestr)
  call tell(who,'startTime('//trim(what)//') = '//trim(timestr))
  call ESMF_TimeGet(ctime,timeStringISOFrac=timestr)
  call tell(who,' currTime('//trim(what)//') = '//trim(timestr))
  call ESMF_TimeGet(etime,timeStringISOFrac=timestr)
  call tell(who,' stopTime('//trim(what)//') = '//trim(timestr))

  call ESMF_TimeIntervalGet(dtime,d=d,h=h,m=m,s=s)
  call tell(who,' timeStep('//trim(what)//') = ', &
  	(/((d*100+h)*100+m)*100+s/),format='(i10.8)')

  if(direction == ESMF_DIRECTION_FORWARD) then
    call tell(who,'direction('//trim(what)//') = ESMF_DIRECTION_FORWARD')
  elseif(direction == ESMF_DIRECTION_REVERSE) then
    call tell(who,'direction('//trim(what)//') = ESMF_DIRECTION_REVERSE')
  else
    call tell(who,'direction('//trim(what)//') = .unknown.')
  endif

!  if(calKindFlag == ESMF_CALKIND_GREGORIAN) then
!    call tell(who,'  calType('//trim(what)//') is ESMF_CALKIND_GREGORIAN')
!  else
!    call tell(who,'  calType('//trim(what)//') is .NOT. ESMF_CALKIND_GREGORIAN')
!  endif

end subroutine myClock_show_

subroutine myClock_TLconfig_(myClock,iymd,ihms,tstep,rc)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 subroutine myClock_TLconfig_
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:	 2011-02-10
!
! abstract: 
!
! program history log:
!   2011-02-10  j guo   - added this document block
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

  use ESMF, only: ESMF_SUCCESS
  use ESMF, only: ESMF_Clock
  use ESMF, only: ESMF_ClockCreate
  use ESMF, only: ESMF_ClockDestroy
  use ESMF, only: ESMF_ClockGet
  use ESMF, only: ESMF_ClockSet
  use ESMF, only: ESMF_DIRECTION_FORWARD
  use ESMF, only: ESMF_Time
  use ESMF, only: ESMF_TimeGet
  use ESMF, only: ESMF_TimeSet
  use ESMF, only: ESMF_TimeInterval
  use ESMF, only: ESMF_TimeIntervalGet
  use ESMF, only: ESMF_TimeIntervalSet
  use ESMF, only: ESMF_Direction_Flag
  use ESMF, only: operator(/=)
  use kinds    , only: i_kind
  use mpeu_util, only: tell,warn,perr,die
  use mpeu_util, only: stdout_open
  use mpeu_util, only: stdout_close
  implicit none
  type(ESMF_clock)    ,intent(inout):: myClock
  integer(kind=i_kind),intent(in ):: iymd	! initial date (YYYYMMDD) of the perturbation state
  integer(kind=i_kind),intent(in ):: ihms	! initial time (HHMMSS) of the perturbation state
  integer(kind=i_kind),intent(out):: tstep	! a single TLM time step in seconds
  integer(kind=i_kind),intent(out):: rc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname_=myname//'::myClock_TLconfig_'
  type(ESMF_TimeInterval):: timeStep
  type(ESMF_Time):: startTime
  type(ESMF_Time):: stopTime
  type(ESMF_Direction_Flag):: direction
_ENTRY_(myname_)

  rc=ESMF_SUCCESS
  tstep = huge(tstep)	! must be defined

  	! Destroy the current clock to get a fresh one.
  call ESMF_ClockDestroy(myClock, rc=rc)
  		if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_ClockDestroy(), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

  	! make a copy of pertmod_clock_forward_ for this TL-AD cycle.
  myClock=ESMF_ClockCreate(pertmod_clock_forward_, rc=rc)
		if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_ClockCreate(), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

  	! Get tstep from myClock
    call ESMF_ClockGet(myClock, timeStep= timeStep, &
  			       direction=direction, rc=rc)
		if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_ClockGet(), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call ESMF_TimeIntervalGet (timeStep, S= tstep, rc=rc)
		if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_TimeIntervalGet(), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

  call myClock_show_(myClock,myname_,'pertmod_clock',show=.false.)

_EXIT_(myname_)
end subroutine myClock_TLconfig_
subroutine myClock_ADconfig_(myClock,iymd,ihms,tstep,rc)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:	 subroutine myClock_ADconfig_
!   prgmmr:	 j guo <jguo@nasa.gov>
!      org:	 NASA/GSFC, Global Modeling and Assimilation Office, 900.3
!     date:	 2011-02-10
!
! abstract: 
!
! program history log:
!   2011-02-10  j guo   - added this document block
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

  use ESMF , only: ESMF_SUCCESS
  use ESMF , only: ESMF_Clock
  use ESMF , only: ESMF_ClockGet
  use ESMF , only: ESMF_ClockSet
  use ESMF , only: ESMF_ClockCreate
  use ESMF , only: ESMF_ClockDestroy
  use ESMF , only: ESMF_DIRECTION_REVERSE
  use ESMF , only: ESMF_Time
  use ESMF , only: ESMF_TimeGet
  use ESMF , only: ESMF_TimeSet
  use ESMF , only: ESMF_TimeInterval
  use ESMF , only: ESMF_TimeIntervalGet
  use ESMF , only: ESMF_TimeIntervalSet
  use ESMF , only: ESMF_Direction_Flag
  use ESMF , only: operator(/=)

  use kinds    , only: i_kind
  use mpeu_util, only: tell,warn,perr,die
  use mpeu_util, only: stdout_open
  use mpeu_util, only: stdout_close
  implicit none
  type(ESMF_clock)    ,intent(inout):: myClock
  integer(kind=i_kind),intent(in ):: iymd	! initial date (YYYYMMDD) of the perturbation state
  integer(kind=i_kind),intent(in ):: ihms	! initial time (HHMMSS) of the perturbation state
  integer(kind=i_kind),intent(out):: tstep	! a single ADM time step in seconds
  integer(kind=i_kind),intent(out):: rc

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  character(len=*),parameter :: myname_=myname//'::myClock_ADconfig_'
  type(ESMF_TimeInterval):: timeStep
  type(ESMF_Time):: startTime
  type(ESMF_Time):: stopTime
  type(ESMF_Direction_Flag):: direction

_ENTRY_(myname_)

  rc=ESMF_SUCCESS
  tstep = huge(tstep)	! default return value is a Not-A-TimeStep

  if(.not.allocated(pertmod_clock_reverse_)) then
    allocate(pertmod_clock_reverse_)
    pertmod_clock_reverse_ = ESMF_ClockCreate(pertmod_clock, rc=rc)
		if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_ClockCreate(1), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

    call ESMF_ClockSet(pertmod_clock_reverse_, &
    	direction = ESMF_DIRECTION_REVERSE, rc=rc)
		if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_ClockSet(REVERSE), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif
  endif

  	! Destroy the current clock to get a fresh one.
  call ESMF_ClockDestroy(myClock, rc=rc)
  		if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_ClockDestroy(), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

  pertmod_clock = ESMF_ClockCreate(pertmod_clock_reverse_, rc=rc)
		if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_ClockCreate(2), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

  	! Get tstep from myClock
  call ESMF_ClockGet(myClock, timeStep = timeStep, &
    				direction=direction, rc=rc)
		if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_ClockGet(timeStep, direction), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

  call ESMF_TimeIntervalGet (timeStep, S= tstep, rc=rc)
		if(rc/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_TimeIntervalGet(), rc =',rc)
		  _EXIT_(myname_)
		  return
		endif

  call myClock_show_(myClock,myname_,'pertmod_clock',show=.false.)

_EXIT_(myname_)
end subroutine myClock_ADconfig_

subroutine myconfig_getTime_(myConfig,atime,label,rc)
  use ESMF, only: ESMF_Config
  use ESMF, only: ESMF_ConfigGetAttribute
  use ESMF, only: ESMF_Time
  use ESMF, only: ESMF_TimeSet
  use ESMF, only: ESMF_SUCCESS
  use mpeu_util, only: tell,perr,die
  implicit none
  type(ESMF_Config)       ,intent(inout):: myConfig
  type(ESMF_Time  )       ,intent(  out):: atime
  character(len=* )       ,intent(in   ):: label
  integer(i_kind),optional,intent(  out):: rc

  integer(i_kind),dimension(2):: datetime
  integer(i_kind):: ier
  integer(i_kind):: yr,mo,dy,hr,mi,sc

  character(len=*),parameter:: myname_=myname//'::myconfig_getTime_'
#ifdef _MYNAME_
#undef _MYNAME_
#endif
#define _MYNAME_ myname_//'/'//trim(label)
_ENTRY_(_MYNAME_)
  if(present(rc)) rc=ESMF_SUCCESS

  datetime=-1
  call ESMF_ConfigGetAttribute(myConfig,datetime, count=2, label=label, rc=ier)
  		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_ConfigGetAttribute("'//trim(label)//'"), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  _EXIT_(_MYNAME_)
		  return
		endif

	  	if(datetime(1)<=0.or.datetime(2)<0) then
		  call perr(myname_,'bad value of "'//trim(label)// &
		  	'", yyyymmdd:hhmmss =',datetime,format='(i10.6)')
		  if(.not.present(rc)) call die(myname_)
		  rc=GENERIC_ERROR
		  _EXIT_(_MYNAME_)
		  return
		endif

  yr=datetime(1)
  dy=mod(yr,100)
  yr=    yr/100
  mo=mod(yr,100)
  yr=    yr/100

  hr=datetime(2)
  sc=mod(hr,100)
  hr=    hr/100
  mi=mod(hr,100)
  hr=    hr/100

  call ESMF_TimeSet(atime,yy=yr,mm=mo,dd=dy,h=hr,m=mi,s=sc,rc=ier)
  		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_TimeSet("'//trim(label)//'"), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  _EXIT_(_MYNAME_)
		  return
		endif

#ifdef DEBUG_VERBOSE
  call tell(myname_,'value of "'//trim(label)// &
	  	'", yyyymmdd:hhmmss =',datetime,format='(i10.6)')
#endif
_EXIT_(_MYNAME_)
end subroutine myconfig_getTime_

subroutine myconfig_getTimeInterval_(myConfig,aStep,label,rc)
  use ESMF, only: ESMF_Config
  use ESMF, only: ESMF_ConfigGetAttribute
  use ESMF, only: ESMF_TimeInterval
  use ESMF, only: ESMF_TimeIntervalSet
  use ESMF, only: ESMF_SUCCESS
  use mpeu_util, only: tell,perr,die
  implicit none
  type(ESMF_Config)       ,intent(inout):: myConfig
  type(ESMF_TimeInterval) ,intent(  out):: aStep
  character(len=*)        ,intent(in   ):: label
  integer(i_kind),optional,intent(  out):: rc

  integer(i_kind),dimension(2):: datetime
  integer(i_kind):: ier
  integer(i_kind):: dtsec

  character(len=*),parameter:: myname_=myname//'::myconfig_getTimeInterval_'
#ifdef _MYNAME_
#undef _MYNAME_
#endif
#define _MYNAME_ myname_//'/'//trim(label)
_ENTRY_(_MYNAME_)
  if(present(rc)) rc=ESMF_SUCCESS

  dtsec=0
  call ESMF_ConfigGetAttribute(myConfig,dtsec,label=label, rc=ier)
  		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_ConfigGetAttribute("'//trim(label)//'"), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  _EXIT_(_MYNAME_)
		  return
		endif

  call ESMF_TimeIntervalSet(aStep,s=dtsec, rc=ier)
  		if(ier/=ESMF_SUCCESS) then
		  call perr(myname_,'ESMF_TimeIntervalSet("'//trim(label)//'"), rc =',ier)
		  if(.not.present(rc)) call die(myname_)
		  rc=ier
		  _EXIT_(_MYNAME_)
		  return
		endif

#ifdef DEBUG_VERBOSE
  call tell(myname_,'"'//trim(label)//'" in seconds =',dtsec)
#endif
_EXIT_(_MYNAME_)
end subroutine myconfig_getTimeInterval_

end module geos_pertmod
