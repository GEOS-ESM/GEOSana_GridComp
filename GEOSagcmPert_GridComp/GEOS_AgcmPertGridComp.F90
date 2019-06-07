! $Id$
!#define DOT_RAND_PERT
!#define DO_DOT_TEST
!#define VERBOSE
!#define DEBUG_METACOMP

!#define DO_ADTEST_
!#define NO_TIME_LOOP_
!#define NO_PHYSCOMP_

#define R8  8

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: 

! GEOS\_AgcmPertGridCompMod -- A Module to combines the linear Supedynamics 
!     and Physics Gridded Components used by GSI in the inner loop of the
!     4D-Var minimization.

! !INTERFACE:

module GEOS_AgcmPertGridCompMod

! !USES:

  use ESMF
  use MAPL_Mod
  use GEOS_PertSharedMod

  use GEOS_DynCorePertGridCompMod,  only:  DYNA_SetServices => SetServices
  use GEOS_PhysicsPertGridCompMod,  only:  PHYS_SetServices => SetServices

  use mpp_domains_mod,  only: DGRID_NE, DGRID_SW
  use mpp_domains_mod,  only: mpp_get_boundary

  use GEOS_DynCorePertGridCompMod, only: fv_atm

  use mpp_mod, only: mpp_pe
  use fv_timing_mod,  only: timing_on, timing_off, timing_prt

#ifdef DO_ADTEST_
  use GEOS_AgcmPertADTester, only: ADTester_test
#endif

  use MAPL_AbstractRegridderMod
  use MAPL_RegridderManagerMod
  use MAPL_TransposeRegridderMod

  use tapenade_iter, only: cp_iter_controls, initialize_cp_iter, finalize_cp_iter

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION:  The component has two children: DYNAmics and PHYSics.
!     Its Run method has 2 phases, corresponding to TLM and Adjoint
!     propagations of the trajectory and perturbation states. 
!     Initialization can be done beginning with either a TLM or
!     an adjoint propagation. 

!     The component's natural grid is the grid of the analysis, from
!     which we assume the operators are being called. Currently, this
!     grid is a Lat-Lon grid and the TLM and Adjoint operate on a
!     Cubed-Sphere grid.

!EOP

  integer,save :: DYNA
  integer,save :: PHYS

  real, parameter :: EPS = MAPL_VIREPS

! Private internal state type for this GridComp
!----------------------------------------------

  type T_AGCMPert_STATE

! These two lists of ESMF_FieldBundles are the "master" copies of snapshots
!  of the trajectory at the end points of each propagation interval (TRAJ(:)) 
!  and of a possible trajectory tendency (DTRAJ(:)) within each interval.

     type (ESMF_FIELDBUNDLE), allocatable :: TRAJ3(:)
     type (ESMF_FIELDBUNDLE), allocatable :: TRAJ2(:)

     type(T_3DVar),  allocatable, dimension(:) :: TRAJwrk3, PERTwrk 
     type(T_2DVar),  allocatable, dimension(:) :: TRAJwrk2 

     type(T_3DVar),  allocatable, dimension(:) :: PERTdel

! The current interval. Each call to AGCM does one full interval.

     integer                              :: CallCounter = 0

! Model time steps per call

     integer                              :: StepsPerCall

! The GCM time step in seconds (== length of interval / StepsPerCall)

     real                                 :: DT

! Handles to the various transforms, which are set in Initialize

     class (AbstractRegridder), pointer :: L2C, C2L
     type (TransposeRegridder) :: L2Ct, C2Lt
     type (ESMF_Clock)                    :: CLOCK
  end type T_AGCMPERT_STATE

! Wrapper for carrying the private internal state in the GC
! ---------------------------------------------------------

  type T_AGCMPert_Wrap
     type (T_AGCMPert_State), pointer :: State
  end type T_AGCMPERT_WRAP

  type(ESMF_VM), save :: VM

  interface echo_time_
         module procedure echo_time0_
         module procedure echo_time1_
  end interface

  real(ESMF_KIND_R8),allocatable,dimension(:,:,:,:),save :: pertWrk_saved_

  integer, save                                   :: MEMTRAJ
  type (ESMF_Time),allocatable, save              :: Time_lst(:)
  character(len=ESMF_MAXSTR), allocatable         :: TRAJ_FILE_TMPL(:)
  integer, allocatable, save                      :: Traj_Orgi_Idx(:), Traj_Skip_Idx(:)

  type(ESMF_FieldBundle)                          :: TRAJwrk3_4gq, TRAJwrk2_4gq
  integer, save                                   :: gq_ntrj
  real, save                                      :: gq_wght

  logical :: verbose

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp)                :: GC  ! gridded component
    integer            , intent(  out) :: RC  ! return code

! !DESCRIPTION:  The SetServices for the Physics GC needs to register its
!   Initialize and Run.  It uses the MAPL_Generic construct for defining 
!   state specs and couplings among its children.  In addition, it creates the   
!   children GCs (SURF, CHEM, RADIATION, MOIST, TURBULENCE) and runs their
!   respective SetServices.

!EOP

    type (MAPL_MetaComp),      pointer  :: MAPL

! ErrLog Variables
! ----------------

    character(len=ESMF_MAXSTR)         :: IAm
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME
    integer                            :: phase0, phase1

! Begin
! -----

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

!BOP

! !IMPORT STATE:
    
    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'U'  ,                                       &
         LONG_NAME  = 'eastward_wind_analysis_increment',          &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'V'  ,                                       &
         LONG_NAME  = 'northward_wind_analysis_increment',         &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'TV'  ,                                      &
         LONG_NAME  = 'virtual_temperature_analysis_increment',    &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'DP',                                        &
         LONG_NAME  = 'pressure_thickness_analysis_increment',     &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'QV'  ,                                      &
         LONG_NAME  = 'specific_humidity_analysis_increment',      &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'QL'  ,                                      &
         LONG_NAME  = 'cloud_liquid_water_analysis_increment',     &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'QI'  ,                                      &
         LONG_NAME  = 'cloud_ice_analysis_increment',              &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'O3'  ,                                      &
         LONG_NAME  = 'ozone_analysis_increment',                  &
         UNITS      = 'ppmv',                                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

! !EXPORT STATE:

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'perturbation_eastward_wind_on_A-grid',      &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V',                                         &
         LONG_NAME  = 'perturbation_northward_wind_on_A-grid',     &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'TV',                                        &
         LONG_NAME  = 'perturbation_air_virtual_temperature',      &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DP',                                        &
         LONG_NAME  = 'perturbation_air_pressure_thickness',       &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'QV',                                        &
         LONG_NAME  = 'perturbation_specific_humidity',            &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'QL',                                        &
         LONG_NAME  = 'perturbation_cloud_liquid_water_mixing_ration',&
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'QI',                                        &
         LONG_NAME  = 'perturbation_cloud_ice_mixing_ration',      &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'O3',                                        &
         LONG_NAME  = 'perturbation_ozone_mole_mixing_ration',     &
         UNITS      = 'ppmv',                                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

!EOP

! Get the component switches
! --------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

    !Get the ObserverMode switch
    call MAPL_GetResource( MAPL, ObserverMode, 'PertObserverMode:', default = 0, rc = STATUS )
    VERIFY_(STATUS)

    !Dynamics switches
    call MAPL_GetResource( MAPL, DO_DYNAMICS, 'DO_DYNAMICS:', default = 1, rc = STATUS )
    VERIFY_(STATUS)

    !Boundary layer switch
    call MAPL_GetResource( MAPL, DO_BL_PHYS, 'DO_TURBULENCE:', default = 1, rc = STATUS )
    VERIFY_(STATUS)

    !Moist physics switch
    call MAPL_GetResource( MAPL, DO_MOIST_PHYS, 'DO_MOIST:', default = 2, rc = STATUS )
    VERIFY_(STATUS)

    !Radiation switch
    call MAPL_GetResource( MAPL, DO_RAD_PHYS, 'DO_RADIATION:', default = 0, rc = STATUS )
    VERIFY_(STATUS)

    !Gravity wave drag switch
    call MAPL_GetResource( MAPL, DO_GWD_PHYS, 'DO_GWAVEDRAG:', default = 0, rc = STATUS )
    VERIFY_(STATUS)

    !Interpolate winds D2A switch
    call MAPL_GetResource( MAPL, DO_PHYS_D2A, 'DO_PHYS_D2A:', default = 0, rc = STATUS )
    VERIFY_(STATUS)

! Switches to control the custom Tapenade checkpointing interface
! ---------------------------------------------------------------

    call MAPL_GetResource( MAPL, cp_iter_controls%cp_i , 'CP_i:'          , default = 0  , RC = STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, cp_iter_controls%cp_nt, 'CallsPerWindow:', default = 0  , RC = STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, cp_iter_controls%cp_gb, 'CP_gb:'         , default = 0.0, RC = STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, cp_iter_controls%cp_nm, 'CP_nm:'         , default = 0  , RC = STATUS )
    VERIFY_(STATUS)
    cp_iter_controls%cp_t = 0

    !If iterative then allocate variable that will hold everything
    call initialize_cp_iter

! Register services for this component
! ------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,  Initialize,  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE, Finalize,    RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,   Run,         RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,   Run,         RC=STATUS )
    VERIFY_(STATUS)

    phase0=phase(import=.true. )
    phase1=phase(import=.false.)

    if(phase0/=phase1) then

      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,   Run,         RC=STATUS )
      VERIFY_(STATUS)
      call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,   Run,         RC=STATUS )
      VERIFY_(STATUS)

      phase0=phase(export=.true. )
      phase1=phase(export=.false.)

      if(phase0/=phase1) then

        call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,   Run,         RC=STATUS )
        VERIFY_(STATUS)
        call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,   Run,         RC=STATUS )
        VERIFY_(STATUS)

        call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,   Run,         RC=STATUS )
        VERIFY_(STATUS)
        call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,   Run,         RC=STATUS )
        VERIFY_(STATUS)

	if(MAPL_AM_I_ROOT()) write(*,'(1x,2a)') trim(Iam),' -- 8 RUN() phases set'

      else
	if(MAPL_AM_I_ROOT()) write(*,'(1x,2a)') trim(Iam),' -- 4 RUN() phases set'
      endif

    else
	if(MAPL_AM_I_ROOT()) write(*,'(1x,2a)') trim(Iam),' -- 2 RUN() phases set'
    endif

! Create childrens gridded components and invoke their SetServices
! ----------------------------------------------------------------

    DYNA = MAPL_AddChild(GC, NAME='DYNAMICS_PERT', SS=DYNA_SetServices, RC=STATUS)
    VERIFY_(STATUS)
#ifdef DEBUG_METACOMP
    call MAPL_MetaCompPrint(GC,header=trim(Iam)//'(), after mapl_addChild(DYNA)')
#endif

#ifndef NO_PHYSCOMP_
    PHYS = MAPL_AddChild(GC, NAME='PHYSICS_PERT',  SS=PHYS_SetServices, RC=STATUS)
    VERIFY_(STATUS)
#ifdef DEBUG_METACOMP
    call MAPL_MetaCompPrint(GC,header=trim(Iam)//'(), after mapl_addChild(PHYS)')
#endif
#endif

! SetServices for children
!-------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS )
    VERIFY_(STATUS)
#ifdef DEBUG_METACOMP
    call MAPL_MetaCompPrint(GC,header=trim(Iam)//'(), after MAPL_GenericSetServices()')
#endif

! Clocks
!-------

    call MAPL_TimerAdd(GC, name="INITIALIZE" , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="--MakeGrid" , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="--ReadTraj" , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="--Transform", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="RUN"        , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="--RUNTLM"   , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="----RunDynaTLM", RC=STATUS)
    VERIFY_(STATUS)
#ifndef NO_PHYSCOMP_
    call MAPL_TimerAdd(GC, name="----RunPhysTLM", RC=STATUS)
    VERIFY_(STATUS)
#endif
    call MAPL_TimerAdd(GC, name="--RUNADJ"   , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="----RunDynaADJ", RC=STATUS)
    VERIFY_(STATUS)
#ifndef NO_PHYSCOMP_
    call MAPL_TimerAdd(GC, name="----RunPhysADJ", RC=STATUS)
    VERIFY_(STATUS)
#endif
    call MAPL_TimerAdd(GC, name="FINALIZE"   , RC=STATUS)
    VERIFY_(STATUS)


! All done
!---------

#ifdef DEBUG_METACOMP
    call MAPL_MetaCompPrint(GC,header=trim(Iam)//'(), exiting ..',deep=.true.)
#endif
    RETURN_(ESMF_SUCCESS)  
  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: Initialize -- Initialize method for the composite Agcm Gridded Component

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer         ,    intent(  out) :: RC     ! Error code

! !DESCRIPTION: 
 

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm 
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (T_AGCMPert_Wrap)              :: WRAP
    type (T_AGCMPert_State),   pointer  :: STATE
    type (MAPL_MetaComp),      pointer  :: MAPL

    type(ESMF_Config)                   :: CF

    type (ESMF_State),         pointer  :: GIM(:)
    type (ESMF_GridComp),      pointer  :: GCS(:)
    type (ESMF_FieldBundle)             :: Bundle
    type (ESMF_TimeInterval)            :: DT
    type (ESMF_Time)                    :: Time
    type (ESMF_Grid)                    :: CUGrid
    type (ESMF_Grid)                    :: LLGrid
    character(len=ESMF_MAXSTR)          :: TRAJ_FILE
    integer                             :: I, J, N, Interval, CallsPerWindow
    type(ESMF_FieldBundle)              :: TRAJwrk3, TRAJwrk2, PERTwrk
    type(ESMF_FieldBundle)              :: TRAJwrk3_tmp, TRAJwrk2_tmp

    type(ESMF_FieldBundle)              :: PERTdel
    type(ESMF_FieldBundle)              :: LLPert
    type(ESMF_Time)                     :: startTime, stopTime
    type(ESMF_TimeInterval)             :: timeStep
    real, pointer                       :: PT(:,:,:)

    integer       :: seed
    integer       :: L, lindx
    real, pointer :: tvi(:,:,:)
    real, pointer :: lons(:,:)
    real, pointer :: lats(:,:)
    integer       :: PERT_PHASE
    integer       :: N_TRAJ, N_TRAJ_SKIP
    integer       :: FREQ_YY, FREQ_MM, FREQ_DD, FREQ_H, FREQ_M, FREQ_S
    integer       :: TRAJ_FREQ
    real          :: TRAJ_FREQ_DT
    integer       :: DO_ESMF_INTERP  ! temporary flag to allow for ease of testing

    character(len=ESMF_MAXSTR) :: PERT_FILE, PERT_FILE_DEF, TILEFILE, lstvars


! =============================================================================

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, config=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Initialize"
#ifdef DEBUG_METACOMP
    call MAPL_MetaCompPrint(GC,header=trim(Iam)//'(), entered ..')
#endif

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)
#ifdef DEBUG_METACOMP
    call MAPL_MetaCompPrint(GC,header=trim(Iam)//'(), after MAPL_GetObjectyFromGC()')
#endif

    call MAPL_TimerOn(MAPL,"INITIALIZE")
    call MAPL_TimerOn(MAPL,"TOTAL")

! Get my children GCs and IMPORT States
!--------------------------------------

    call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, RC=STATUS )
    VERIFY_(STATUS)
#ifdef DEBUG_METACOMP
    call MAPL_MetaCompPrint(GC,header=trim(Iam)//'(), after MAPL_Get()')
#endif

! Create the children's cube-sphere grids
!----------------------------------------

    call MAPL_TimerOn(MAPL,"--MakeGrid")

    call MAPL_GridCreate(GCS(DYNA), rc=status)
    VERIFY_(STATUS)
#ifndef NO_PHYSCOMP_
    call MAPL_GridCreate(GCS(PHYS), srcGC=GCS(DYNA), rc=status)
    VERIFY_(STATUS)
#endif

! Create my own grid
!-------------------

! Get the cubed grid from one of the children
!   and the lat-lon grid from self.
!--------------------------------------------

    call ESMF_GridCompGet ( GCS(DYNA), grid=CUGrid, RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_GridCompGet ( GC       , grid=LLGrid, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOff(MAPL,"--MakeGrid")

! Allocate this instance of the private internal state
!   and put it in wrapper.
! ----------------------------------------------------

    allocate( STATE, stat=status )
    VERIFY_(STATUS)

    wrap%state => STATE
 
! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC,'PERTstate', wrap, STATUS )
    VERIFY_(STATUS)

! We are cloning the global clock by re-creating it
! leaving the alarms intentionally (not to interfere with History)
!----------------------------------------------------------------
    call echo_time_(Clock,'inside init before extract curr time')
    call ESMF_ClockGet(Clock, CurrTime=Time, &
         timeStep=timeStep, StartTime=StartTime, stopTime=stopTime, RC=STATUS)
    VERIFY_(STATUS)
    call echo_time_(Time,'inside init after extract curr time')

    State%Clock = ESMF_ClockCreate(name=trim(COMP_NAME)//"Clock", &
         timeStep=timeStep, startTime=startTime, stopTime=stopTime, RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_ClockSet(State%Clock, CurrTime=Time, rc=status )
    VERIFY_(STATUS)


! The name of the file containing the trajectory and 
!  parameters describing the analysis window. Defaults
!  assume 30 minute model time step and 12 hour window.
!------------------------------------------------------

    call gq_(CF, gq_ntrj, gq_wght)
    allocate(TRAJ_FILE_TMPL(gq_ntrj))
    call get_TRAJ_FILE_TMPL_(CF, gq_ntrj, TRAJ_FILE_TMPL)

    call MAPL_GetResource ( MAPL, STATE%DT,            'RUN_DT:',                       __RC__ )
    call MAPL_GetResource ( MAPL, CallsPerWindow,      'CallsPerWindow:',   default=24, __RC__ )
    call MAPL_GetResource ( MAPL, STATE%StepsPerCall,  'StepsPerCall:',     default=1,  __RC__ )
    call MAPL_GetResource ( MAPL, TILEFILE,            'TILEFILE:',         default='', __RC__ ) !JK Neutral TILE file
    call MAPL_GetResource ( MAPL, MEMTRAJ,             'MEMTRAJ:',          default=1,  rc=STATUS )  
    VERIFY_(STATUS)                                           
    call MAPL_GetResource ( MAPL, TRAJ_FREQ,           'TRAJ_FREQ:',        default=0,   rc=STATUS )
    VERIFY_(STATUS)                                           
    call MAPL_GetResource ( MAPL, DO_ESMF_INTERP,      'DO_ESMF_INTERP:',   default=0,   rc=STATUS )
    VERIFY_(STATUS)                                           


! Switches for test cases
! -----------------------

    !Switch for baroclinic wave test
    call MAPL_GetResource( MAPL, DynTestCase, 'DynamicsTestCase:', default=0, rc = STATUS )
    VERIFY_(STATUS)

    !Switch and constants for tracer advection test
    if (DynTestCase == 3 .or. DynTestCase == 4) then
       call MAPL_GetResource( MAPL, TrT_ver, 'TrTest_ver:', default=0, rc = STATUS )
       VERIFY_(STATUS)
       call MAPL_GetResource( MAPL, TrT_h0, 'TrTest_h0:', default=0.0, rc = STATUS )
       VERIFY_(STATUS)
       call MAPL_GetResource( MAPL, TrT_R, 'TrTest_R:', default=0.0, rc = STATUS )
       VERIFY_(STATUS)
       call MAPL_GetResource( MAPL, TrT_latc, 'TrTest_lat:', default=0.0, rc = STATUS )
       VERIFY_(STATUS)
       call MAPL_GetResource( MAPL, TrT_lonc, 'TrTest_lon:', default=0.0, rc = STATUS )
       VERIFY_(STATUS)
    endif

! For now we have this limitation, which we assert here
!------------------------------------------------------

    ASSERT_(State%StepsPerCall==1)

! Initialize the base state for the entire window by reading
!   the trajectory written during the outer loop segment into
!   the Traj list in the State. Each bundle in the list contains a
!   snapshot at the begining of each subinterval of the window,
!   of a state (U, V, PT, PE, and the tracers), on the
!   cubed sphere and at the inner loop resolution.
!-----------------------------------------------------------------

    call MAPL_TimerOn (MAPL,"--ReadTraj")

    If (DO_MOIST_PHYS.eq.0) Then
        NUM_GCM3Dvars = 8
        NUM_GCM2Dvars = 12
    Else
        NUM_GCM3Dvars = 11
        NUM_GCM2Dvars = 16
    Endif

    If (DO_RAD_PHYS .eq. 1) then
       NUM_GCM2Dvars = 24
    Endif
    

! Make bundles for the working trajectory and perturbation
!   These will be put in the children's IMPORT States. The
!   "del" bundles are only a tricky way of allocating the
!   state 3DVar copies. Note we create these whether we
!   need them or not. Another bit of sloppiness is to read
!   the TRAJwrk3 bundle. This is to load the PHIS.
!---------------------------------------------------------
 
    TRAJwrk3 = ESMF_FieldBundleCreate(name="TRAJWRK3",  __RC__)
    call ESMF_FieldBundleSet(TRAJwrk3, grid=CUGrid, __RC__ )
    TRAJwrk2 = ESMF_FieldBundleCreate(name="TRAJWRK2",  __RC__)
    call ESMF_FieldBundleSet(TRAJwrk2, grid=CUGrid, __RC__ )
    PERTwrk  = ESMF_FieldBundleCreate(name="PERTWRK",   __RC__)
    call ESMF_FieldBundleSet(PERTwrk, grid=CUGrid, __RC__ )
    PERTdel  = ESMF_FieldBundleCreate(name="PERTDEL",   __RC__)
    call ESMF_FieldBundleSet(PERTdel, grid=CUGrid, __RC__ )

    call Resolve_template_ ( TRAJ_FILE, TRAJ_FILE_TMPL(1), Time, trim(Iam), STATUS ); VERIFY_(STATUS)

    if (DO_MOIST_PHYS == 0) then
       allocate(STATE%TRAJwrk3(NUM_GCM3DvarsRead),stat=STATUS)
       VERIFY_(STATUS)
    else
        allocate(STATE%TRAJwrk3(NUM_GCM3DvarsReadMoist),stat=STATUS)
       VERIFY_(STATUS)
    endif
    allocate(STATE%PERTwrk(NUM_GCM3Dvars),stat=STATUS)
    VERIFY_(STATUS)
    allocate(STATE%TRAJwrk2(NUM_GCM2Dvars),stat=STATUS)
    VERIFY_(STATUS)
    allocate(STATE%PERTdel(NUM_GCM3Dvars),stat=STATUS)
    VERIFY_(STATUS)

    call MAPL_CFIORead(TRAJ_FILE, Time, PERTwrk,  noread=.true.,     __RC__)
    call MAPL_CFIORead(TRAJ_FILE, Time, PERTdel,  noread=.true.,     __RC__)

    if (DO_MOIST_PHYS == 0) then
       lstvars = vec2list(GCM3DvarsRead)
    else
       lstvars = vec2list(GCM3DvarsReadMoist)
    endif
    call MAPL_CFIORead(TRAJ_FILE, Time, TRAJwrk3, only_vars=lstvars, __RC__)
    lstvars = vec2list(GCM2Dvars(1:NUM_GCM2Dvars))
    call MAPL_CFIORead(TRAJ_FILE, Time, TRAJwrk2, only_vars=lstvars, __RC__)

    if ( gq_ntrj > 1 ) then

       TRAJwrk3_4gq = ESMF_FieldBundleCreate(name="TRAJWRK3_TMP",  __RC__)
       call ESMF_FieldBundleSet(TRAJwrk3_4gq, grid=CUGrid, __RC__ )
       TRAJwrk2_4gq = ESMF_FieldBundleCreate(name="TRAJWRK2_TMP",  __RC__)
       call ESMF_FieldBundleSet(TRAJwrk2_4gq, grid=CUGrid, __RC__ )

       call sax_Bundle_(TRAJwrk3, gq_wght)
       call sax_Bundle_(TRAJwrk2, gq_wght)

       do j=2,gq_ntrj

          call Resolve_template_ ( TRAJ_FILE, TRAJ_FILE_TMPL(J), Time, trim(Iam), STATUS ); VERIFY_(STATUS)

          if (DO_MOIST_PHYS == 0) then
             lstvars = vec2list(GCM3DvarsRead)
          else
             lstvars = vec2list(GCM3DvarsReadMoist)
          endif
          call MAPL_CFIORead(TRAJ_FILE, Time, TRAJwrk3_4gq, only_vars=lstvars, __RC__)
          lstvars = vec2list(GCM2Dvars(1:NUM_GCM2Dvars))
          call MAPL_CFIORead(TRAJ_FILE, Time, TRAJwrk2_4gq, only_vars=lstvars, __RC__)

          call saxpy_Bundle_(TRAJwrk3, TRAJwrk3_4gq, gq_wght)
          call saxpy_Bundle_(TRAJwrk2, TRAJwrk2_4gq, gq_wght)

       enddo

    endif

    if (DO_MOIST_PHYS == 0) then
       do i=1,NUM_GCM3DvarsRead
          State%TRAJwrk3(i)%name = GCM3DvarsRead(i)
       end do
    else
       do i=1,NUM_GCM3DvarsReadMoist
          State%TRAJwrk3(i)%name = GCM3DvarsReadMoist(i)
       end do
    endif
    do i=1,NUM_GCM3Dvars
       State%PERTwrk(i)%name  = GCM3Dvars(i)
       State%PERTdel(i)%name  = GCM3Dvars(i)
    end do

    do i=1,NUM_GCM2Dvars                                                                                 
       State%TRAJwrk2(i)%name = GCM2Dvars(i)
    end do

    call PERT_RefVarsFromBundle   (TrajWrk3, State%TRAJwrk3, rc=STATUS); VERIFY_(STATUS)
    call PERT_RefVarsFromBundle_2D(TrajWrk2, State%TRAJwrk2, rc=STATUS); VERIFY_(STATUS) 
    call PERT_RefVarsFromBundle   (PERTwrk,  State%PERTwrk,  rc=STATUS); VERIFY_(STATUS)
    call PERT_RefVarsFromBundle   (PERTdel,  State%PERTdel,  rc=STATUS); VERIFY_(STATUS)

! Put the working trajectory and perturbation Bundles
!   in the children's imports
!----------------------------------------------------

    do i=1,size(GIM)
       call MAPL_StateAdd(GIM(i), TRAJwrk3, rc=STATUS); VERIFY_(STATUS)
       call MAPL_StateAdd(GIM(i), PERTwrk,  rc=STATUS); VERIFY_(STATUS)
       call MAPL_StateAdd(GIM(i), TRAJwrk2, rc=STATUS); VERIFY_(STATUS)
    end do

! Create a bundle list containing the trajectory at the
!  end points of all intervals in the window.
!------------------------------------------------------
    call ESMF_TimeIntervalSet(DT, s=NINT(State%DT), rc=STATUS)
    VERIFY_(STATUS)

    allocate(TIME_lst(CallsPerWindow+1),stat=STATUS)
    VERIFY_(STATUS)

    if (TRAJ_FREQ.eq.0) then
       TRAJ_FREQ_DT = STATE%DT
    else
       TRAJ_FREQ_DT = real(TRAJ_FREQ)
    endif

    N_TRAJ_SKIP = int(TRAJ_FREQ_DT/STATE%DT)

    if (N_TRAJ_SKIP.le.0) then
        if(MAPL_AM_I_ROOT()) write(*,*) "The number of Traj skip is zero."
        call exit(0)
    endif

    call Get_Traj_Size(CallsPerWindow+1, N_TRAJ_SKIP, N_Traj)

    allocate(Traj_Orgi_Idx(CallsPerWindow+1),stat=STATUS)
    VERIFY_(STATUS)
    allocate(Traj_Skip_Idx(N_Traj),stat=STATUS)
    VERIFY_(STATUS)

    call Get_Traj_Idx(CallsPerWindow+1, N_TRAJ_SKIP, N_TRAJ)

    IF (MEMTRAJ .eq. 1) then
       allocate(STATE%Traj3(N_Traj),stat=STATUS)
       VERIFY_(STATUS)
       allocate(STATE%Traj2(N_Traj),stat=STATUS)
       VERIFY_(STATUS)
       do i=1, N_Traj
          STATE%Traj3(i) = ESMF_FieldBundleCreate(name="Trajectory", __RC__)
          call ESMF_FieldBundleSet(STATE%Traj3(i), grid=CUGrid, __RC__)
          STATE%Traj2(i) = ESMF_FieldBundleCreate(name="Trajectory_P2", __RC__)
          call ESMF_FieldBundleSet(STATE%Traj2(i), grid=CUGrid,    __RC__)
       enddo
    ELSE
       allocate(STATE%Traj3(3),stat=STATUS)
       VERIFY_(STATUS)
       allocate(STATE%Traj2(3),stat=STATUS)
       VERIFY_(STATUS)
       do i=1, 3
          STATE%Traj3(i) = ESMF_FieldBundleCreate(name="Trajectory", __RC__)
          call ESMF_FieldBundleSet(STATE%Traj3(i), grid=CUGrid, __RC__)
          STATE%Traj2(i) = ESMF_FieldBundleCreate(name="Trajectory_P2", __RC__)
          call ESMF_FieldBundleSet(STATE%Traj2(i), grid=CUGrid,    __RC__)
       enddo
    ENDIF

! Read into memory the non-linear trajectory at the beginning of
!  each subinterval in the window.
!-----------------------------------------------------------------

    do i=1, CallsPerWindow+1
       TIME_lst(i)  = Time
       Time = Time + DT*State%StepsPerCall
    end do

    do i=1, N_Traj

       IF (MEMTRAJ .eq. 1) then
          call Read_Traj(STATE%Traj3(I), STATE%Traj2(I), TIME_lst(Traj_Skip_Idx(i)), IAm)
    
          ! Also convert the true potential temperature in the trajectory file 
          !  to FV's scaled PT.
          !-------------------------------------------------------------------
  
          call ESMFL_BundleGetPointerToData(STATE%Traj3 (I), 'PT', PT, rc=STATUS)
          VERIFY_(STATUS)

          PT = PT/(MAPL_P00**MAPL_KAPPA)
       ENDIF

    end do

    call MAPL_TimerOff(MAPL,"--ReadTraj")


! Call Initialize for every Child
!--------------------------------

    call MAPL_TimerOff(MAPL,"TOTAL")

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, State%Clock,  RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn (MAPL,"TOTAL")

!  Create transforms between cubed and lat-lon grids.
!---------------------------------------------------
    call MAPL_TimerOn (MAPL,"--Transform")

    if (DO_ESMF_INTERP==1) then
       if(MAPL_AM_I_ROOT()) write(*,*) "Using ESMF TL/AD Interpolation Routines."
       state%L2C => new_regridder_manager%make_regridder(LLGrid, CUGrid, REGRID_METHOD_BILINEAR, &
            & rc=status); VERIFY_(status)
       state%C2L => new_regridder_manager%make_regridder(CUGrid, LLGrid, REGRID_METHOD_BILINEAR, &
            & rc=status); VERIFY_(status)
    else
       state%L2C => regridder_manager%make_regridder(LLGrid, CUGrid, REGRID_METHOD_BILINEAR, &
            & rc=status); VERIFY_(status)
       state%C2L => regridder_manager%make_regridder(CUGrid, LLGrid, REGRID_METHOD_BILINEAR, &
            & rc=status); VERIFY_(status)
    endif
    state%L2Ct = TransposeRegridder(state%L2C)
    state%C2Lt = TransposeRegridder(state%C2L)

    call MAPL_TimerOff(MAPL,"--Transform")

! All Done
!---------

    call MAPL_TimerOff(MAPL,"TOTAL")
    call MAPL_TimerOff(MAPL,"INITIALIZE")

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize



!BOP

! !IROUTINE: Run

! !DESCRIPTION: This is the first Run stage of FV. It is the container
!    for the dycore calculations. Subroutines from the core are 
!    invoked to do most of the work. A second run method, descibed below,
!    adds the import tendencies from external sources to the FV 
!    variables. 
!
!    In addition to computing and adding all dynamical contributions
!    to the FV variables (i.e., winds, pressures, and temperatures),
!    this method advects an arbitrary number of  tracers. These appear
!    in a ``Friendly'' bundle in the IMPORT state and are updated with
!    the advective tendency.
!
!
! !INTERFACE:

subroutine Run(gc, import, export, clock, rc)

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: gc
  type (ESMF_State),   intent(inout) :: import
  type (ESMF_State),   intent(inout) :: export
  type (ESMF_Clock),   intent(inout) :: clock
  integer          ,   intent(  out) :: rc  ! return code

!EOP

    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: IAm
    character(len=ESMF_MAXSTR)         :: COMP_NAME
    
! !Local Variables:
  

    type (T_AGCMPert_Wrap)             :: WRAP
    type (T_AGCMPert_State),   pointer :: STATE
    type (MAPL_MetaComp),      pointer :: MAPL 
  
    type (ESMF_GridComp),      pointer :: GCS(:)
    type (ESMF_State),         pointer :: GIM(:)
    type (ESMF_State),         pointer :: GEX(:)

    integer                            :: ADD_INCS
    integer                            :: CallsPerWindow
    logical                            :: AddIncs
    logical                            :: MakeExp
    integer                            :: Phase, currPhase
    integer                            :: Interval
    integer                            :: icall
    integer                            :: n, nc, nn
    integer                            :: i
    type(T_3DVAR), pointer             :: Traj3(:), Pert(:), Dtraj(:), Dpert(:)
    type(T_3DVAR)                      :: varsGSI(NUM_GSIvars)

    type(T_3DVAR), allocatable         :: varsTRAJ3(:)
    type(T_2DVAR), allocatable         :: varsTRAJ2(:)

    character(len=ESMF_MAXSTR)         :: PHASE_NAME
    character(len=ESMF_MAXSTR)         :: Timer
    class (AbstractRegridder), pointer :: TransIn, TransOut
    class (AbstractRegridder), pointer :: TransInt, TransOutt

    real, pointer                      :: PT(:,:,:)
    real :: d1, d2
    character(len=ESMF_MAXSTR)         :: lstvars

    integer                            :: interval_cur, interval_new
    integer,save                       :: interval_old_cur

! Begin
!------

    Iam = "Run"
    call ESMF_GridCompGet( GC, name=COMP_NAME, currentPhase=currPhase, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the generic state
! -----------------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

    if (DO_MOIST_PHYS == 0) then
       allocate(varsTRAJ3(NUM_GCM3DvarsRead))
    else
       allocate(varsTRAJ3(NUM_GCM3DvarsReadMoist))
    endif
    allocate(varsTRAJ2(NUM_GCM2Dvars))

    call MAPL_TimerOn(MAPL,"TOTAL")

    call timing_on('TOTAL')
    
! Retrieve the pointer to the internal state
! -----------------------------------------

    call ESMF_UserCompGetInternalState(gc, 'PERTstate', wrap, status)
    VERIFY_(STATUS)

    State => wrap%state

    ASSERT_(associated(State))

    select case(phase_opAdjoint(currPhase))
      case(.false.)
         PHASE_NAME = "TLM"
	 phase      = TLMphase
      case(.true.)
         PHASE_NAME = "ADJ"
	 phase      = ADJphase
      case default
         ASSERT_(.false.)
    end select

    call MAPL_TimerOn(MAPL,"--RUN"//PHASE_NAME)

! Handy temporary to hold pointers
!---------------------------------

    do i=1,size(varsGSI)
       varsGSI(i)%name = trim(GSIvars(i))
    end do


    if (DO_MOIST_PHYS == 0) then
       do i=1,size(varsTRAJ3)
          varsTRAJ3(i)%name = trim(GCM3DvarsRead(i))
       end do
    else
       do i=1,size(varsTRAJ3)
          varsTRAJ3(i)%name = trim(GCM3DvarsReadMoist(i))
       end do
    endif

    do i=1,size(varsTRAJ2)
       varsTRAJ2(i)%name = trim(GCM2Dvars(i))
    end do


! If we are starting with an ADJ, we pretend we just 
!  finished a TLM pass by resetting the call counter.
!  If we start with a TLM, we leave it at zero.
!----------------------------------------------------

    CallsPerWindow = size(TIME_lst)-1


    if(State%CallCounter==0 .and. phase==ADJphase) then
       State%CallCounter = CallsPerWindow
    endif

! Intervals within the window are numbered from 1 to CallsPerWindow
!------------------------------------------------------------------

    if(phase==TLMphase) then
       Interval = mod(State%CallCounter,CallsPerWindow) + 1

       if(Interval==1) then
          call ESMF_ClockSet( State%clock, direction=ESMF_DIRECTION_FORWARD, rc=status )
          VERIFY_(STATUS)

          do i=1,size(State%PERTwrk)
             State%PERTwrk(i)%X = 0.0
          end do
       end if

       transIn => state%L2C
       transOut => state%C2L

       transInt => state%L2Ct
       transOutt => state%C2Lt
    else
       Interval = CallsPerWindow - mod(State%CallCounter,CallsPerWindow)

       if(Interval==CallsPerWindow) then
          call ESMF_ClockSet( State%Clock, direction=ESMF_DIRECTION_REVERSE, rc=status )
          VERIFY_(STATUS)

          do i=1,size(State%PERTwrk)
             State%PERTwrk(i)%X = 0.0
          end do
       end if

       transIn => state%C2Lt
       transOut => state%L2Ct
       TransInt => State%C2L
       TransOutt => State%L2C

    end if

! Update the call counter
!------------------------

    State%CallCounter = State%CallCounter + 1

! Adding import-state increments to the perturbation is controlled by the
!  ADD_INCS resource.
!  If this resource is 0, increments will be added according to currPhase;
!  if it is 1, increments will be added only during the TLM phase;
!  if it is 2, increments will be added only during the ADJ phase;
!  if it is 3, increments will be added during both phases;
!  and for the DEFAULT case -1 increments will be added during the first 
!  interval of the TLM phase and at every interval of the ADJ phase. 
!------------------------------------------------------------------------

    call MAPL_GetResource ( MAPL, ADD_INCS, 'ADD_INCS:', default=-1, rc=STATUS )
    VERIFY_(STATUS)

    select case(ADD_INCS)

    case(-1)
       AddIncs = phase==ADJphase .or. Interval==1
       MakeExp = .true.

    case(0)
       AddIncs = phase_addImport(currphase)
       MakeExp = phase_getExport(currphase)
       
    case(1)
       AddIncs = phase==TLMphase
       MakeExp = .true.

    case(2)
       AddIncs = phase==ADJphase
       MakeExp = .true.

    case(3)
       AddIncs = .true.
       MakeExp = .true.

    case default
       ASSERT_(.false.)

    end select

    if(MAPL_AM_I_ROOT()) then
      print'(1x,2a,i2,a,b3.3,a,i6,4x,2l2)', trim(Iam), &
      	': phase, interval, Addincs, MakeExp = ', currphase,'(',currphase-1,')', interval, AddIncs, makeExp
      call flush(6)
    endif

! Get children and their im/ex states from my generic state.
!----------------------------------------------------------

    call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX, RC=STATUS )
    VERIFY_(STATUS)

! The Import State increments use the same names as the GSI state variable.
!  Note that this is an increment not a tendency
!  and is added only once per Interval (Call) before going into the
!  model's time loop for both TLM and ADJ. TransIn has the inverse
!  and transpose info for the GSI2GCM (the "T" operator) call.
!--------------------------------------------------------------------
    if(phase==TLMphase) then
      if (Interval.eq.1) then
          Interval_old_cur = -99
      endif
    else
      if (Interval.eq.CallsPerWindow) then
          Interval_old_cur = -99
      endif
    endif

    if(phase==TLMphase) then
       IF (MEMTRAJ .ne. 0) then
          Interval_cur = Traj_Orgi_Idx(Interval)
          Interval_new = Traj_Orgi_Idx(Interval+1)
          call PERT_RefVarsFromBundle   (State%Traj3(Interval_cur), varsTRAJ3, rc=STATUS)
          VERIFY_(STATUS)
          call PERT_RefVarsFromBundle_2D(State%Traj2(Interval_cur), varsTRAJ2, rc=STATUS)
          VERIFY_(STATUS)
       ELSE          
          Interval_cur = Traj_Skip_Idx(Traj_Orgi_Idx(Interval))
          Interval_new = Traj_Skip_Idx(Traj_Orgi_Idx(Interval+1))
          if (Interval_old_cur .ne. Interval_cur) then
             call Read_Traj(State%Traj3(1), State%Traj2(1), Time_lst(Interval_cur), IAm)
             call ESMFL_BundleGetPointerToData(State%Traj3(1), 'PT', PT, rc=STATUS)
             PT = PT/(MAPL_P00**MAPL_KAPPA)
          else
             State%Traj3(1) = State%Traj3(3)
             State%Traj2(1) = State%Traj2(3)
          endif

          if (Interval_new.ne.Interval_cur) then
             call Read_Traj(State%Traj3(2), State%Traj2(2), Time_lst(Interval_new), IAm)
             call ESMFL_BundleGetPointerToData(State%Traj3(2), 'PT', PT, rc=STATUS)
             PT = PT/(MAPL_P00**MAPL_KAPPA)
          endif

          call PERT_RefVarsFromBundle   (State%Traj3(1), varsTRAJ3, rc=STATUS)
          VERIFY_(STATUS)
          call PERT_RefVarsFromBundle_2D(State%Traj2(1), varsTRAJ2, rc=STATUS)
          VERIFY_(STATUS)
       ENDIF
    else
       IF (MEMTRAJ .ne. 0) then
          Interval_cur = Traj_Orgi_Idx(Interval)
          Interval_new = Traj_Orgi_Idx(Interval+1)
          call PERT_RefVarsFromBundle   (State%Traj3(Interval_new), varsTRAJ3, rc=STATUS)
          VERIFY_(STATUS)
          call PERT_RefVarsFromBundle_2D(State%Traj2(Interval_new), varsTRAJ2, rc=STATUS)
          VERIFY_(STATUS)
       ELSE
          Interval_cur = Traj_Skip_Idx(Traj_Orgi_Idx(Interval))
          Interval_new = Traj_Skip_Idx(Traj_Orgi_Idx(Interval+1))
          if (Interval_old_cur .ne. Interval_cur) then
             call Read_Traj(State%Traj3(1), State%Traj2(1), Time_lst(Interval_cur), IAm)
             call ESMFL_BundleGetPointerToData(State%Traj3(1), 'PT', PT, rc=STATUS)
             PT = PT/(MAPL_P00**MAPL_KAPPA)
          else
             State%Traj3(1) = State%Traj3(3)
             State%Traj2(1) = State%Traj2(3)
          endif

          if (Interval_new.ne.Interval_cur) then
             call Read_Traj(State%Traj3(2), State%Traj2(2), Time_lst(Interval_new), IAm)
             call ESMFL_BundleGetPointerToData(State%Traj3(2), 'PT', PT, rc=STATUS)
             PT = PT/(MAPL_P00**MAPL_KAPPA)
          endif

          if (Interval_new.eq.Interval_cur) then
             call PERT_RefVarsFromBundle   (State%Traj3(1), varsTRAJ3, rc=STATUS)
             VERIFY_(STATUS)
             call PERT_RefVarsFromBundle_2D(State%Traj2(1), varsTRAJ2, rc=STATUS)
             VERIFY_(STATUS)
          else
             call PERT_RefVarsFromBundle   (State%Traj3(2), varsTRAJ3, rc=STATUS)
             VERIFY_(STATUS)
             call PERT_RefVarsFromBundle_2D(State%Traj2(2), varsTRAJ2, rc=STATUS)
             VERIFY_(STATUS)
          endif
       ENDIF
    endif

    if(AddIncs) then
       call PERT_RefVarsFromState(IMPORT, varsGSI, rc=status )
       VERIFY_(STATUS)

       do i=1,size(State%PERTdel)
          State%PERTdel(i)%X  = 0.
       end do

       call timing_on('GSI2GCM')
       call GSI2GCM(varsGSI, State%PERTdel, varsTRAJ3, varsTRAJ2, TransIn, rc=STATUS)
       VERIFY_(STATUS)
       call timing_off('GSI2GCM')

       do i=1,size(State%PERTwrk)
          State%PERTwrk(i)%X  = State%PERTwrk(i)%X  + State%PERTdel(i)%X
       end do

#ifdef DO_DOT_TEST
!---dot test for TransIn, only if VarsGSI has been associated to IMPORT ----------
    call Dot_TransIn(VarsGSI,State%PERTwrk,VarsTRAJ3,varsTRAJ2,TransIn,TransInt,PHASE_NAME,rc)
    VERIFY_(status)
!---------------------------------------------------------------------------------
#endif 
    end if

! At this point, we need the trajectory at the beginning
!  of the interval in the working trajectory, regardless
!  of phase.
!-------------------------------------------------------

#ifndef NO_TIME_LOOP_
       IF (MEMTRAJ .ne. 0) THEN
          call PERT_CopyVarsFromBundle   (State%Traj3(Interval_cur), State%TRAJwrk3, rc=STATUS)
          VERIFY_(STATUS)
          call PERT_CopyVarsFromBundle_2D(State%Traj2(Interval_cur), State%TRAJwrk2, rc=STATUS)
          VERIFY_(STATUS)
       ELSE
          call PERT_CopyVarsFromBundle   (State%Traj3(1), State%TRAJwrk3, rc=STATUS)
          VERIFY_(STATUS)
          call PERT_CopyVarsFromBundle_2D(State%Traj2(1), State%TRAJwrk2, rc=STATUS)
          VERIFY_(STATUS)
       ENDIF
#endif /* NO_TIME_LOOP_ */

#ifdef DO_ADTEST_
    icall=0
    if(phase==ADJphase) icall=1
    call adtester_test(State%PERTwrk,icall,interval=interval,phase=phase)
#endif
! Loop over the GCM steps per interval calling all the
!  children each time. All children update the 
!  working trajectory and perturbation. We also apply an
!  optional "physics" tendency to the trajectory. Note that
!  during the ADJ phase the order of calling for the children
!  is reversed.
!------------------------------------------------------------

#ifndef NO_TIME_LOOP_

    call timing_on('MODEL TOTAL')

    TIME_LOOP: do n=1,State%StepsPerCall

       !Controls for iterative traj storing
       if (cp_iter_controls%cp_i .ne. 0) then
          cp_iter_controls%cp_t = cp_iter_controls%cp_t + 1 !Advance timestep counter
       endif

       CHILDREN: do nc = 1, size(GCS)

          if(phase==TLMphase) then
             nn = nc
          else
             nn = size(GCS) + 1 - nc
          end if

          if(nn==DYNA) then
             Timer = "----RunDyna"//PHASE_NAME
          else
             Timer = "----RunPhys"//PHASE_NAME
          end if

          call MAPL_TimerOn(MAPL,Timer)
          call ESMF_GridCompRun (GCS(nn), importState=GIM(nn), exportState=GEX(nn), &
                                 clock=State%Clock, phase=phase, userRC=STATUS )
          VERIFY_(STATUS)
          call MAPL_TimerOff(MAPL,Timer)

       end do CHILDREN

! Advance the local clock
!------------------------

       call ESMF_ClockAdvance(State%Clock, RC=STATUS )
       VERIFY_(STATUS)

       if (cp_iter_controls%cp_i .ne. 0 .and. cp_iter_controls%cp_t==cp_iter_controls%cp_nt) then
          if (phase==ADJphase) cp_iter_controls%cp_i = cp_iter_controls%cp_i + 1 !Advance iteration counter
          cp_iter_controls%cp_t = 0 !Reset time step counter
       endif

    end do TIME_LOOP

    call timing_off('MODEL TOTAL')

#endif /* NO_TIME_LOOP_ */

#ifdef DO_ADTEST_
    icall=mod(icall+1,2)	! toggle icall values in 0 and 1
    call adtester_test(State%PERTwrk,icall,interval=interval,phase=phase)
#endif
! Fill the EXPORT perturbation.For the ADJphase, this requires the trajectory
!  at the beginning of the interval. For the TLMphase we need the trajectory 
!  at the end of the interval.
!-----------------------------------------------------------------------------

    if(phase==TLMphase) then
      IF (MEMTRAJ .ne. 0) THEN
         call PERT_RefVarsFromBundle   (State%Traj3(Interval_new), varsTRAJ3, rc=STATUS)
         VERIFY_(STATUS)
         call PERT_RefVarsFromBundle_2D(State%Traj2(Interval_new), varsTRAJ2, rc=STATUS)
         VERIFY_(STATUS)
      ELSE
         if (Interval_new.eq.Interval_cur) then
            call PERT_RefVarsFromBundle   (State%Traj3(1), varsTRAJ3, rc=STATUS)
            VERIFY_(STATUS)
            call PERT_RefVarsFromBundle_2D(State%Traj2(1), varsTRAJ2, rc=STATUS)
            VERIFY_(STATUS)
         else
            call PERT_RefVarsFromBundle   (State%Traj3(2), varsTRAJ3, rc=STATUS)
            VERIFY_(STATUS)
            call PERT_RefVarsFromBundle_2D(State%Traj2(2), varsTRAJ2, rc=STATUS)
            VERIFY_(STATUS)
         endif
      ENDIF
    else
      IF (MEMTRAJ .ne. 0) THEN
         call PERT_RefVarsFromBundle   (State%Traj3(Interval_cur), varsTRAJ3, rc=STATUS)
         VERIFY_(STATUS)
         call PERT_RefVarsFromBundle_2D(State%Traj2(Interval_cur), varsTRAJ2, rc=STATUS)
         VERIFY_(STATUS)
      ELSE
         call PERT_RefVarsFromBundle   (State%Traj3(1), varsTRAJ3, rc=STATUS)
         VERIFY_(STATUS)
         call PERT_RefVarsFromBundle_2D(State%Traj2(1), varsTRAJ2, rc=STATUS)
         VERIFY_(STATUS)
      ENDIF
    endif

    if(MakeExp) then
! Put the EXPORT pointers of the perturbation in VarsGSI.
!--------------------------------------------------------

      call PERT_RefVarsFromState(EXPORT, varsGSI, rc=status )
      VERIFY_(STATUS)

      call pertWrk_save_(State%PERTwrk)
      call timing_on('GCM2GSI')
      call GSI2GCM(varsGSI, State%PERTwrk, varsTRAJ3, varsTRAJ2, TransOut, rc=STATUS)
      VERIFY_(STATUS)
      call timing_off('GCM2GSI')
      call pertWrk_restore_(State%PERTwrk)

#ifdef DO_DOT_TEST
!---dot test lines for TransOut, only if VarsGSI has been associated to EXPORT ------
    call Dot_TransOut(VarsGSI,State%PERTwrk,VarsTRAJ3,VarsTRAJ2,TransOut,TransOutt,PHASE_NAME,rc)
!------------------------------------------------------------------------------------
#endif 
    endif	! (MakeExp)

    IF (MEMTRAJ .eq. 0) then
       State%Traj3(3) = State%Traj3(1)
       State%Traj2(3) = State%Traj2(1)
       Interval_old_cur = Interval_cur
    ENDIF

! All done
! --------
    call MAPL_TimerOff(MAPL,"--RUN"//PHASE_NAME)
    call MAPL_TimerOff(MAPL,"TOTAL")

#ifdef VERBOSE
    if(MAPL_AM_I_ROOT()) print *, 'Done with ', trim(PHASE_NAME)
#endif /* VERBOSE */

    deallocate(varsTRAJ3)
    deallocate(varsTRAJ2)

    call timing_off('TOTAL')

    if(phase==TLMphase .and. .not. BothPhases) then
      if (Interval.eq.CallsPerWindow) then
          call timing_prt( mpp_pe() )
      endif
    elseif (phase==ADJphase) then
      if (Interval.eq.1) then
          call timing_prt( mpp_pe() )
      endif
    endif

    RETURN_(ESMF_SUCCESS)
  end subroutine Run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !IROUTINE: Finalize -- Finalize method for the composite Agcm Gridded Component

! !INTERFACE:

  subroutine Finalize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer         ,    intent(  out) :: RC     ! Error code

! !DESCRIPTION: 
 

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm 
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (T_AGCMPert_Wrap)              :: WRAP
    type (T_AGCMPert_State),   pointer  :: STATE
    type (MAPL_MetaComp),      pointer  :: MAPL

! =============================================================================

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Finalize"

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"FINALIZE")

! Call Finalize for every Child
!--------------------------------

    call MAPL_GenericFinalize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")

! Retrieve the pointer to the internal state
! ------------------------------------------

    call ESMF_UserCompGetInternalState(gc, 'PERTstate', wrap, status)
    VERIFY_(STATUS)

    State => wrap%state

    ASSERT_(associated(State))

! Destroy and deallocate al resources and memory 
!  associated with this instances internal state.
!-----------------------------------------------

    call PertState_DeepDestroy(State,rc=STATUS) 
    VERIFY_(STATUS)

    deallocate(State)
    deallocate(TIME_lst)
    deallocate(TRAJ_FILE_TMPL)
    deallocate(Traj_Orgi_Idx)
    deallocate(Traj_Skip_Idx)

!Destroy iterative traj container
!--------------------------------

    if (cp_iter_controls%cp_i .ne. 0) call finalize_cp_iter

! All done
! --------

    call MAPL_TimerOff(MAPL,"TOTAL")
    call MAPL_TimerOff(MAPL,"FINALIZE")


    call WRITE_PARALLEL(" ")
    call WRITE_PARALLEL(" Times for "//trim(COMP_NAME))
    !call MAPL_ProfWrite(MAPL%TIMES,RC=STATUS)
    !VERIFY_(STATUS)
    call WRITE_PARALLEL(" ")


    RETURN_(ESMF_SUCCESS)
  end subroutine Finalize


  subroutine pertWrk_save_(pertWrk)
    implicit none
    type(T_3DVar), dimension(:), intent(INOUT) :: pertWrk

    integer:: mx,my,mz,iv,mv

    mv=size(pertWrk)
    mx=0; my=0; mz=0
    if(mv>0) then
      mx=size(pertwrk(1)%x,1)
      my=size(pertwrk(1)%x,2)
      mz=size(pertwrk(1)%x,3)
    endif

    allocate(pertWrk_saved_(mx,my,mz,mv))

    do iv=1,mv
      pertWrk_saved_(:,:,:,iv) = pertwrk(iv)%x(:,:,:)
    enddo
  end subroutine pertWrk_save_

  subroutine pertWrk_restore_(pertWrk)
    implicit none
    type(T_3DVar), dimension(:), intent(INOUT) :: pertWrk

    integer:: iv,mv

    mv=size(pertwrk)
    do iv=1,mv
      pertWrk(iv)%x(:,:,:) = pertWrk_saved_(:,:,:,iv)
    enddo
    deallocate(pertWrk_saved_)
  end subroutine pertWrk_restore_

!=============================================================================

  subroutine PertState_DeepDestroy(State,RC)

    type (T_AGCMPert_State),  intent(INOUT) :: STATE
    integer,                  intent(  OUT) :: RC

    integer :: i, STATUS
    character(len=ESMF_MAXSTR) :: IAm="PertState_DeepDestroy" 

    !DEEP DESTORY RANK 2 TRAJECTORY
       do i=1,size(State%Traj2)-1
          call Bundle_DeepDestroy_2D(State%Traj2(i),rc=STATUS)
          VERIFY_(STATUS)
       end do   
       call Bundle_DeepDestroy_2D(State%Traj2(size(State%Traj2)),rc=STATUS)
       VERIFY_(STATUS)

    if (gq_ntrj > 1) then
        call Bundle_DeepDestroy_2D(TRAJwrk2_4gq,rc=STATUS)
        VERIFY_(STATUS)
    endif

       do i=1,size(State%TrajWrk2)
          deallocate(State%TrajWrk2(i)%X)
       end do


       do i=1,size(State%Traj3)-1
          call Bundle_DeepDestroy(State%Traj3(i),rc=STATUS)
          VERIFY_(STATUS)
       end do

       call Bundle_DeepDestroy(State%Traj3(size(State%Traj3)),rc=STATUS)
       VERIFY_(STATUS)

    if (gq_ntrj > 1) then
        call Bundle_DeepDestroy(TRAJwrk3_4gq,rc=STATUS)
        VERIFY_(STATUS)
    endif

!----
    do i=1,size(State%TrajWrk3)
       deallocate(State%TrajWrk3(i)%X)
    end do
       deallocate(State%TrajWrk3)
!----
    do i=1,size(State%PertWrk)
       deallocate(State%PertWrk(i)%X)
    end do
       deallocate(State%PertWrk)
!----
    do i=1,size(State%PertDel)
       deallocate(State%PertDel(i)%X)
    end do
       deallocate(State%PertDel)
!----
    do i=1,size(State%TrajWrk2)
       deallocate(State%TrajWrk2(i)%X)
    end do
       deallocate(State%TrajWrk2)
!----

    call regridder_manager%delete_regridder(state%L2C)
    call regridder_manager%delete_regridder(state%C2L)

    RETURN_(ESMF_SUCCESS)
  end subroutine PertState_DeepDestroy


!!!===================================================================

  subroutine Bundle_DeepDestroy(Bundle,RC)
    type (ESMF_FieldBundle),  intent(INOUT) :: BUNDLE
    integer,              intent(  OUT) :: RC

    integer           :: n, nq, STATUS
    real, pointer     :: X(:,:,:)
    type (ESMF_Field) :: Field

    character(len=ESMF_MAXSTR) :: IAm="Bundle_DeepDestroy" 

    call ESMF_FieldBundleGet (Bundle, fieldCount=NQ, RC=STATUS )
    VERIFY_(STATUS)

    do N=1,NQ
       call ESMF_FieldBundleGet(Bundle, n, Field, RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_FieldGet(Field, 0, X, RC=STATUS)
       VERIFY_(STATUS)

       if(associated(X)) deallocate(X) 

       call ESMF_FieldDestroy(Field, RC=STATUS) 
       VERIFY_(STATUS)
    end do

    call ESMF_FieldBundleDestroy(Bundle, RC=STATUS) 
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine Bundle_DeepDestroy

  subroutine Bundle_DeepDestroy_2D(Bundle,RC)
    type (ESMF_FieldBundle),  intent(INOUT) :: BUNDLE
    integer,              intent(  OUT) :: RC

    integer           :: n, nq, STATUS
    real, pointer     :: X(:,:)
    type (ESMF_Field) :: Field

    character(len=ESMF_MAXSTR) :: IAm="Bundle_DeepDestroy_2D" 

    call ESMF_FieldBundleGet (Bundle, fieldCount=NQ, RC=STATUS )
    VERIFY_(STATUS)

    do N=1,NQ
       call ESMF_FieldBundleGet(Bundle, n, Field, RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_FieldGet(Field, 0, X, RC=STATUS)
       VERIFY_(STATUS)

       if(associated(X)) deallocate(X) 

       call ESMF_FieldDestroy(Field, RC=STATUS) 
       VERIFY_(STATUS)
    end do

    call ESMF_FieldBundleDestroy(Bundle, RC=STATUS) 
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine Bundle_DeepDestroy_2D


! This is the complete GSI to GCM Operator.


  subroutine GSI2GCM(VarsGSI, VarsGCM, VarsTRAJ3, varsTRAJ2, TR, rc)
    type(T_3Dvar),                  intent(INOUT) :: VarsGSI(:)
    type(T_3Dvar),                  intent(INOUT) :: VarsGCM(:)
    type(T_3Dvar), target,          intent(IN   ) :: VarsTRAJ3(:)
    type(T_2Dvar), target,          intent(IN   ) :: VarsTRAJ2(:)
    class (AbstractRegridder),      intent(in) :: TR
    integer,     optional,          intent(  OUT) :: rc

! ErrLog vars

    integer                    :: STATUS
    character(len=ESMF_MAXSTR) :: IAm="GSI2GCM"

! Local vars

    integer                    :: IM, JM, LM, i,j,k,ii
    logical                    :: Transpose
    logical                    :: FillGSI
    logical                    :: Inverse
    character*60               :: GridTypeIn

    real, allocatable, target  :: TVP(:,:,:)
    real, allocatable, target  :: UAP(:,:,:)
    real, allocatable, target  :: VAP(:,:,:)

    real, allocatable  ::  BX(:,:,:),  CX(:,:,:)
    real, allocatable  ::  AY(:,:,:),  BY(:,:,:),  CY(:,:,:)
    real, allocatable  :: PKP(:,:,:), PKT(:,:,:)
    real, allocatable  :: PET(:,:,:), PKE(:,:,:)
    real, allocatable  :: PrT(:,:,:), DPT_Ps(:,:,:)

    real, pointer      :: PTP(:,:,:), PTT(:,:,:)
    real, pointer      :: QVP(:,:,:), QVT(:,:,:)
    real, pointer      :: DPP(:,:,:)

    real, pointer      :: PS(:,:)

    real(8),  pointer     :: ak(:), bk(:)


! Begin
!------

! Dimensions of GCM variables
!----------------------------

    IM = size(varsTRAJ3(1)%X,1)
    JM = size(varsTRAJ3(1)%X,2)
    LM = size(varsTRAJ3(1)%X,3)

! Set the modes of behavior. Inverse=T means GCM to GSI,
!  Transpose=T uses the transpose of all operators, and
!  FillGSI means that the GSI variables are the outputs,
!  which is the case when Inverse-Transpose are T-F or F-T;
!  for T-T or F-F the GCM variables are the output.
!----------------------------------------------------------

    
    transpose = TR%isTranspose()
    block
       type(RegridderSpec) :: spec
       spec = TR%get_spec()
       if (transpose) then
          call spec%get_grid_type(OutputGridType=GridTypeIn,rc=status)
          VERIFY_(status)
       else
          call spec%get_grid_type(InputGridType=GridTypeIn,rc=status)
          VERIFY_(status)
       end if
    end block
    FillGSI   = GridTypeIn == 'Cubed-Sphere'

! Allocate space for temporary arrays on the GCM grid
!----------------------------------------------------

    allocate(UAP(-2:IM+3,-2:JM+4,LM), stat=status); VERIFY_(STATUS)
    allocate(VAP(-2:IM+4,-2:JM+3,LM), stat=status); VERIFY_(STATUS)
    allocate(TVP(IM,JM,1:LM)        , stat=status); VERIFY_(STATUS)
    allocate(PKT(IM,JM,1:LM)        , stat=status); VERIFY_(STATUS)
    allocate(PKE(IM,JM,0:LM)        , stat=status); VERIFY_(STATUS)
    allocate(PET(IM,JM,0:LM)        , stat=status); VERIFY_(STATUS)
    allocate(BX (IM,JM,1:LM)        , stat=status); VERIFY_(STATUS)
    allocate(CX (IM,JM,1:LM)        , stat=status); VERIFY_(STATUS)
    allocate(PKP(IM,JM,1:LM)        , stat=status); VERIFY_(STATUS)
    allocate(AY (IM,JM,1:LM)        , stat=status); VERIFY_(STATUS)
    allocate(BY (IM,JM,1:LM)        , stat=status); VERIFY_(STATUS)
    allocate(CY (IM,JM,1:LM)        , stat=status); VERIFY_(STATUS)

    allocate(PrT(IM,JM,0:LM)        , stat=status); VERIFY_(STATUS)
    allocate(DPT_Ps(IM,JM,1:LM)        , stat=status); VERIFY_(STATUS)
    allocate(ak(size(pert_ak)),      stat=status); VERIFY_(STATUS)
    allocate(bk(size(pert_bk)),      stat=status); VERIFY_(STATUS)

! Aliases for the thermodynamic trajectory varibles
!--------------------------------------------------

    PTT => varsTRAJ3(PERT_GetIndex(varsTRAJ3,"PT"))%X
    QVT => varsTRAJ3(PERT_GetIndex(varsTRAJ3,"QV"))%X

    PS  => varsTRAJ2(PERT_GetIndex2D(varsTRAJ2,"PS"))%X

    if (DynTestCase > 0) PS = 1.e5

! Compute reference pressure from Ps
!-----------------------------------
    ak = pert_ak
    bk = pert_bk
    do i = 0,lm
       PrT(:,:,i) = AK(i+1) + BK(i+1)*PS(:,:)
    enddo

! Compute reference pressure thickness
!-------------------------------------
    DPT_Ps(:,:,1:LM) = PrT(:,:,1:LM) - PrT(:,:,0:LM-1)

! Compute trajectory factors based on DPT. The trajectory
!  is always on the GCM's Cubed sphere grid.
!--------------------------------------------------------

    call GetPressVarsTraj(DPT_Ps, PET, PKE, PKT, BX, CX)

! These factors are used only on the GSI to GCM direction
!--------------------------------------------------------

    AY  = 1.0 + EPS*QVT
    BY  = PTT*EPS/AY
    CY  = PTT/PKT
    AY  = 1.0/(PKT*AY)

! Aliases for the perturbation variables involved in the
!  creation of ancillary GCM variables from the GSI set.
!-------------------------------------------------------

    PTP => VarsGCM(PERT_GetIndex(VarsGCM,"PT"))%X
    QVP => VarsGCM(PERT_GetIndex(VarsGCM,"QV"))%X
    DPP => VarsGCM(PERT_GetIndex(VarsGCM,"DP"))%X

! Create perturbation PK from DP and trajectory factors
!------------------------------------------------------

    select type (TR)
    type is (TransposeRegridder) ! not a transpose
       transpose = .true.
    class default
       transpose = .false.
    end select

    if(.not.Transpose) then
       call DP2PK( DPT_Ps, PKE, PET, BX, CX, DPP, PKP, Transpose)
    endif

! Do the transform. The change of variables is always on
!  the cube side. 
!-------------------------------------------------------

    A2MorM2A: if(FillGSI) then
       if(.not.Transpose) then
          TVP = PTT*PKT*EPS*QVP + (1.0 + EPS*QVT)*(PKP*PTT + PKT*PTP)
          if (DynTestCase > 0) TVP = PKT*PTP
       else
          TVP =       AY*PTP
          QVP = QVP - BY*PTP
          PKP =     - CY*PTP
       end if

       call RegridVars(FillGSI,rc=status)
       VERIFY_(status)
    else
       call RegridVars(FillGSI,rc=status)
       VERIFY_(status)

       if(.not.Transpose) then
          PTP = AY*TVP - BY*QVP - CY*PKP
       else
          PTP =     (1.0 + EPS*QVT)*PKT*TVP
          QVP = QVP +      EPS*PKT *PTT*TVP
          PKP =     (1.0 + EPS*QVT)*PTT*TVP
       end if
    end if A2MorM2A

! Update DP influences from perturbation PK and trajectory factors
!-----------------------------------------------------------------

    if(Transpose) then
       call DP2PK( DPT_Ps, PKE, PET, BX, CX,  DPP, PKP, Transpose)
    endif

! Clean up
!---------

    deallocate(UAP)
    deallocate(VAP)

    deallocate(BX )
    deallocate(CX )       
    deallocate(TVP)

    deallocate(PKE)                 
    deallocate(PKT)
    deallocate(PET)
    deallocate(PKP)

    deallocate(AY )
    deallocate(BY )
    deallocate(CY )

    deallocate(PrT)
    deallocate(DPT_Ps) 
    deallocate(ak) 
    deallocate(bk)

    RETURN_(ESMF_SUCCESS)

  contains

    subroutine RegridVars(OutputIsLL,rc)
      logical, intent(IN) :: OutputIsLL
      integer, optional, intent(OUT) :: rc

      integer       :: status
      character(len=ESMF_MAXSTR) :: Iam ="RegridVars"
      integer       :: i, nu, nv
      real, pointer :: VarGCM (:,:,:)
      real, pointer :: VarGCMu(:,:,:)
      real, pointer :: VarGCMv(:,:,:)
      
      do i=1,size(VarsGSI)
         if(associated(varsGSI(i)%X)) then
            select case(VarsGSI(i)%name)
            case("TV")
               VarGCM  => TVP
!#ifdef _PROPER_WINDS_
            case("U")
               VarGCMu => UAP(1:IM,1:JM,:)
               cycle
            case("V")
               VarGCMv => VAP(1:IM,1:JM,:)
               cycle
!#endif /* _PROPER_WINDS_ */
            case default
               VarGCM  => VarsGCM(PERT_GetIndex(VarsGCM,VarsGSI(i)%name))%X
            end select

            if (OutputIsLL) then
               call tr%regrid(varGCM, VarsGSI(i)%X, rc=status)
               VERIFY_(STATUS)

               call Flip3(VarsGSI(i)%X)
            else
               call Flip3(VarsGSI(i)%X)

               call tr%regrid(VarsGSI(i)%X, VarGCM, rc=status)
               VERIFY_(STATUS)
            end if
         else
            ASSERT_(OutputIsLL) ! GSI pointer can only be null when OutputIsLL
         end if
      end do

      nu = PERT_GetIndex(VarsGSI,"U")
      nv = PERT_GetIndex(VarsGSI,"V")

!#ifdef _PROPER_WINDS_
      DO_UV: if(associated(VarsGSI(nu)%X) .and. associated(VarsGSI(nv)%X)) then
         if(OutputIsLL) then 
            UAP(1:IM,1:JM,:) = VarsGCM(PERT_GetIndex(VarsGCM,"U"))%X
            VAP(1:IM,1:JM,:) = VarsGCM(PERT_GetIndex(VarsGCM,"V"))%X

            call ReStaggerWinds(UAP,VAP,D2A=.true.)

            call tr%regrid(VarGCMu, VarGCMv, VarsGSI(nu)%X, VarsGSI(nv)%X, rotate=.true., RC=STATUS)
            VERIFY_(STATUS)

            call Flip3( VarsGSI(nu)%X )
            call Flip3( VarsGSI(nv)%X )
         else
            call Flip3( VarsGSI(nu)%X   )
            call Flip3( VarsGSI(nv)%X   )

            call tr%regrid(VarsGSI(nu)%X, VarsGSI(nv)%X, VarGCMu, VarGCMv, rotate=.true., RC=STATUS)
            VERIFY_(STATUS)

            call ReStaggerWinds(UAP,VAP,D2A=.false.)

            VarsGCM(PERT_GetIndex(VarsGCM,"U"))%X = UAP(1:IM,1:JM,:)
            VarsGCM(PERT_GetIndex(VarsGCM,"V"))%X = VAP(1:IM,1:JM,:)
         end if
      end if DO_UV
      RETURN_(ESMF_SUCCESS)
!#endif /* _PROPER_WINDS_ */

    end subroutine RegridVars
    
    subroutine ReStaggerWinds(U, V, D2A)
      real,                      intent(INOUT) :: U(-2:,-2:,:)
      real,                      intent(INOUT) :: V(-2:,-2:,:)
      logical,                   intent(IN   ) :: D2A

      real ::  ebuffery(ubound(V,2)-3,size(V,3))
      real ::  wbuffery(ubound(V,2)-3,size(V,3))

      real ::  nbufferx(ubound(U,1)-3,size(U,3))
      real ::  sbufferx(ubound(U,1)-3,size(U,3))

      integer :: im, jm

      IM = ubound(U,1)-3
      JM = ubound(V,2)-3

      if(D2A) then
        call mpp_get_boundary(U, V, FV_Atm(1)%domain, ebuffery=ebuffery,&
                                                      nbufferx=nbufferx,&
                                                      wbuffery=wbuffery,&
                                                      sbufferx=sbufferx,&
                                                      gridtype=DGRID_NE)

        U(1:IM,JM+1,:) = nbufferx
        V(IM+1,1:JM,:) = ebuffery

        U(1:IM,1:JM,:) = 0.5*(U(1:IM,1:JM,:) + U(1:IM,2:JM+1,:))
        V(1:IM,1:JM,:) = 0.5*(V(1:IM,1:JM,:) + V(2:IM+1,1:JM,:))
      else
        U(1:IM,JM+1,:) = 0.5*(U(1:IM,  JM,:))
        U(1:IM,2:JM,:) = 0.5*(U(1:IM,2:JM,:) + U(1:IM,1:JM-1,:))
        U(1:IM,1   ,:) = 0.5*(U(1:IM,1   ,:))
      
        V(IM+1,1:JM,:) = 0.5*(V(  IM,1:JM,:))
        V(2:IM,1:JM,:) = 0.5*(V(2:IM,1:JM,:) + V(1:IM-1,1:JM,:))
        V(1   ,1:JM,:) = 0.5*(V(1   ,1:JM,:))

        call mpp_get_boundary(U,V,FV_Atm(1)%domain,ebuffery=ebuffery,&
                                                   nbufferx=nbufferx,&
                                                   wbuffery=wbuffery,&
                                                   sbufferx=sbufferx,&
                                                   gridtype=DGRID_NE)

         U(1:IM,   1,:) = U(1:IM,   1,:) + sbufferx
         V(1   ,1:JM,:) = V(1   ,1:JM,:) + wbuffery  
    endif

      return
    end subroutine ReStaggerWinds

    subroutine Flip3(A)
      real, intent(INOUT) :: A(:,:,:)

      real    :: B(size(A,1),size(A,2))
      integer :: LM, L
#ifdef _DOFLIP_
      LM = size(A,3)

      do L=1,LM/2
         B             = A(:,:,L)
         A(:,:,     L) = A(:,:,LM+1-L)
         A(:,:,LM+1-L) = B
      end do
#endif /* _DOFLIP_ */

      return
    end subroutine Flip3

  end subroutine GSI2GCM

  subroutine Dot_TransIn(VarsGSI, VarsGCM, VarsTRAJ3, varsTRAJ2, TR, TRt, PHASE_NAME, rc)
    type(T_3Dvar),                  intent(INOUT) :: VarsGSI(:)
    type(T_3Dvar),                  intent(INOUT) :: VarsGCM(:)
    type(T_3Dvar), target,          intent(IN   ) :: VarsTRAJ3(:)
    type(T_2Dvar), target,          intent(IN   ) :: VarsTRAJ2(:)
    class (AbstractRegridder), intent(in) :: TR
    class (AbstractRegridder), intent(in) :: TRt
    character(len=ESMF_MAXSTR)                    :: PHASE_NAME
    integer,     optional,          intent(  OUT) :: rc

! ErrLog vars

    integer                    :: STATUS
    character(len=ESMF_MAXSTR) :: IAm="DotTransIn"

! Local vars

    integer                    :: IM, JM, LM
    logical                    :: Transpose
    logical                    :: FillGSI
    logical                    :: Inverse
    character*60               :: GridTypeIn
    character*3                :: Dot_GridType

    integer                            :: i,i_,j,k
    integer                            :: img,jmg,kmg
    type(T_3DVAR)                      :: varsGSI_dot(NUM_GSIvars)
    type(T_3DVAR), allocatable         :: varsGCM_dot(:)

    real, allocatable :: varsGCM_random_(:,:,:,:),varsGSI_random_(:,:,:,:)
    real(R8) :: dot_d1(1),dot_d2(1),d1_gsum,d2_gsum,d1_tmp,d2_tmp
    real(R8) :: dot_d1_uv(1),dot_d2_uv(1),d1_gsum_uv,d2_gsum_uv
    real(R8) :: dot_d1_tv(1),dot_d2_tv(1),d1_gsum_tv,d2_gsum_tv
    real(R8) :: dot_d1_o3(1),dot_d2_o3(1),d1_gsum_o3,d2_gsum_o3
    real(R8) :: dot_d1_qi(1),dot_d2_qi(1),d1_gsum_qi,d2_gsum_qi
    real(R8) :: dot_d1_ql(1),dot_d2_ql(1),d1_gsum_ql,d2_gsum_ql

    allocate(varsGCM_dot(NUM_GCM3Dvars))

    varsGSI_dot    = varsGSI
    varsGCM_dot    = varsGCM

    allocate(VarsGCM_random_(size(varsGCM_dot(1)%x,1),&
                             size(varsGCM_dot(1)%x,2),&
                             size(varsGCM_dot(1)%x,3),&
                             size(varsGCM_dot) ))
    allocate(VarsGSI_random_(size(varsGSI_dot(1)%x,1),&
                             size(varsGSI_dot(1)%x,2),&
                             size(varsGSI_dot(1)%x,3),&
                             size(varsGSI_dot) ))
    varsGCM_random_ = 0.0
    varsGSI_random_ = 0.0

    if (PHASE_NAME.eq.'TLM') then
    do i=1,size(varsGCM_dot)
#ifdef DOT_RAND_PERT
       varsGCM_dot(i)%x = 0.0
       CALL RANDOM_NUMBER( varsGSI_dot(i)%x )
       VarsGSI_random_(:,:,:,i) = varsGSI_dot(i)%x ! save for dot operation later
#else
       VarsGSI_random_(:,:,:,i) = varsGSI(i)%x     ! save for dot operation later
#endif
    enddo
    endif

    if (PHASE_NAME.eq.'ADJ') then
    do i=1,size(varsGSI_dot)
#ifdef DOT_RAND_PERT
       varsGSI_dot(i)%x = 0.0
       CALL RANDOM_NUMBER( varsGCM_dot(i)%x )
       VarsGCM_random_(:,:,:,i) = varsGCM_dot(i)%x ! save for dot operation later
#else
       VarsGCM_random_(:,:,:,i) = varsGCM(i)%x     ! save for dot operation later
#endif
    enddo
    endif

    call GSI2GCM(varsGSI_dot, varsGCM_dot, varsTRAJ3, varsTRAJ2, TR, rc=STATUS)
    VERIFY_(STATUS)

    dot_d1=0
    dot_d1_uv=0
    dot_d1_tv=0
    dot_d1_o3=0
    dot_d1_qi=0
    dot_d1_ql=0

    do i=1,size(varsGSI)
    if (PHASE_NAME.eq.'TLM') then
        kmg=size(varsGCM_dot(i)%X,3)
        jmg=size(varsGCM_dot(i)%X,2)
        img=size(varsGCM_dot(i)%X,1)
    else
        kmg=size(varsGSI_dot(i)%X,3)
        jmg=size(varsGSI_dot(i)%X,2)
        img=size(varsGSI_dot(i)%X,1)
    endif
    if(associated(varsGSI(i)%X)) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d1 = dot_d1 + varsGCM_dot(i)%X(i_,j,k)*varsGCM_dot(i)%X(i_,j,k)
          else
          dot_d1 = dot_d1 + varsGSI_dot(i)%X(i_,j,k)*varsGSI_dot(i)%X(i_,j,k)
          endif
       enddo
       enddo
       enddo
    endif
    if(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."TV".or.&
                                    VarsGSI(i)%name.eq."DP".or.&
                                    VarsGSI(i)%name.eq."QV"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d1_tv = dot_d1_tv + varsGCM_dot(i)%X(i_,j,k)*varsGCM_dot(i)%X(i_,j,k)
          else
          dot_d1_tv = dot_d1_tv + varsGSI_dot(i)%X(i_,j,k)*varsGSI_dot(i)%X(i_,j,k)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."U".or.&
                                        VarsGSI(i)%name.eq."V"    &
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d1_uv = dot_d1_uv + varsGCM_dot(i)%X(i_,j,k)*varsGCM_dot(i)%X(i_,j,k)
          else
          dot_d1_uv = dot_d1_uv + varsGSI_dot(i)%X(i_,j,k)*varsGSI_dot(i)%X(i_,j,k)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."O3"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d1_o3 = dot_d1_o3 + varsGCM_dot(i)%X(i_,j,k)*varsGCM_dot(i)%X(i_,j,k)
          else
          dot_d1_o3 = dot_d1_o3 + varsGSI_dot(i)%X(i_,j,k)*varsGSI_dot(i)%X(i_,j,k)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."QI"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d1_qi = dot_d1_qi + varsGCM_dot(i)%X(i_,j,k)*varsGCM_dot(i)%X(i_,j,k)
          else
          dot_d1_qi = dot_d1_qi + varsGSI_dot(i)%X(i_,j,k)*varsGSI_dot(i)%X(i_,j,k)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."QL"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d1_ql = dot_d1_ql + varsGCM_dot(i)%X(i_,j,k)*varsGCM_dot(i)%X(i_,j,k)
          else
          dot_d1_ql = dot_d1_ql + varsGSI_dot(i)%X(i_,j,k)*varsGSI_dot(i)%X(i_,j,k)
          endif
       enddo
       enddo
       enddo
    endif
    enddo

    do i=1,size(varsGSI_dot)
       if (PHASE_NAME.eq.'TLM') then
           varsGSI_dot(i)%X=0.
       else
           varsGCM_dot(i)%X=0.
       endif
    enddo

    call GSI2GCM(varsGSI_dot, varsGCM_dot, varsTRAJ3, varsTRAJ2, TRt, rc=STATUS)
    VERIFY_(STATUS)

    dot_d2=0
    dot_d2_uv=0
    dot_d2_tv=0
    dot_d2_o3=0
    dot_d2_qi=0
    dot_d2_ql=0

    do i=1,size(varsGSI)
    if (PHASE_NAME.eq.'TLM') then
        kmg=size(varsGSI_dot(i)%X,3)
        jmg=size(varsGSI_dot(i)%X,2)
        img=size(varsGSI_dot(i)%X,1)
    else
        kmg=size(varsGCM_dot(i)%X,3)
        jmg=size(varsGCM_dot(i)%X,2)
        img=size(varsGCM_dot(i)%X,1)
    endif
    if( associated(varsGSI(i)%X) ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d2=dot_d2+ varsGSI_dot(i)%X(i_,j,k)*varsGSI_random_(i_,j,k,i)
          else
          dot_d2=dot_d2+ varsGCM_dot(i)%X(i_,j,k)*varsGCM_random_(i_,j,k,i)
          endif
       enddo
       enddo
       enddo
    endif
    if(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."TV".or.&
                                    VarsGSI(i)%name.eq."DP".or.&
                                    VarsGSI(i)%name.eq."QV"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d2_tv=dot_d2_tv+ varsGSI_dot(i)%X(i_,j,k)*varsGSI_random_(i_,j,k,i)
          else
          dot_d2_tv=dot_d2_tv+ varsGCM_dot(i)%X(i_,j,k)*varsGCM_random_(i_,j,k,i)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."U".or.&
                                        VarsGSI(i)%name.eq."V"    &
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d2_uv=dot_d2_uv+ varsGSI_dot(i)%X(i_,j,k)*varsGSI_random_(i_,j,k,i)
          else
          dot_d2_uv=dot_d2_uv+ varsGCM_dot(i)%X(i_,j,k)*varsGCM_random_(i_,j,k,i)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."O3"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d2_o3=dot_d2_o3+ varsGSI_dot(i)%X(i_,j,k)*varsGSI_random_(i_,j,k,i)
          else
          dot_d2_o3=dot_d2_o3+ varsGCM_dot(i)%X(i_,j,k)*varsGCM_random_(i_,j,k,i)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."QI"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d2_qi=dot_d2_qi+ varsGSI_dot(i)%X(i_,j,k)*varsGSI_random_(i_,j,k,i)
          else
          dot_d2_qi=dot_d2_qi+ varsGCM_dot(i)%X(i_,j,k)*varsGCM_random_(i_,j,k,i)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."QL"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d2_ql=dot_d2_ql+ varsGSI_dot(i)%X(i_,j,k)*varsGSI_random_(i_,j,k,i)
          else
          dot_d2_ql=dot_d2_ql+ varsGCM_dot(i)%X(i_,j,k)*varsGCM_random_(i_,j,k,i)
          endif
       enddo
       enddo
       enddo
    endif
    enddo

    d1_gsum=0.
    d2_gsum=0.
    d1_gsum_tv=0.
    d2_gsum_tv=0.
    d1_gsum_uv=0.
    d2_gsum_uv=0.
    d1_gsum_o3=0.
    d2_gsum_o3=0.
    d1_gsum_qi=0.
    d2_gsum_qi=0.
    d1_gsum_ql=0.
    d2_gsum_ql=0.

    call ESMF_VmGetCurrent(vm)
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d1), &
                                recvbuf=d1_gsum, cnt= 1, rc=status)
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d2), &
                                recvbuf=d2_gsum, cnt= 1, rc=status)
    IF (d2_gsum.ne.0.) then
    IF (MAPL_AM_I_ROOT()) &
    print *, 'dot TransIn  (gsum error):',trim(PHASE_NAME),d1_gsum,(d1_gsum-d2_gsum)/d2_gsum
    ENDIF
!tv
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d1_tv), &
                                recvbuf=d1_gsum_tv, cnt= 1, rc=status)
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d2_tv), &
                                recvbuf=d2_gsum_tv, cnt= 1, rc=status)
    if(d2_gsum_tv.ne.0.) then
    IF (MAPL_AM_I_ROOT()) &
    print *, 'dot product tv ',trim(PHASE_NAME),d1_gsum_tv,(d1_gsum_tv-d2_gsum_tv)/d2_gsum_tv
    endif
!uv
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d1_uv), &
                                recvbuf=d1_gsum_uv, cnt= 1, rc=status)
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d2_uv), &
                                recvbuf=d2_gsum_uv, cnt= 1, rc=status)
    if(d2_gsum_uv.ne.0.) then
    IF (MAPL_AM_I_ROOT()) &
    print *, 'dot product uv ',trim(PHASE_NAME),d1_gsum_uv,(d1_gsum_uv-d2_gsum_uv)/d2_gsum_uv
    endif
!o3
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d1_o3), &
                                recvbuf=d1_gsum_o3, cnt= 1, rc=status)
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d2_o3), &
                                recvbuf=d2_gsum_o3, cnt= 1, rc=status)
    if(d2_gsum_o3.ne.0.) then
    IF (MAPL_AM_I_ROOT()) &
    print *, 'dot product o3 ',trim(PHASE_NAME),d1_gsum_o3,(d1_gsum_o3-d2_gsum_o3)/d2_gsum_o3
    endif
!qi
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d1_qi), &
                                recvbuf=d1_gsum_qi, cnt= 1, rc=status)
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d2_qi), &
                                recvbuf=d2_gsum_qi, cnt= 1, rc=status)
    if(d2_gsum_qi.ne.0.) then
    IF (MAPL_AM_I_ROOT()) &
    print *, 'dot product qi ',trim(PHASE_NAME),d1_gsum_qi,(d1_gsum_qi-d2_gsum_qi)/d2_gsum_qi
    endif
!ql
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d1_ql), &
                                recvbuf=d1_gsum_ql, cnt= 1, rc=status)
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d2_ql), &
                                recvbuf=d2_gsum_ql, cnt= 1, rc=status)
    if(d2_gsum_ql.ne.0.) then
    IF (MAPL_AM_I_ROOT()) &
    print *, 'dot product ql ',trim(PHASE_NAME),d1_gsum_ql,(d1_gsum_ql-d2_gsum_ql)/d2_gsum_ql
    endif

    deallocate(VarsGCM_random_,VarsGSI_random_)  
    deallocate(varsGCM_dot)

  end subroutine Dot_TransIn

  subroutine Dot_TransOut(VarsGSI, VarsGCM, VarsTRAJ3, varsTRAJ2, TR, TRt, PHASE_NAME, rc)
    type(T_3Dvar),                  intent(INOUT) :: VarsGSI(:)
    type(T_3Dvar),                  intent(INOUT) :: VarsGCM(:)
    type(T_3Dvar), target,          intent(IN   ) :: VarsTRAJ3(:)
    type(T_2Dvar), target,          intent(IN   ) :: VarsTRAJ2(:)
    class (AbstractRegridder), intent(in) :: TR
    class (AbstractRegridder), intent(in) :: TRt
    character(len=ESMF_MAXSTR)                    :: PHASE_NAME
    integer,     optional,          intent(  OUT) :: rc

! ErrLog vars

    integer                    :: STATUS
    character(len=ESMF_MAXSTR) :: IAm="DotTransOut"

! Local vars

    integer                    :: IM, JM, LM
    logical                    :: Transpose
    logical                    :: FillGSI
    logical                    :: Inverse
    character*60               :: GridTypeIn
    character*3                :: Dot_GridType

    integer                            :: i,i_,j,k
    type(T_3DVAR)                      :: varsGSI_dot(NUM_GSIvars)
    type(T_3DVAR), allocatable         :: varsGCM_dot(:)

    real, allocatable :: varsGCM_random_(:,:,:,:),varsGSI_random_(:,:,:,:)
    real(R8) :: dot_d1(1),dot_d2(1),d1_gsum,d2_gsum,d1_tmp,d2_tmp
    real(R8) :: dot_d1_uv(1),dot_d2_uv(1),d1_gsum_uv,d2_gsum_uv
    real(R8) :: dot_d1_tv(1),dot_d2_tv(1),d1_gsum_tv,d2_gsum_tv
    real(R8) :: dot_d1_o3(1),dot_d2_o3(1),d1_gsum_o3,d2_gsum_o3
    real(R8) :: dot_d1_qi(1),dot_d2_qi(1),d1_gsum_qi,d2_gsum_qi
    real(R8) :: dot_d1_ql(1),dot_d2_ql(1),d1_gsum_ql,d2_gsum_ql

    integer img,jmg,kmg

    allocate(varsGCM_dot(NUM_GCM3Dvars))

    varsGSI_dot    = varsGSI
    varsGCM_dot    = varsGCM

    allocate(VarsGCM_random_(size(varsGCM_dot(1)%x,1),&
                             size(varsGCM_dot(1)%x,2),&
                             size(varsGCM_dot(1)%x,3),&
                             size(varsGCM_dot) ))
    allocate(VarsGSI_random_(size(varsGSI_dot(1)%x,1),&
                             size(varsGSI_dot(1)%x,2),&
                             size(varsGSI_dot(1)%x,3),&
                             size(varsGSI_dot) ))
    varsGCM_random_ = 0.0
    varsGSI_random_ = 0.0

    if (PHASE_NAME.eq.'TLM') then
    do i=1,size(varsGCM_dot)
#ifdef DOT_RAND_PERT
       CALL RANDOM_NUMBER( varsGCM_dot(i)%x )
       varsGSI_dot(i)%x    = 0.0
       varsGCM_random_(:,:,:,i) = varsGCM_dot(i)%x ! save for dot product later
#else
       varsGCM_random_(:,:,:,i) = varsGCM(i)%x     ! save for dot product later
#endif
    enddo
    endif

    if (PHASE_NAME.eq.'ADJ') then
    do i=1,size(varsGSI_dot)
#ifdef DOT_RAND_PERT
       CALL RANDOM_NUMBER( varsGSI_dot(i)%x )
       varsGCM_dot(i)%x    = 0.0
       varsGSI_random_(:,:,:,i) = varsGSI_dot(i)%x ! save for dot product later
#else
       varsGSI_random_(:,:,:,i) = varsGSI(i)%x      ! save for dot product later
#endif
    enddo
    endif

    call GSI2GCM(varsGSI_dot, varsGCM_dot, varsTRAJ3, varsTRAJ2, TR, rc=STATUS)
    VERIFY_(STATUS)

    dot_d2=0.d0
    dot_d2_uv=0.d0
    dot_d2_tv=0.d0
    dot_d2_o3=0.d0
    dot_d2_qi=0.d0
    dot_d2_ql=0.d0

    do i=1,size(varsGSI)
    if (PHASE_NAME.eq.'TLM') then
       kmg = size(varsGSI_dot(i)%X,3) 
       jmg = size(varsGSI_dot(i)%X,2) 
       img = size(varsGSI_dot(i)%X,1) 
    else
       kmg = size(varsGCM_dot(i)%X,3) 
       jmg = size(varsGCM_dot(i)%X,2) 
       img = size(varsGCM_dot(i)%X,1) 
    endif
    if( associated(varsGSI(i)%X) ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d2 = dot_d2 + varsGSI_dot(i)%X(i_,j,k)*varsGSI_dot(i)%X(i_,j,k)
          else
          dot_d2 = dot_d2 + varsGCM_dot(i)%X(i_,j,k)*varsGCM_dot(i)%X(i_,j,k)
          endif
       enddo
       enddo
       enddo
    endif
    if(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."TV".or.&
                                    VarsGSI(i)%name.eq."DP".or.&
                                    VarsGSI(i)%name.eq."QV"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d2_tv = dot_d2_tv + varsGSI_dot(i)%X(i_,j,k)*varsGSI_dot(i)%X(i_,j,k)
          else
          dot_d2_tv = dot_d2_tv + varsGCM_dot(i)%X(i_,j,k)*varsGCM_dot(i)%X(i_,j,k)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."U".or.&
                                        VarsGSI(i)%name.eq."V"    &
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d2_uv = dot_d2_uv + varsGSI_dot(i)%X(i_,j,k)*varsGSI_dot(i)%X(i_,j,k)
          else
          dot_d2_uv = dot_d2_uv + varsGCM_dot(i)%X(i_,j,k)*varsGCM_dot(i)%X(i_,j,k)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."O3"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d2_o3 = dot_d2_o3 + varsGSI_dot(i)%X(i_,j,k)*varsGSI_dot(i)%X(i_,j,k)
          else
          dot_d2_o3 = dot_d2_o3 + varsGCM_dot(i)%X(i_,j,k)*varsGCM_dot(i)%X(i_,j,k)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."QI"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d2_qi = dot_d2_qi + varsGSI_dot(i)%X(i_,j,k)*varsGSI_dot(i)%X(i_,j,k)
          else
          dot_d2_qi = dot_d2_qi + varsGCM_dot(i)%X(i_,j,k)*varsGCM_dot(i)%X(i_,j,k)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."QL"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d2_ql = dot_d2_ql + varsGSI_dot(i)%X(i_,j,k)*varsGSI_dot(i)%X(i_,j,k)
          else
          dot_d2_ql = dot_d2_ql + varsGCM_dot(i)%X(i_,j,k)*varsGCM_dot(i)%X(i_,j,k)
          endif
       enddo
       enddo
       enddo
    endif
    enddo

    do i=1,size(varsGSI_dot)
       if (PHASE_NAME.eq.'TLM') then
           varsGCM_dot(i)%X=0.
       else
           varsGSI_dot(i)%X=0.
       endif
    enddo

    call GSI2GCM(varsGSI_dot, varsGCM_dot, varsTRAJ3, varsTRAJ2, TRt, rc=STATUS)
    VERIFY_(STATUS)

    dot_d1=0.d0
    dot_d1_uv=0.d0
    dot_d1_tv=0.d0
    dot_d1_o3=0.d0
    dot_d1_qi=0.d0
    dot_d1_ql=0.d0

    do i=1,size(varsGSI)
    if (PHASE_NAME.eq.'TLM') then
       kmg = size(varsGCM_dot(i)%X,3) 
       jmg = size(varsGCM_dot(i)%X,2) 
       img = size(varsGCM_dot(i)%X,1) 
    else
       kmg = size(varsGSI_dot(i)%X,3) 
       jmg = size(varsGSI_dot(i)%X,2) 
       img = size(varsGSI_dot(i)%X,1) 
    endif
    if(associated(varsGSI(i)%X)) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d1 = dot_d1 + varsGCM_dot(i)%X(i_,j,k)*varsGCM_random_(i_,j,k,i)
          else
          dot_d1 = dot_d1 + varsGSI_dot(i)%X(i_,j,k)*varsGSI_random_(i_,j,k,i)
          endif
       enddo
       enddo
       enddo
    endif
    if(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."TV".or.&
                                    VarsGSI(i)%name.eq."DP".or.&
                                    VarsGSI(i)%name.eq."QV"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d1_tv=dot_d1_tv+ varsGCM_dot(i)%X(i_,j,k)*varsGCM_random_(i_,j,k,i)
          else
          dot_d1_tv=dot_d1_tv+ varsGSI_dot(i)%X(i_,j,k)*varsGSI_random_(i_,j,k,i)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."U".or.&
                                        VarsGSI(i)%name.eq."V"    &
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d1_uv=dot_d1_uv+ varsGCM_dot(i)%X(i_,j,k)*varsGCM_random_(i_,j,k,i)
          else
          dot_d1_uv=dot_d1_uv+ varsGSI_dot(i)%X(i_,j,k)*varsGSI_random_(i_,j,k,i)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."O3"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d1_o3=dot_d1_o3+ varsGCM_dot(i)%X(i_,j,k)*varsGCM_random_(i_,j,k,i)
          else
          dot_d1_o3=dot_d1_o3+ varsGSI_dot(i)%X(i_,j,k)*varsGSI_random_(i_,j,k,i)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."QI"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d1_qi=dot_d1_qi+ varsGCM_dot(i)%X(i_,j,k)*varsGCM_random_(i_,j,k,i)
          else
          dot_d1_qi=dot_d1_qi+ varsGSI_dot(i)%X(i_,j,k)*varsGSI_random_(i_,j,k,i)
          endif
       enddo
       enddo
       enddo
    elseif(associated(varsGSI(i)%X).and.VarsGSI(i)%name.eq."QL"&
    ) then
       do k =1,kmg
       do j =1,jmg
       do i_=1,img
          if (PHASE_NAME.eq.'TLM') then
          dot_d1_ql=dot_d1_ql+ varsGCM_dot(i)%X(i_,j,k)*varsGCM_random_(i_,j,k,i)
          else
          dot_d1_ql=dot_d1_ql+ varsGSI_dot(i)%X(i_,j,k)*varsGSI_random_(i_,j,k,i)
          endif
       enddo
       enddo
       enddo
    endif
    enddo

    d1_gsum=0.d0
    d2_gsum=0.d0
    d1_gsum_tv=0.d0
    d2_gsum_tv=0.d0
    d1_gsum_uv=0.d0
    d2_gsum_uv=0.d0
    d1_gsum_o3=0.d0
    d2_gsum_o3=0.d0
    d1_gsum_qi=0.d0
    d2_gsum_qi=0.d0
    d1_gsum_ql=0.d0
    d2_gsum_ql=0.d0

    call ESMF_VmGetCurrent(vm)
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d1), &
                                recvbuf=d1_gsum, cnt= 1, rc=status)
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d2), &
                                recvbuf=d2_gsum, cnt= 1, rc=status)
    IF (d2_gsum>0.d0) then
    IF (MAPL_AM_I_ROOT()) &
    print *, 'dot TransOut (gsum error):',trim(PHASE_NAME),d1_gsum,(d1_gsum-d2_gsum)/d2_gsum
    ENDIF
    IF (MAPL_AM_I_ROOT()) print *, ' '

!tv
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d1_tv), &
                                recvbuf=d1_gsum_tv, cnt= 1, rc=status)
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d2_tv), &
                                recvbuf=d2_gsum_tv, cnt= 1, rc=status)
    if(d2_gsum_tv>0.d0) then
    IF (MAPL_AM_I_ROOT()) &
    print *, 'dot product tv ',trim(PHASE_NAME),d1_gsum_tv,(d1_gsum_tv-d2_gsum_tv)/d2_gsum_tv
    endif
!uv
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d1_uv), &
                                recvbuf=d1_gsum_uv, cnt= 1, rc=status)
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d2_uv), &
                                recvbuf=d2_gsum_uv, cnt= 1, rc=status)
    if(d2_gsum_uv>0.d0) then
    IF (MAPL_AM_I_ROOT()) &
    print *, 'dot product uv ',trim(PHASE_NAME),d1_gsum_uv,(d1_gsum_uv-d2_gsum_uv)/d2_gsum_uv
    endif
!o3
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d1_o3), &
                                recvbuf=d1_gsum_o3, cnt= 1, rc=status)
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d2_o3), &
                                recvbuf=d2_gsum_o3, cnt= 1, rc=status)
    if(d2_gsum_o3>0.d0) then
    IF (MAPL_AM_I_ROOT()) &
    print *, 'dot product o3 ',trim(PHASE_NAME),d1_gsum_o3,(d1_gsum_o3-d2_gsum_o3)/d2_gsum_o3
    endif
!qi
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d1_qi), &
                                recvbuf=d1_gsum_qi, cnt= 1, rc=status)
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d2_qi), &
                                recvbuf=d2_gsum_qi, cnt= 1, rc=status)
    if(d2_gsum_qi>0.d0) then
    IF (MAPL_AM_I_ROOT()) &
    print *, 'dot product qi ',trim(PHASE_NAME),d1_gsum_qi,(d1_gsum_qi-d2_gsum_qi)/d2_gsum_qi
    endif
!ql
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d1_ql), &
                                recvbuf=d1_gsum_ql, cnt= 1, rc=status)
    call MAPL_CommsAllReduceSum(vm, sendbuf=sum(dot_d2_ql), &
                                recvbuf=d2_gsum_ql, cnt= 1, rc=status)
    if(d2_gsum_ql>0.d0) then
    IF (MAPL_AM_I_ROOT()) &
    print *, 'dot product ql ',trim(PHASE_NAME),d1_gsum_ql,(d1_gsum_ql-d2_gsum_ql)/d2_gsum_ql
    endif
    IF (MAPL_AM_I_ROOT()) print *, ' '
    deallocate(VarsGCM_random_,VarsGSI_random_)  
    deallocate(varsGCM_dot)

  end subroutine Dot_TransOut

  subroutine resolve_template_ ( fname, ftmpl, TIME, msg, status )
  use m_StrTemplate
  implicit none
    character(len=ESMF_MAXSTR),intent(out) :: fname
    character(len=ESMF_MAXSTR),intent(in)  :: ftmpl
    type(ESMF_Time)                        :: Time
    character(len=*),          intent(in)  :: msg
    integer,intent(out) :: status
    
    integer :: YY, MM, DD, HH, MN, SC, nymd, nhms
    integer :: rc

    call ESMF_TimeGet(Time, yy=YY, mm=MM, dd=DD, h=HH, m=MN, s=SC, rc=rc)
    if(MAPL_AM_I_ROOT()) write(6,'(a,1x,i4.4,4(a,i2.2))') trim(msg)//" The current TIME is ", &
                                        YY, "/", MM, "/", DD, " ", HH, ":", MN
    nymd = yy*10000 + mm*100 + dd
    nhms = hh*10000 + mn*100 + sc
    call StrTemplate ( fname, ftmpl, 'GRADS', nymd=nymd, nhms=nhms,stat=STATUS )
  end subroutine resolve_template_

  subroutine gq_ (CF,ntrj,wght)
    implicit none

    type(ESMF_Config), intent(inout) :: CF
    integer,           intent(  out) :: ntrj
    real,              intent(  out) :: wght

    character(len=ESMF_MAXSTR), parameter :: IAm="gq_"

    character(len=ESMF_MAXSTR) :: charbuf
    integer :: STATUS, RC

    ! initialize
    ntrj = 0
    wght = 0.0

    ! go to the desired label in resource file
    CALL ESMF_ConfigFindLabel(CF, label = 'TRAJ_FILE::', __RC__ )
    do
      CALL ESMF_ConfigNextLine    (CF,          __RC__)
      CALL ESMF_ConfigGetAttribute(CF, charbuf, __RC__)
      if ( trim(charbuf) == '::' ) exit
      ntrj = ntrj + 1
    end do

    ! get the wght's for quadrature based on number of points
    select case ( ntrj )
      case (1)
        wght = 1.0
      case (2)
        wght = 0.5
    end select

    if ( MAPL_am_I_ROOT() ) then
       print*,'In subroutine gq_ ntrj = ', ntrj
       print*,'In subroutine gq_ wght = ', wght
    endif

    return
  end subroutine gq_

  subroutine get_TRAJ_FILE_TMPL_ (CF, ntrj, TMPL)
    implicit none

    type(ESMF_Config),           intent(inout) :: CF
    integer,                     intent(in   ) :: ntrj
    character(len=ESMF_MAXSTR),  intent(inout) :: TMPL(ntrj)

    character(len=ESMF_MAXSTR), parameter :: IAm = "get_TRAJ_FILE_TMPL_"
    
    integer :: i, STATUS, RC
    character(len=ESMF_MAXSTR) :: charbuf
    
    CALL ESMF_ConfigFindLabel(CF, label = 'TRAJ_FILE::', __RC__)
    do i=1,ntrj
      CALL ESMF_ConfigNextLine    (CF,          __RC__)
      CALL ESMF_ConfigGetAttribute(CF, charbuf, __RC__)
      TMPL(i) = trim(adjustl(charbuf))
    enddo

    return
  end subroutine get_TRAJ_FILE_TMPL_

  subroutine sax_Bundle_ (Bundle, alpha)
    implicit none

    type(ESMF_FieldBundle), intent(inout) :: Bundle
    real,                   intent(in  )  :: alpha

    type(ESMF_Field) :: Field
    integer :: nq, rank
    integer :: i, j, STATUS, RC
    real, pointer, dimension(:,:)   :: X2 => null()
    real, pointer, dimension(:,:,:) :: X3 => null()

    character(len=ESMF_MAXSTR), parameter :: IAm="sax_Bundle_"

    call ESMF_FieldBundleGet(Bundle, fieldCount=nq, __RC__)

    do i=1,nq

      call ESMF_FieldBundleGet(Bundle, i, Field, __RC__)
      call ESMF_FieldGet(Field, dimCount=rank,   __RC__)

      if ( rank == 2 ) then
        call ESMF_FieldGet(Field, farrayPtr=X2, __RC__)
        X2 = alpha * X2
      elseif ( rank == 3 ) then
        call ESMF_FieldGet(Field, farrayPtr=X3, __RC__)
        X3 = alpha * X3
      endif

    enddo

    return
  end subroutine sax_Bundle_

  subroutine saxpy_Bundle_ (BundleX, BundleY, alpha)
    implicit none

    type(ESMF_FieldBundle), intent(inout) :: BundleX, BundleY
    real,                   intent(in   ) :: alpha

    type(ESMF_Field) :: FieldX, FieldY
    integer :: nqX, nqY, rankX, rankY
    character(len=ESMF_MAXSTR) :: nameX, nameY
    integer :: i, j, STATUS, RC
    real, pointer, dimension(:,:)   :: X2 => null()
    real, pointer, dimension(:,:)   :: Y2 => null()
    real, pointer, dimension(:,:,:) :: X3 => null()
    real, pointer, dimension(:,:,:) :: Y3 => null()

    character(len=ESMF_MAXSTR), parameter :: IAm="saxpy_Bundle_"

    call ESMF_FieldBundleGet(BundleX, fieldCount=nqX, __RC__)
    call ESMF_FieldBundleGet(BundleY, fieldCount=nqY, __RC__)

    if ( nqX /= nqY ) STATUS = 100
    VERIFY_(STATUS)

    do i=1,nqX

      call ESMF_FieldBundleGet(BundleX, i, FieldX,           __RC__)
      call ESMF_FieldGet(FieldX, name=nameX, dimCount=rankX, __RC__)

      call ESMF_FieldBundleGet(BundleY, fieldName=trim(nameX), field=FieldY, __RC__)
      call ESMF_FieldGet(FieldY, dimCount=rankY,             __RC__)

      if ( rankX /= rankY ) STATUS = 200
      VERIFY_(STATUS)

      if ( rankX == 2 ) then
        call ESMF_FieldGet(FieldX, farrayPtr=X2, __RC__)
        call ESMF_FieldGet(FieldY, farrayPtr=Y2, __RC__)
        X2 = X2 + alpha * Y2
      elseif ( rankX == 3 ) then
        call ESMF_FieldGet(FieldX, farrayPtr=X3, __RC__)
        call ESMF_FieldGet(FieldY, farrayPtr=Y3, __RC__)
        X3 = X3 + alpha * Y3
      endif

    enddo

    return
  end subroutine saxpy_Bundle_

  subroutine echo_time0_ ( CLOCK, msg )
   implicit none
   type(ESMF_Clock) :: Clock
   character(len=*) :: msg

   type(ESMF_Time)  :: Time
   integer :: YY, MM, DD, HH, MN, SC, nymd, nhms
   integer :: rc,status
   character(len=*), parameter:: Iam ='echo_time0_'

   call ESMF_ClockGet(Clock, CurrTime=Time, __RC__)
   call echo_time1_( TIME, msg )
  end subroutine echo_time0_

  subroutine echo_time1_ ( TIME, msg )
  use m_StrTemplate
  implicit none
    type(ESMF_Time)  :: Time
    character(len=*) :: msg
    
    integer :: YY, MM, DD, HH, MN, SC, nymd, nhms
    integer :: rc

    call ESMF_TimeGet(Time, yy=YY, mm=MM, dd=DD, h=HH, m=MN, s=SC, rc=rc)
    if(MAPL_AM_I_ROOT()) write(6,'(a,1x,i4.4,4(a,i2.2))') trim(msg)//" The current TIME is ", &
                                        YY, "/", MM, "/", DD, " ", HH, ":", MN
  end subroutine echo_time1_

  subroutine Read_Traj(TRAJ_wrk3_, TRAJ_wrk2_, Time, IAm)

  integer                            :: STATUS
  integer                            :: rc  ! return code

  type(ESMF_FieldBundle)             :: TRAJ_wrk3_, TRAJ_wrk2_
  type (ESMF_Time)                   :: Time
  character(len=ESMF_MAXSTR)         :: lstvars
  character(len=ESMF_MAXSTR)         :: TRAJ_FILE
  character(len=ESMF_MAXSTR)         :: IAm

  real, pointer                      :: PT(:,:,:)
  integer                            :: j

  call Resolve_template_ ( TRAJ_FILE, TRAJ_FILE_TMPL(1), Time, trim(Iam), STATUS ); VERIFY_(STATUS)

  if (DO_MOIST_PHYS == 0) then
     lstvars = vec2list(GCM3DvarsRead)
  else
     lstvars = vec2list(GCM3DvarsReadMoist)
  endif
  call MAPL_CFIORead( TRAJ_FILE, Time, Traj_wrk3_, only_vars=lstvars, rc=STATUS)
  lstvars = vec2list(GCM2Dvars(1:NUM_GCM2Dvars))
  call MAPL_CFIORead( TRAJ_FILE, Time, Traj_wrk2_, only_vars=lstvars, rc=STATUS)

  if ( gq_ntrj > 1) then
     call sax_Bundle_(Traj_wrk3_, gq_wght)
     call sax_Bundle_(Traj_wrk2_, gq_wght)
     do j=2,gq_ntrj
        call Resolve_template_ ( TRAJ_FILE, TRAJ_FILE_TMPL(J), Time, trim(Iam), STATUS ); VERIFY_(STATUS)
        if (DO_MOIST_PHYS == 0) then
           lstvars = vec2list(GCM3DvarsRead)
        else
           lstvars = vec2list(GCM3DvarsReadMoist)
        endif
        call MAPL_CFIORead(TRAJ_FILE, Time, TRAJwrk3_4gq, only_vars=lstvars, __RC__)
        lstvars = vec2list(GCM2Dvars(1:NUM_GCM2Dvars))
        call MAPL_CFIORead(TRAJ_FILE, Time, TRAJwrk2_4gq, only_vars=lstvars, __RC__)

        call saxpy_Bundle_(Traj_wrk3_, TRAJwrk3_4gq, gq_wght)
        call saxpy_Bundle_(Traj_wrk2_, TRAJwrk2_4gq, gq_wght)
     enddo

  endif

  end subroutine Read_Traj

  subroutine Get_Traj_Size(n, n_skip, m_out)
  integer n_remain, n, m, m_out, l, n_skip

  if (n_skip > 1) then
     m = int((n-1)/n_skip)
     if (m.le.0) then
        if(MAPL_AM_I_ROOT()) write(*,*) "Traj_freq is larger than total integration time."
        call exit(0)
     endif
     l = n_skip*m+1
     if (l.eq.n) then
        m_out=m+1
     else
        n_remain=n-l
        m_out=m +n -m*n_skip
        if(MAPL_AM_I_ROOT()) write(*,*) "Traj_freq does not fit with GEOS_AgcmPert job segment."
        call exit(0)
     endif
  else
     m_out=n
  endif
  end subroutine Get_Traj_Size

  subroutine Get_Traj_Idx(n, n_skip, m_in)
  integer, intent(in) :: n, n_skip, m_in 
  integer n_remain, m, l, i, j, icount

  if (n_skip > 1) then
      m = int((n-1)/n_skip)
      l = n_skip*m+1
      
      if (l.eq.n) then
      j=1
      icount=1
      Traj_Skip_Idx(icount)= j
      do i=1,m
         j=j+n_skip
         icount=icount+1
         Traj_Skip_Idx(icount)= j
      enddo
      endif
      if (l.ne.n) then
      n_remain = n-l
      j = 1
      icount = 1
      Traj_Skip_Idx(icount)= j
      do i=1,m
         j=j+n_skip
         icount=icount+1
         Traj_Skip_Idx(icount)= j
      enddo
      do i=j+1,n
         icount=icount+1
         Traj_Skip_Idx(icount)= i
      enddo
      endif
      do i=2,m_in
         do j=Traj_Skip_Idx(i-1),Traj_Skip_Idx(i)-1
            Traj_Orgi_Idx(j)=i-1
         enddo
      enddo
      Traj_Orgi_Idx(n) = m_in
  else
      do i=1,n
         Traj_Orgi_Idx(i)= i
         Traj_Skip_Idx(i)= i
      enddo
  endif
  end subroutine Get_Traj_Idx

end module GEOS_AgcmPertGridCompMod
