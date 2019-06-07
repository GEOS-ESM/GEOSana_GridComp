!  $Id$

! VERIFY_ and RETURN_ macros for error handling.

#include "MAPL_Generic.h"

#define MAPL_FieldBundleGetPointer ESMFL_BundleGetPointerToData

!=============================================================================
!BOP

! !MODULE: GEOS_RadiationPert -- A Module to compute perturbation radiation processes. 
!          Includes linearized longwave (irrad) and shortwave (solar) radiation.
!
!          Fluxes are updated only every hour (DT=3600s), in the last time step of the hour. 
!          Those fluxes are stored in an internal state and applied to the temperature at 
!          every time step.
!
!          Info on TRACE and OVERCAST used in IRRAD
!          If trace = .true., absorption due to n2o, ch4, cfcs, and the 
!          two minor co2 bands in the window region is included.
!          Otherwise, absorption in those minor bands is neglected.
!
!          If overcast=.true. , the layer cloud cover is either 0 or 1.
!          If overcast=.false., the cloud cover can be anywhere between 0 and 1.
!          Computation is faster for the .true. option than the .false. option.
!
!          List of simplifications
!          - No rain and snow trajecotry or perturbations
!          - No N2O, CH4, CFC trajectory or perturbations
!          - No pressure perturbation 



! !INTERFACE:

module GEOS_RadiationPertGridCompMod

! ESMF and MAPL components
  use ESMF
  use MAPL_MOD
  use GEOS_Mod
  use GEOS_PertSharedMod
  use AeroOptPropTableMod

! Radiation constants
  use rad_constants
  use irrad_constants
  use sorad_constants

! Module to couple cloud variables
  use cloudradcouple

! Main radiation TL/AD subroutines
  use IRRAD_TL
  use IRRAD_AD
  use SORAD_TL
  use SORAD_AD

! Timers from core
  use fv_timing_mod,  only: timing_on, timing_off

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION:
! 
!   {\tt GEOS_RadiationPertGridCompMod} implements linearised radiative processes in GEOS-5. 
!   

!EOP

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code
!EOP


!=============================================================================
!
! ErrLog Variables

    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: IAm
    character(len=ESMF_MAXSTR)              :: COMP_NAME
    
!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,   name="RUNTLM"       ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="RUNADJ"       ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="-RADIATION"   ,RC=STATUS)
    VERIFY_(STATUS)

! The same run method is registered twice
!  It queries the phase to do TLM or ADJ
! ---------------------------------------

    call MAPL_GridCompSetEntryPoint (gc, ESMF_METHOD_RUN,  Run,        rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint (gc, ESMF_METHOD_RUN,  Run,        rc=status)
    VERIFY_(STATUS)

! Tangent linear internal state
! -----------------------------

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'IRFLX_TL',                                  &
        LONG_NAME  = 'pert_net_downward_longwave_flux_in_air_tl', &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                                                       &
        SHORT_NAME = 'IRDFDTS_TL',                                                                      &
        LONG_NAME  = 'pert_sensitivity_of_net_downward_longwave_flux_in_air_to_surface_temperature_tl', &
        UNITS      = 'W m-2',                                                                           &
        DIMS       = MAPL_DimsHorzVert,                                                                 &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS                                        )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'SOFLX_TL',                                  &
        LONG_NAME  = 'pert_net_downward_shortwave_flux_in_air_tl',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

! Adjoint internal state
! ----------------------

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'IRFLX_AD',                                  &
        LONG_NAME  = 'pert_net_downward_longwave_flux_in_air_ad', &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                                                       &
        SHORT_NAME = 'IRDFDTS_AD',                                                                      &
        LONG_NAME  = 'pert_sensitivity_of_net_downward_longwave_flux_in_air_to_surface_temperature_ad', &
        UNITS      = 'W m-2',                                                                           &
        DIMS       = MAPL_DimsHorzVert,                                                                 &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS                                        )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                 &
        SHORT_NAME = 'SOFLX_AD',                                  &
        LONG_NAME  = 'pert_net_downward_shortwave_flux_in_air_ad',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)


! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)
     
! Return
! ------
    RETURN_(ESMF_SUCCESS)
     
  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Run(gc, import, export, clock, rc)

  !ARGUMENTS
  type(ESMF_GridComp), intent(inout) :: gc
  type (ESMF_State),   intent(inout) :: import
  type (ESMF_State),   intent(inout) :: export
  type (ESMF_Clock),   intent(inout) :: clock
  integer,  optional,  intent(  out) :: rc 

  !MODEL SETUP VARIABLES
  type (MAPL_MetaComp), pointer      :: MAPL
  type (ESMF_Config)                 :: CF
  type (ESMF_State)                  :: INTERNAL
  integer                            :: STATUS
  character(len=ESMF_MAXSTR)         :: IAm
  character(len=ESMF_MAXSTR)         :: COMP_NAME
  character(len=ESMF_MAXSTR)         :: PHASE_NAME
  integer                            :: PHASE
  real(8)                            :: DT
  type (ESMF_FieldBundle)            :: Traj3, Traj2, Pert
  type (ESMF_Time)                   :: CURRENTTIME
  integer                            :: DOY, YY, MM, DD, HH, MN, SEC

  !LOCAL VARIABLES
  integer                            :: RUNalarm
  integer                            :: IM, JM, LM, I, J, L, AI

  !POINTERS TO THE REFERENCE AND PERTURBATION TRAJECTORIES.
  real*4,  pointer, dimension(:,:,:)  :: UT4, PTT4, QVT4, O3T4, CFLST4, CFCNT4
  real*4,  pointer, dimension(:,:,:)  ::      PTP4, QVP4, O3P4, CFLSP4, CFCNP4
  real*4,  pointer, dimension(:,:,:)  :: DU001T4, DU002T4, DU003T4, DU004T4, DU005T4
  real*4,  pointer, dimension(:,:,:)  :: DU001P4, DU002P4, DU003P4, DU004P4, DU005P4
  real*4,  pointer, dimension(:,:,:)  :: QLST4, QCNT4
  real*4,  pointer, dimension(:,:,:)  :: QIP4, QLP4
  real*4,  pointer, dimension(:,:)    :: PS4, TS4
  real*4,  pointer, dimension(:,:)    :: EMIS4, DELT4, COSZ4, SLR4
  real*4,  pointer, dimension(:,:)    :: RGBUV4, RGFUV4, RGBIR4, RGFIR4

  real(8), pointer, dimension(:,:,:)  :: UT, PTT, QVT, O3T, CFLST, CFCNT
  real(8), pointer, dimension(:,:,:)  ::     PTP, QVP, O3P, CFLSP, CFCNP
  real(8), pointer, dimension(:,:,:)  :: DU001T, DU002T, DU003T, DU004T, DU005T
  real(8), pointer, dimension(:,:,:)  :: DU001P, DU002P, DU003P, DU004P, DU005P
  real(8), pointer, dimension(:,:,:)  :: QLST, QCNT
  real(8), pointer, dimension(:,:,:)  :: QIP, QLP
  real(8), pointer, dimension(:,:)    :: PS, TS
  real(8), pointer, dimension(:,:)    :: EMIS, DELT, COSZ, SLR
  real(8), pointer, dimension(:,:)    :: RGBUV, RGFUV, RGBIR, RGFIR

  !LOCAL REFERENCE AND PERTURBATION TRAJECTORY VARIABLES
  real(8), allocatable, dimension(:,:,:) :: TT, PLE, PLESO, PLOf
  real(8), allocatable, dimension(:,:,:) :: TP
  real(8), allocatable, dimension(:,:)   :: T2MT, T2MP
  real(8), allocatable, dimension(:,:,:) :: DTDtP, DELTAP
  real(8), allocatable, dimension(:,:,:) :: PKT, PET, PKE, BX, CX

  !REFERENCE PRESSURE, AK & BK
  real(8), allocatable, dimension(:)  :: VERT_AK(:), VERT_BK(:)
  real(8), allocatable, dimension(:)  :: PREF

  !ATMOSPHERIC GASES
  real(8), allocatable, dimension(:,:,:) :: O3mmT, O3mmP
  real*4                                 :: CO2_4
  real(8)                                :: CO2
  real(8), allocatable, dimension(:,:,:) :: N2OT,CH4T,CFC11T,CFC12T,CFC22T

  !CLOUD VARIABLES
  real(8), allocatable, dimension(:,:,:) :: fQi, ILSF, ICNF, LLSF, LCNF
  real(8), allocatable, dimension(:,:,:) :: QILST, QLLST, QICNT, QLCNT, QIT, QLT
  real(8), allocatable, dimension(:,:,:) :: QILSP, QLLSP, QICNP, QLCNP


  integer :: levs925, LCLDMH, LCLDLM
  real(8), allocatable, dimension(:,:)   :: TEMPOR
  real(8), allocatable, dimension(:,:,:) :: RAD_QLT, RAD_QIT, RAD_CFT, RAD_RLT, RAD_RIT
  real(8), allocatable, dimension(:,:,:) :: RAD_QLP, RAD_QIP, RAD_CFP, RAD_RLP, RAD_RIP

  real(8), allocatable, dimension(:,:,:,:) :: CWCT, REFFT
  real(8), allocatable, dimension(:,:,:,:) :: CWCP, REFFP

  !SURFACE
  integer, parameter                       :: NS = 1
  real(8), allocatable, dimension(:,:,:)   :: FST, TGT, TVT
  real(8), allocatable, dimension(:,:,:,:) :: EGT, EVT, RVT

  !AEROSOLS
  integer, parameter                       :: NA = 0
  character(len=ESMF_MAXSTR), pointer      :: AEROSOLS(:)
  real(8)                                  :: X
  integer                                  :: OFFSET
  real(8), allocatable, dimension(:,:,:,:) :: AEROT, AEROP
  real(8), allocatable, dimension(:,:,:,:) :: SAEROT, SAEROP

  !AEROSOLS - MieTable
  real*4,  allocatable, dimension(:,:,:)     :: RH
  real*4,  allocatable, dimension(:,:,:,:,:) :: TAUA_IR_C, SSAA_IR_C, ASYA_IR_C
  real*4,  allocatable, dimension(:,:,:,:,:) :: TAUA_SO_C, SSAA_SO_C, ASYA_SO_C

  !AEROSOLS - IRRAD
  integer, parameter                         :: NBCHOU_IR = 10
  real(8), allocatable, dimension(:,:,:,:)   :: TAUA_IRT, SSAA_IRT, ASYA_IRT
  real(8), allocatable, dimension(:,:,:,:)   :: TAUA_IRP, SSAA_IRP, ASYA_IRP

  !AEROSOLS - SORAD
  integer, parameter                         :: NBCHOU_SO = 8
  real(8), allocatable, dimension(:,:,:,:)   :: TAUA_SOT, SSAA_SOT, ASYA_SOT
  real(8), allocatable, dimension(:,:,:,:)   :: TAUA_SOP, SSAA_SOP, ASYA_SOP

  !IRRAD CONSTANTS
  logical, parameter :: TRACE = .false., OVERCAST = .false. 
  
  !SOLAR VARIABLES
  real*4  :: SC, HK_4(8)
  real(8) :: HK(8), HK_IR_TEMP(3,10), HK_UV_TEMP(5)

  real(8), allocatable, dimension(:,:) :: COSZTMP

  character(len=ESMF_MAXPATHLEN) :: SolCycFileName

! Fluxes
! ------
  !LOCAL IRRAD FLUXES OUTPUT
  real(8), allocatable, dimension(:,:,:) :: IR_UFLXT, IR_DFLXT,           IR_DFDTsT
  real(8), allocatable, dimension(:,:,:) :: IR_UFLXP, IR_DFLXP, IR_NFLXP, IR_DFDTsP
  real(8), allocatable, dimension(:,:,:) :: IR_FLXP

  !LOCAL SORAD FLUXES OUTPUT
  real(8), allocatable, dimension(:,:,:) :: SO_NFLXT
  real(8), allocatable, dimension(:,:,:) :: SO_NFLXP
  real(8), allocatable, dimension(:,:,:) :: SO_FLXP

  !INTERNAL FLUXES FOR LINEAR MODEL
  real*4,  pointer,     dimension(:,:,:) :: IR_NFLXP_TL_INT, SO_NFLXP_TL_INT, IR_DFDTsP_TL_INT
  real*4,  pointer,     dimension(:,:,:) :: IR_NFLXP_AD_INT, SO_NFLXP_AD_INT, IR_DFDTsP_AD_INT
  
! Constants
! ---------
  !Real 8 constants  
  real(8) :: MAPL8_GRAV, MAPL8_P00, MAPL8_RGAS, MAPL8_CP, MAPL8_KAPPA

  !Real 8 constants - Irrad
  real(8)                              :: w11_8, w12_8, w13_8, p11_8, p12_8
  real(8)                              :: p13_8, dwe_8, dpe_8
  real(8), allocatable, dimension(:,:) :: aib_ir_8, awb_ir_8, aiw_ir_8
  real(8), allocatable, dimension(:,:) :: aww_ir_8, aig_ir_8, awg_ir_8
  real(8), allocatable, dimension(:)   :: xkw_8, xke_8, aw_8
  real(8), allocatable, dimension(:)   :: bw_8, pm_8
  real(8), allocatable, dimension(:,:) :: fkw_8, gkw_8, cb_8, dcb_8
  real(8), allocatable, dimension(:,:) :: c1_8, c2_8, c3_8
  real(8), allocatable, dimension(:,:) :: oo1_8, oo2_8, oo3_8
  real(8), allocatable, dimension(:,:) :: h11_8, h12_8, h13_8
  real(8), allocatable, dimension(:,:) :: h21_8, h22_8, h23_8
  real(8), allocatable, dimension(:,:) :: h81_8, h82_8, h83_8

  !Real 8 constants - Sorad
  real(8)                                :: aib_nir_8, aib_uv_8
  real(8), allocatable, dimension(:)     :: wk_uv_8, zk_uv_8, ry_uv_8
  real(8), allocatable, dimension(:)     :: xk_ir_8, ry_ir_8
  real(8), allocatable, dimension(:,:)   :: cah_8, coa_8
  real(8), allocatable, dimension(:)     :: aig_uv_8, awg_uv_8, arg_uv_8
  real(8), allocatable, dimension(:)     :: awb_uv_8, arb_uv_8
  real(8), allocatable, dimension(:,:)   :: awb_nir_8, arb_nir_8
  real(8), allocatable, dimension(:,:)   :: aia_nir_8, awa_nir_8, ara_nir_8
  real(8), allocatable, dimension(:,:)   :: aig_nir_8, awg_nir_8, arg_nir_8
  real(8), allocatable, dimension(:,:,:) :: caib_8
  real(8), allocatable, dimension(:,:)   :: caif_8

  !Debubbing/printing
  type (ESMF_VM) :: vm
  integer :: proc_id

! Begin
! -----

! Debugging options
! -----------------

    call ESMF_VMGetCurrent(VM, rc=status)
    call ESMF_VmGet(VM, localPet=proc_id, rc=status)

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet( GC, name=COMP_NAME, CONFIG=CF, currentPhase=Phase, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the generic state
! -----------------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

! Turn on MAPL timers for rad
! ---------------------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"RADIATION")
    call MAPL_TimerOn(MAPL,PHASE_NAME)

! Begin calculations if we are doing radiation
! --------------------------------------------

    !ONLY DO RADIATION PHYSICS IF SELECTED IN AGCM_apert
    IF (DO_RAD_PHYS .ne. 0) THEN !Do radiation

       !Turn on TLM/ADM timers
       if (PHASE==ADJphase) then         
          call timing_on('RADIA_ADM')             
       elseif (PHASE==TLMphase) then         
          call timing_on('RADIA_TLM')            
       endif

! Get the model time step
! -----------------------

       call MAPL_GetResource(MAPL, DT, Label="RUN_DT:", RC=STATUS)
       VERIFY_(STATUS)

! Get parameters from generic state. The RUNALARM is used to control
!  the calling of the full transfer calculation
! ------------------------------------------------------------------

       call MAPL_Get(MAPL, IM                  = IM,                                     &
                           JM                  = JM,                                     &
                           LM                  = LM,                                     &
                           INTERNAL_ESMF_STATE = INTERNAL,                               &
                     RC=STATUS )    

! Get the current model time
! --------------------------

       call ESMF_ClockGet(CLOCK, currTIME=CURRENTTIME,       RC=STATUS)
       VERIFY_(STATUS)

!Get the Solar Constant
! ---------------------

       call MAPL_GetResource( MAPL, SC, 'SOLAR_CONSTANT:', RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetResource( MAPL, SolCycFileName, "SOLAR_CYCLE_FILE_NAME:", DEFAULT='/dev/null', RC=STATUS)
       VERIFY_(STATUS) 

! Pointers to Internals, tl/ad separated to avoid any dot product test issue 
! --------------------------------------------------------------------------

      call MAPL_GetPointer(INTERNAL, IR_NFLXP_TL_INT , 'IRFLX_TL'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, SO_NFLXP_TL_INT , 'SOFLX_TL'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, IR_DFDTsP_TL_INT, 'IRDFDTS_TL', RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(INTERNAL, IR_NFLXP_AD_INT , 'IRFLX_AD'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, SO_NFLXP_AD_INT , 'SOFLX_AD'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL, IR_DFDTsP_AD_INT, 'IRDFDTS_AD', RC=STATUS); VERIFY_(STATUS)

! Set TLM/ADM phase
! -----------------

       select case(phase)
       case(TLMphase)
          PHASE_NAME = "RUNTLM"
       case(ADJphase)
          PHASE_NAME = "RUNADJ"
       case default
          ASSERT_(.false.)
       end select
    
! Get pointer to reference and perturbation trajectory
! ----------------------------------------------------

       call ESMF_StateGet(Import, "TRAJWRK3", Traj3, rc=status)
       VERIFY_(STATUS)
       call ESMF_StateGet(Import, 'TRAJWRK2', Traj2, rc=STATUS)
       VERIFY_(STATUS)
       call ESMF_StateGet(Import, "PERTWRK" , Pert , rc=STATUS)
       VERIFY_(STATUS)
 
! Allocate memory for local arrays used at every time step
! --------------------------------------------------------

       !Real 8 holders of ref and pert trajectories
       allocate(PTT(IM,JM,LM),                 stat=status); VERIFY_(STATUS)
       allocate(PTP(IM,JM,LM),                 stat=status); VERIFY_(STATUS)
       allocate(PS(IM,JM),                     stat=status); VERIFY_(STATUS)
       allocate(SLR(IM,JM),                    stat=status); VERIFY_(STATUS)
       allocate(DELT(IM,JM),                   stat=status); VERIFY_(STATUS)

       !Local trajectory
       allocate(TT(IM,JM,LM),                  stat=status); VERIFY_(STATUS)
       allocate(TP(IM,JM,LM),                  stat=status); VERIFY_(STATUS)
       allocate(PLE(IM,JM,0:LM),               stat=status); VERIFY_(STATUS)
       allocate(PLOf(IM,JM,LM),                stat=status); VERIFY_(STATUS)
       allocate(DELTAP(IM,JM,LM),              stat=status); VERIFY_(STATUS)
       allocate(VERT_AK(size(pert_ak)),        stat=status); VERIFY_(STATUS)
       allocate(VERT_BK(size(pert_bk)),        stat=status); VERIFY_(STATUS)
       allocate(PREF(0:size(VERT_AK)-1),       stat=status); VERIFY_(STATUS)
       allocate(DTDtP(IM,JM,LM),               stat=status); VERIFY_(STATUS)
       allocate(PKT(IM,JM,  LM),               stat=status); VERIFY_(STATUS)
       allocate(PET(IM,JM,0:LM),               stat=status); VERIFY_(STATUS)
       allocate(BX (IM,JM,  LM),               stat=status); VERIFY_(STATUS)
       allocate(CX (IM,JM,  LM),               stat=status); VERIFY_(STATUS)
       allocate(PKE(IM,JM,0:LM),               stat=status); VERIFY_(STATUS)

       !IR fluxes
       allocate(IR_UFLXT(IM,JM,0:LM),          stat=status); VERIFY_(STATUS)
       allocate(IR_UFLXP(IM,JM,0:LM),          stat=status); VERIFY_(STATUS)
       allocate(IR_DFLXT(IM,JM,0:LM),          stat=status); VERIFY_(STATUS)
       allocate(IR_DFLXP(IM,JM,0:LM),          stat=status); VERIFY_(STATUS)
       allocate(IR_NFLXP(IM,JM,0:LM),          stat=status); VERIFY_(STATUS)
       allocate(IR_DFDTsT(IM,JM,0:LM),         stat=status); VERIFY_(STATUS)
       allocate(IR_DFDTsP(IM,JM,0:LM),         stat=status); VERIFY_(STATUS)
       allocate(IR_FLXP(IM,JM,0:LM),           stat=status); VERIFY_(STATUS)

       !SO fluxes
       allocate(SO_NFLXT(IM,JM,0:LM),          stat=status); VERIFY_(STATUS)
       allocate(SO_NFLXP(IM,JM,0:LM),          stat=status); VERIFY_(STATUS)
       allocate(SO_FLXP(IM,JM,0:LM),           stat=status); VERIFY_(STATUS)

! Get the trajectory variables needed even if the alarm is not ringing
! --------------------------------------------------------------------

       call MAPL_FieldBundleGetPointer(Traj3, "PT"    , PTT4    , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Pert , "PT"    , PTP4    , rc=STATUS)
       VERIFY_(STATUS)

       call MAPL_FieldBundleGetPointer(Traj2, "PS"    , PS4     , rc=STATUS)
       VERIFY_(STATUS)

       call MAPL_FieldBundleGetPointer(Traj2, "SLR"   , SLR4    , rc=STATUS)
       VERIFY_(STATUS)  
  
       call MAPL_FieldBundleGetPointer(Traj2, "DELTIRD"  , DELT4   , rc=STATUS)
       VERIFY_(STATUS)

       !Move to real 8
       PTT  = dble(PTT4)
       PTP  = dble(PTP4)
       PS   = dble(PS4)
       SLR  = dble(SLR4)
       DELT = dble(DELT4)

! Trajectory calculations done at every time step
! -------------------------------------------------

       !Real 8 versions of MAPL_Constants
       MAPL8_GRAV  = dble(MAPL_GRAV)
       MAPL8_P00   = dble(MAPL_P00)
       MAPL8_RGAS  = dble(MAPL_RGAS)
       MAPL8_CP    = dble(MAPL_CP)
       MAPL8_KAPPA = dble(MAPL_KAPPA)

       !ak and bk pressure variables
       VERT_AK = dble(pert_ak)
       VERT_BK = dble(pert_bk)

       !Reference pressure caclulation (Pa)
       DO L = 0,LM
          PREF(L) = VERT_AK(L+1) + VERT_BK(L+1)*MAPL8_P00
       enddo

       !Pressure at the half levels from Ps (Pa)
       DO L = 0,LM
          PLE(:,:,L) = VERT_AK(L+1) + VERT_BK(L+1)*PS(:,:)
       enddo
 
       !Pressure thickness (Pa)
       DELTAP = ( PLE(:,:,1:LM) - PLE(:,:,0:LM-1) )

       !Pressure (hPa) at the full levels
       PLOf(:,:,1:LM) = 0.01 * 0.5 * (PLE(:,:,0:LM-1) +  PLE(:,:,1:LM  ) )

       !Get PKT, P to the power kappa
       call GetPressVarsTrajR8(DELTAP, PET, PKE, PKT, BX, CX)

       !FV scaled potential temperature to air temperature (K)
       TT = PTT*PKT


! Pre flux update perturbation calculations done at every time step
! -------------------------------------------------------------------

       if (PHASE==TLMphase) then 

          !FV scaled potential temperature to air temperature
          TP = PTP*PKT

       elseif (PHASE==ADJphase) then

          !Adjoint of temperature to FV scaled potential temperture
          TP = PTP/PKT

          !Adjoint of unweighted temperature tendency calculation
          DTDtP = DT*TP/DELTAP

          !Adjoint of weighted temperature tendency calculation
          IR_FLXP = 0.0
          IR_FLXP(:,:,0:LM-1) = IR_FLXP(:,:,0:LM-1) + (MAPL8_GRAV*DTDtP/MAPL8_CP)
          IR_FLXP(:,:,1:LM)   = IR_FLXP(:,:,1:LM)   - (MAPL8_GRAV*DTDtP/MAPL8_CP)
          SO_FLXP = 0.0
          SO_FLXP(:,:,0:LM-1) = SO_FLXP(:,:,0:LM-1) + (MAPL8_GRAV*DTDtP/MAPL8_CP)
          SO_FLXP(:,:,1:LM)   = SO_FLXP(:,:,1:LM)   - (MAPL8_GRAV*DTDtP/MAPL8_CP)

          !Adjoint of final total flux calculation
          IR_NFLXP  = dble(IR_NFLXP_AD_INT)
          IR_DFDTsP = dble(IR_DFDTsP_AD_INT)
          SO_NFLXP  = dble(SO_NFLXP_AD_INT)

          !Adjoint of final total flux calculation - IR part
          DO L = LM,0,-1
             IR_NFLXP(:,:,L) = IR_NFLXP(:,:,L) + IR_FLXP(:,:,L)
             IR_DFDTsP(:,:,L) = IR_DFDTsP(:,:,L) + DELT*IR_FLXP(:,:,L)
             IR_FLXP(:,:,L) = 0.0
          END DO

          !Adjoint of final total flux caclulation - SO part
          DO L = LM,0,-1
             SO_NFLXP(:,:,L) = SO_NFLXP(:,:,L) + SLR*SO_FLXP(:,:,L)
             SO_FLXP(:,:,L) = 0.0
          END DO

          !Place back into internal state to be added at next time
          IR_NFLXP_AD_INT = real(IR_NFLXP,4)
          IR_DFDTsP_AD_INT = real(IR_DFDTsP,4)
          SO_NFLXP_AD_INT = real(SO_NFLXP,4)

       endif

! Determine if it is time to re-do linear flux caclulations
! ---------------------------------------------------------

       RunAlarm = 0

       !Get the modle time
       call ESMF_TimeGet (CurrentTime, YY=YY, MM=MM, DD=DD, H=HH, M=MN, S=SEC, RC=STATUS )
       VERIFY_(STATUS)

       !Turn on switch if we are in the time step immediately after the hour.
       IF (PHASE == TLMPhase) THEN

          if (MN == 0) then
             RunAlarm = 1
          endif

       ELSEIF (PHASE == ADJPhase) THEN

          !If going backwards the reported time is that at `end` of time step
          if ( MN*60 + SEC == nint(DT) ) then
             RunAlarm = 1
          endif

       ENDIF

       if (RUNALARM == 1) then

! Allocate memory for local variables used during flux calulation
! ---------------------------------------------------------------

          !Real 8 holders for pert and traj
          allocate(UT(IM,JM,LM),                      stat=status); VERIFY_(STATUS)
          allocate(QVT(IM,JM,LM),                     stat=status); VERIFY_(STATUS)
          allocate(O3T(IM,JM,LM),                     stat=status); VERIFY_(STATUS)
          allocate(CFLST(IM,JM,LM),                   stat=status); VERIFY_(STATUS)
          allocate(CFCNT(IM,JM,LM),                   stat=status); VERIFY_(STATUS)
          allocate(QVP(IM,JM,LM),                     stat=status); VERIFY_(STATUS)
          allocate(O3P(IM,JM,LM),                     stat=status); VERIFY_(STATUS)
          allocate(CFLSP(IM,JM,LM),                   stat=status); VERIFY_(STATUS)
          allocate(CFCNP(IM,JM,LM),                   stat=status); VERIFY_(STATUS)
          allocate(QLST(IM,JM,LM),                    stat=status); VERIFY_(STATUS)
          allocate(QCNT(IM,JM,LM),                    stat=status); VERIFY_(STATUS)
          allocate(QIP(IM,JM,LM),                     stat=status); VERIFY_(STATUS)
          allocate(QLP(IM,JM,LM),                     stat=status); VERIFY_(STATUS)
          allocate(TS(IM,JM),                         stat=status); VERIFY_(STATUS)
          allocate(EMIS(IM,JM),                       stat=status); VERIFY_(STATUS)
          allocate(COSZ(IM,JM),                       stat=status); VERIFY_(STATUS)
          allocate(RGBUV(IM,JM),                      stat=status); VERIFY_(STATUS)
          allocate(RGFUV(IM,JM),                      stat=status); VERIFY_(STATUS)
          allocate(RGBIR(IM,JM),                      stat=status); VERIFY_(STATUS)
          allocate(RGFIR(IM,JM),                      stat=status); VERIFY_(STATUS)

          !Temperatures
          allocate(T2MT(IM,JM),                       stat=status); VERIFY_(STATUS)
          allocate(T2MP(IM,JM),                       stat=status); VERIFY_(STATUS)
          allocate(PLESO(IM,JM,LM+1),                 stat=status); VERIFY_(STATUS)
          allocate(COSZTMP(IM,JM),                    stat=status); VERIFY_(STATUS)

          !Atmospheric gases
          allocate(O3mmT(IM,JM,LM),                   stat=status); VERIFY_(STATUS)
          allocate(O3mmP(IM,JM,LM),                   stat=status); VERIFY_(STATUS)
          allocate(N2OT(IM,JM,LM),                    stat=status); VERIFY_(STATUS)
          allocate(CH4T(IM,JM,LM),                    stat=status); VERIFY_(STATUS)
          allocate(CFC11T(IM,JM,LM),                  stat=status); VERIFY_(STATUS)
          allocate(CFC12T(IM,JM,LM),                  stat=status); VERIFY_(STATUS)
          allocate(CFC22T(IM,JM,LM),                  stat=status); VERIFY_(STATUS)

          !Cloud variables
          allocate(QILST(IM,JM,LM),                   stat=status); VERIFY_(STATUS)
          allocate(QLLST(IM,JM,LM),                   stat=status); VERIFY_(STATUS)
          allocate(QICNT(IM,JM,LM),                   stat=status); VERIFY_(STATUS)
          allocate(QLCNT(IM,JM,LM),                   stat=status); VERIFY_(STATUS)
          allocate(QILSP(IM,JM,LM),                   stat=status); VERIFY_(STATUS)
          allocate(QLLSP(IM,JM,LM),                   stat=status); VERIFY_(STATUS)
          allocate(QICNP(IM,JM,LM),                   stat=status); VERIFY_(STATUS)
          allocate(QLCNP(IM,JM,LM),                   stat=status); VERIFY_(STATUS)
          allocate(QIT(IM,JM,LM),                     stat=status); VERIFY_(STATUS)
          allocate(QLT(IM,JM,LM),                     stat=status); VERIFY_(STATUS)
          allocate(fQi(IM,JM,LM),                     stat=status); VERIFY_(STATUS)
          allocate(ILSF(IM,JM,LM),                    stat=status); VERIFY_(STATUS)
          allocate(ICNF(IM,JM,LM),                    stat=status); VERIFY_(STATUS)
          allocate(LLSF(IM,JM,LM),                    stat=status); VERIFY_(STATUS)
          allocate(LCNF(IM,JM,LM),                    stat=status); VERIFY_(STATUS)
          allocate(CWCT(IM,JM,LM,4),                  stat=status); VERIFY_(STATUS)
          allocate(CWCP(IM,JM,LM,4),                  stat=status); VERIFY_(STATUS)
          allocate(REFFT(IM,JM,LM,4),                 stat=status); VERIFY_(STATUS)
          allocate(REFFP(IM,JM,LM,4),                 stat=status); VERIFY_(STATUS)
          allocate(TEMPOR(IM,JM),                     stat=status); VERIFY_(STATUS)
          allocate(RAD_QLT(IM,JM,LM),                 stat=status); VERIFY_(STATUS)
          allocate(RAD_QLP(IM,JM,LM),                 stat=status); VERIFY_(STATUS)
          allocate(RAD_QIT(IM,JM,LM),                 stat=status); VERIFY_(STATUS)
          allocate(RAD_QIP(IM,JM,LM),                 stat=status); VERIFY_(STATUS)
          allocate(RAD_CFT(IM,JM,LM),                 stat=status); VERIFY_(STATUS)
          allocate(RAD_CFP(IM,JM,LM),                 stat=status); VERIFY_(STATUS)
          allocate(RAD_RLT(IM,JM,LM),                 stat=status); VERIFY_(STATUS)
          allocate(RAD_RLP(IM,JM,LM),                 stat=status); VERIFY_(STATUS)
          allocate(RAD_RIT(IM,JM,LM),                 stat=status); VERIFY_(STATUS)
          allocate(RAD_RIP(IM,JM,LM),                 stat=status); VERIFY_(STATUS)
 
          !Surface radiative properties
          allocate(FST(IM,JM,ns),                     stat=status); VERIFY_(STATUS)
          allocate(TGT(IM,JM,ns),                     stat=status); VERIFY_(STATUS)
          allocate(TVT(IM,JM,ns),                     stat=status); VERIFY_(STATUS)
          allocate(EGT(IM,JM,ns,10),                  stat=status); VERIFY_(STATUS)
          allocate(EVT(IM,JM,ns,10),                  stat=status); VERIFY_(STATUS)
          allocate(RVT(IM,JM,ns,10),                  stat=status); VERIFY_(STATUS)

          allocate(aib_ir_8(size(aib_ir,1),size(aib_ir,2)),        stat=status); VERIFY_(STATUS)
          allocate(awb_ir_8(size(awb_ir,1),size(awb_ir,2)),        stat=status); VERIFY_(STATUS)
          allocate(aiw_ir_8(size(aiw_ir,1),size(aiw_ir,2)),        stat=status); VERIFY_(STATUS)
          allocate(aww_ir_8(size(aww_ir,1),size(aww_ir,2)),        stat=status); VERIFY_(STATUS)
          allocate(aig_ir_8(size(aig_ir,1),size(aig_ir,2)),        stat=status); VERIFY_(STATUS)
          allocate(awg_ir_8(size(awg_ir,1),size(awg_ir,2)),        stat=status); VERIFY_(STATUS)
          allocate(xkw_8(size(xkw,1)),                             stat=status); VERIFY_(STATUS)
          allocate(xke_8(size(xke,1)),                             stat=status); VERIFY_(STATUS)
          allocate(aw_8(size(aw,1)),                               stat=status); VERIFY_(STATUS)
          allocate(bw_8(size(bw,1)),                               stat=status); VERIFY_(STATUS)
          allocate(pm_8(size(pm,1)),                               stat=status); VERIFY_(STATUS)
          allocate(fkw_8(size(fkw,1),size(fkw,2)),                 stat=status); VERIFY_(STATUS)
          allocate(gkw_8(size(gkw,1),size(gkw,2)),                 stat=status); VERIFY_(STATUS)
          allocate(cb_8(size(cb,1),size(cb,2)),                    stat=status); VERIFY_(STATUS)
          allocate(dcb_8(size(dcb,1),size(dcb,2)),                 stat=status); VERIFY_(STATUS)
          allocate(c1_8(size(c1,1),size(c1,2)),                    stat=status); VERIFY_(STATUS)
          allocate(c2_8(size(c2,1),size(c2,2)),                    stat=status); VERIFY_(STATUS)
          allocate(c3_8(size(c3,1),size(c3,2)),                    stat=status); VERIFY_(STATUS)
          allocate(oo1_8(size(oo1,1),size(oo1,2)),                 stat=status); VERIFY_(STATUS)
          allocate(oo2_8(size(oo2,1),size(oo2,2)),                 stat=status); VERIFY_(STATUS)
          allocate(oo3_8(size(oo3,1),size(oo3,2)),                 stat=status); VERIFY_(STATUS)
          allocate(h11_8(size(h11,1),size(h11,2)),                 stat=status); VERIFY_(STATUS)
          allocate(h12_8(size(h12,1),size(h12,2)),                 stat=status); VERIFY_(STATUS)
          allocate(h13_8(size(h13,1),size(h13,2)),                 stat=status); VERIFY_(STATUS)
          allocate(h21_8(size(h21,1),size(h21,2)),                 stat=status); VERIFY_(STATUS)
          allocate(h22_8(size(h22,1),size(h22,2)),                 stat=status); VERIFY_(STATUS)
          allocate(h23_8(size(h23,1),size(h23,2)),                 stat=status); VERIFY_(STATUS)
          allocate(h81_8(size(h81,1),size(h81,2)),                 stat=status); VERIFY_(STATUS)
          allocate(h82_8(size(h82,1),size(h82,2)),                 stat=status); VERIFY_(STATUS)
          allocate(h83_8(size(h83,1),size(h83,2)),                 stat=status); VERIFY_(STATUS)

          allocate(wk_uv_8(size(wk_uv,1)),                         stat=status); VERIFY_(STATUS)
          allocate(zk_uv_8(size(zk_uv,1)),                         stat=status); VERIFY_(STATUS)
          allocate(ry_uv_8(size(ry_uv,1)),                         stat=status); VERIFY_(STATUS)
          allocate(xk_ir_8(size(xk_ir,1)),                         stat=status); VERIFY_(STATUS)
          allocate(ry_ir_8(size(ry_ir,1)),                         stat=status); VERIFY_(STATUS)
          allocate(cah_8(size(cah,1),size(cah,2)),                 stat=status); VERIFY_(STATUS)
          allocate(coa_8(size(coa,1),size(coa,2)),                 stat=status); VERIFY_(STATUS)
          allocate(aig_uv_8(size(aig_uv,1)),                       stat=status); VERIFY_(STATUS)
          allocate(awg_uv_8(size(awg_uv,1)),                       stat=status); VERIFY_(STATUS)
          allocate(arg_uv_8(size(arg_uv,1)),                       stat=status); VERIFY_(STATUS)
          allocate(awb_uv_8(size(awb_uv,1)),                       stat=status); VERIFY_(STATUS)
          allocate(arb_uv_8(size(arb_uv,1)),                       stat=status); VERIFY_(STATUS)
          allocate(awb_nir_8(size(awb_nir,1),size(awb_nir,2)),     stat=status); VERIFY_(STATUS)
          allocate(arb_nir_8(size(arb_nir,1),size(arb_nir,2)),     stat=status); VERIFY_(STATUS)
          allocate(aia_nir_8(size(aia_nir,1),size(aia_nir,2)),     stat=status); VERIFY_(STATUS)
          allocate(awa_nir_8(size(awa_nir,1),size(awa_nir,2)),     stat=status); VERIFY_(STATUS)
          allocate(ara_nir_8(size(ara_nir,1),size(ara_nir,2)),     stat=status); VERIFY_(STATUS)
          allocate(aig_nir_8(size(aig_nir,1),size(aig_nir,2)),     stat=status); VERIFY_(STATUS)
          allocate(awg_nir_8(size(awg_nir,1),size(awg_nir,2)),     stat=status); VERIFY_(STATUS)
          allocate(arg_nir_8(size(arg_nir,1),size(arg_nir,2)),     stat=status); VERIFY_(STATUS)
          allocate(caib_8(size(caib,1),size(caib,2),size(caib,3)), stat=status); VERIFY_(STATUS)
          allocate(caif_8(size(caif,1),size(caif,2)),              stat=status); VERIFY_(STATUS)

          allocate(TAUA_SOT(IM,JM,LM,nbchou_so),      stat=status); VERIFY_(STATUS)
          allocate(TAUA_SOP(IM,JM,LM,nbchou_so),      stat=status); VERIFY_(STATUS)
          allocate(SSAA_SOT(IM,JM,LM,nbchou_so),      stat=status); VERIFY_(STATUS)
          allocate(SSAA_SOP(IM,JM,LM,nbchou_so),      stat=status); VERIFY_(STATUS)
          allocate(ASYA_SOT(IM,JM,LM,nbchou_so),      stat=status); VERIFY_(STATUS)
          allocate(ASYA_SOP(IM,JM,LM,nbchou_so),      stat=status); VERIFY_(STATUS)
          allocate(TAUA_IRT(IM,JM,LM,nbchou_ir),      stat=status); VERIFY_(STATUS)
          allocate(TAUA_IRP(IM,JM,LM,nbchou_ir),      stat=status); VERIFY_(STATUS)
          allocate(SSAA_IRT(IM,JM,LM,nbchou_ir),      stat=status); VERIFY_(STATUS)
          allocate(SSAA_IRP(IM,JM,LM,nbchou_ir),      stat=status); VERIFY_(STATUS)
          allocate(ASYA_IRT(IM,JM,LM,nbchou_ir),      stat=status); VERIFY_(STATUS)
          allocate(ASYA_IRP(IM,JM,LM,nbchou_ir),      stat=status); VERIFY_(STATUS)

          if (NA > 0) then !Doing aerosols so allocate their space
             allocate(DU001T(IM,JM,LM),                  stat=status); VERIFY_(STATUS)
             allocate(DU002T(IM,JM,LM),                  stat=status); VERIFY_(STATUS)
             allocate(DU003T(IM,JM,LM),                  stat=status); VERIFY_(STATUS)
             allocate(DU004T(IM,JM,LM),                  stat=status); VERIFY_(STATUS)
             allocate(DU005T(IM,JM,LM),                  stat=status); VERIFY_(STATUS)
             allocate(DU001P(IM,JM,LM),                  stat=status); VERIFY_(STATUS)
             allocate(DU002P(IM,JM,LM),                  stat=status); VERIFY_(STATUS)
             allocate(DU003P(IM,JM,LM),                  stat=status); VERIFY_(STATUS)
             allocate(DU004P(IM,JM,LM),                  stat=status); VERIFY_(STATUS)
             allocate(DU005P(IM,JM,LM),                  stat=status); VERIFY_(STATUS)

             allocate(AEROSOLS(na),                      stat=status); VERIFY_(STATUS)
             allocate(AEROT(IM,JM,LM,na),                stat=status); VERIFY_(STATUS)
             allocate(AEROP(IM,JM,LM,na),                stat=status); VERIFY_(STATUS)
             allocate(SAEROT(IM,JM,LM,na),               stat=status); VERIFY_(STATUS)
             allocate(SAEROP(IM,JM,LM,na),               stat=status); VERIFY_(STATUS)
             allocate(RH(IM,JM,LM),                      stat=status); VERIFY_(STATUS)

             allocate(TAUA_IR_C(IM,JM,LM,nbchou_ir,na),  stat=status); VERIFY_(STATUS)
             allocate(SSAA_IR_C(IM,JM,LM,nbchou_ir,na),  stat=status); VERIFY_(STATUS)
             allocate(ASYA_IR_C(IM,JM,LM,nbchou_ir,na),  stat=status); VERIFY_(STATUS)

             allocate(TAUA_SO_C(IM,JM,LM,nbchou_ir,na),  stat=status); VERIFY_(STATUS)
             allocate(SSAA_SO_C(IM,JM,LM,nbchou_ir,na),  stat=status); VERIFY_(STATUS)
             allocate(ASYA_SO_C(IM,JM,LM,nbchou_ir,na),  stat=status); VERIFY_(STATUS)

          endif

! Get the rest of the trajectory used in flux calculation
! -------------------------------------------------------

          !U wind (T)
          call MAPL_FieldBundleGetPointer(Traj3, "U"     , UT4     , rc=STATUS)
          VERIFY_(STATUS)
      
          !Specific humidity (T & P)
          call MAPL_FieldBundleGetPointer(Traj3, "QV"    , QVT4    , rc=STATUS)
          VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Pert , "QV"    , QVP4    , rc=STATUS)
          VERIFY_(STATUS)

          !Ozone (T & P)
          call MAPL_FieldBundleGetPointer(Traj3, "O3"    , O3T4    , rc=STATUS)
          VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Pert , "O3"    , O3P4    , rc=STATUS)
          VERIFY_(STATUS)

          !Total cloud liquid ice and water (P)
          call MAPL_FieldBundleGetPointer(Pert , "QI"    , QIP4    , rc=STATUS)
          VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Pert , "QL"    , QLP4    , rc=STATUS)
          VERIFY_(STATUS)

          !Cloud totals for large scale and convective (T)
          call MAPL_FieldBundleGetPointer(Traj3, "QLS"   , QLST4   , rc=STATUS)
          VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Traj3, "QCN"   , QCNT4   , rc=STATUS)
          VERIFY_(STATUS)

          !Anvil fraction (T & P)
          call MAPL_FieldBundleGetPointer(Traj3, "CFCN"  , CFCNT4  , rc=STATUS)
          VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Pert , "CFCN"  , CFCNP4  , rc=STATUS)
          VERIFY_(STATUS)

          !Cloud fraction (T & P) 
          call MAPL_FieldBundleGetPointer(Traj3, "CFLS"  , CFLST4  , rc=STATUS)
          VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Pert , "CFLS"  , CFLSP4  , rc=STATUS)
          VERIFY_(STATUS)

          !Here we read the aerosols if there are any
          if (NA > 0) then
             !Aerosols
             call MAPL_FieldBundleGetPointer(Traj3, "DU001"  , DU001T4  , rc=STATUS)
             VERIFY_(STATUS)
             call MAPL_FieldBundleGetPointer(Pert , "DU001"  , DU001P4  , rc=STATUS)
             VERIFY_(STATUS)

             call MAPL_FieldBundleGetPointer(Traj3, "DU002"  , DU002T4  , rc=STATUS)
             VERIFY_(STATUS)
             call MAPL_FieldBundleGetPointer(Pert , "DU002"  , DU002P4  , rc=STATUS)
             VERIFY_(STATUS)

             call MAPL_FieldBundleGetPointer(Traj3, "DU003"  , DU003T4  , rc=STATUS)
             VERIFY_(STATUS)
             call MAPL_FieldBundleGetPointer(Pert , "DU003"  , DU003P4  , rc=STATUS)
             VERIFY_(STATUS)

             call MAPL_FieldBundleGetPointer(Traj3, "DU004"  , DU004T4  , rc=STATUS)
             VERIFY_(STATUS)
             call MAPL_FieldBundleGetPointer(Pert , "DU004"  , DU004P4  , rc=STATUS)
             VERIFY_(STATUS)

             call MAPL_FieldBundleGetPointer(Traj3, "DU005"  , DU005T4  , rc=STATUS)
             VERIFY_(STATUS)
             call MAPL_FieldBundleGetPointer(Pert , "DU005"  , DU005P4  , rc=STATUS)
             VERIFY_(STATUS)
          endif

          !Surface variables (T)
          call MAPL_FieldBundleGetPointer(Traj2, "TS"    , TS4     , rc=STATUS)
          VERIFY_(STATUS)

          !Radiation variables
          call MAPL_FieldBundleGetPointer(Traj2, "EMIS"     , EMIS4   , rc=STATUS)
          VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Traj2, "RGBUV"    , RGBUV4  , rc=STATUS)
          VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Traj2, "RGFUV"    , RGFUV4  , rc=STATUS)
          VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Traj2, "RGBIR"    , RGBIR4  , rc=STATUS)
          VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Traj2, "RGFIR"    , RGFIR4  , rc=STATUS)
          VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Traj2, "COSZ"     , COSZ4   , rc=STATUS)
          VERIFY_(STATUS)

          !Move to real 8 containers
          UT     = dble(UT4)
          QVT    = dble(QVT4)
          QVP    = dble(QVP4)
          O3T    = dble(O3T4)
          O3P    = dble(O3P4)
          QIP    = dble(QIP4)
          QLP    = dble(QLP4)
          QLST   = dble(QLST4)
          QCNT   = dble(QCNT4)
          CFCNT  = dble(CFCNT4)
          CFCNP  = dble(CFCNP4)
          CFLST  = dble(CFLST4)
          CFLSP  = dble(CFLSP4)
          TS     = dble(TS4)
          EMIS   = dble(EMIS4)
          RGBUV  = dble(RGBUV4)
          RGFUV  = dble(RGFUV4)
          RGBIR  = dble(RGBIR4)
          RGFIR  = dble(RGFIR4)
          COSZ   = dble(COSZ4)

          if (NA > 0) then !Doing aerosols
             DU001T = dble(DU001T4)
             DU001P = dble(DU001P4)
             DU002T = dble(DU002T4)
             DU002P = dble(DU002P4)
             DU003T = dble(DU003T4)
             DU003P = dble(DU003P4)
             DU004T = dble(DU004T4)
             DU004P = dble(DU004P4)
             DU005T = dble(DU005T4)
             DU005P = dble(DU005P4)
          endif

! Trajectory calulations down when alarm in ringing
! -------------------------------------------------

          !2m temperature
          T2MT = TT(:,:,LM)*(0.5*(1.0 + PLE(:,:,LM-1)/PLE(:,:,LM)))**(-MAPL8_KAPPA)

          !Pressure in mb for SORAD
          PLESO(:,:,1:LM+1) = 0.01 * PLE(:,:,0:LM)

          !Caclulate TEMPOR for RAD-Cloud coupling
          levs925  = max(1,count(PREF < 92500.))
          TEMPOR = 0.
          do L = levs925, LM
             where (UT(:,:,L).gt.4.) TEMPOR(:,:) = 1.
          end do

          !Convert ozone from ppmv to kgkg 
          O3mmT = O3T / 1.e6

          !Grab CO2 from the table
          call MAPL_GetResource( MAPL, CO2_4, 'CO2:', RC=STATUS)
          VERIFY_(STATUS)

          if (CO2_4 < 0.0) then
             call ESMF_TimeGet (CURRENTTIME, YY=YY, DayOfYear=DOY, RC=STATUS)
             VERIFY_(STATUS)
             CO2_4 = GETCO2Pert(YY,DOY)
          endif
          CO2 = dble(CO2_4)

          !Initialize other gases, not generally available in the linear model.
          N2OT   = 0.0
          CH4T   = 0.0
          CFC11T = 0.0
          CFC12T = 0.0
          CFC22T = 0.0

          !Split the large and anvil cloud into liquid water and liquid ice parts
          DO I = 1,IM
             DO J = 1,JM
                DO L = 1,LM
                   call IceFraction( TT(I,J,L), fQi(I,J,L) )
                enddo
             enddo
          enddo
          QILST = QLST * fQi
          QLLST = QLST * (1-fQi)
          QICNT = QCNT * fQi
          QLCNT = QCNT * (1-fQi)

          !Fraction of Largescale to convective cloud in reference.
          ILSF = 0.0
          ICNF = 0.0
          LLSF = 0.0
          LCNF = 0.0
          DO I = 1,IM
             DO J = 1,JM
                DO L = 1,LM
                   if ( QILST(I,J,L) + QICNT(I,J,L) .gt. 0.0 ) then
                      ILSF(I,J,L) = QILST(I,J,L) / ( QILST(I,J,L) + QICNT(I,J,L) )
                      ICNF(I,J,L) = QICNT(I,J,L) / ( QILST(I,J,L) + QICNT(I,J,L) )
                   endif
                   if ( QLLST(I,J,L) + QLCNT(I,J,L) .gt. 0.0 ) then
                      LLSF(I,J,L) = QLLST(I,J,L) / ( QLLST(I,J,L) + QLCNT(I,J,L) )
                      LCNF(I,J,L) = QLCNT(I,J,L) / ( QLLST(I,J,L) + QLCNT(I,J,L) )
                   endif
                enddo
             enddo
          enddo

          !Initialize RAD/CLOUD variables
          RAD_CFT = 0.0
          RAD_QLT = 0.0
          RAD_QIT = 0.0
          RAD_RLT = 0.0
          RAD_RIT = 0.0

          !Call RADcouple to produce cloud-rad variables
          !This is normally done in cloudnew in the moist component.
          DO I = 1,IM
             DO J = 1,JM
                DO L = 1,LM
                   call RADCOUPLE( TT(I,J,L), PLOf(I,J,L),                                   & 
                                   CFLST(I,J,L), CFCNT(I,J,L),                               &
                                   QLLST(I,J,L), QILST(I,J,L), QLCNT(I,J,L), QICNT(I,J,L),   & 
                                   RAD_QLT(I,J,L), RAD_QIT(I,J,L),                           &  
                                   RAD_CFT(I,J,L),                                           &
                                   RAD_RLT(I,J,L), RAD_RIT(I,J,L),                           &
                                   TEMPOR(I,J)                                               )
                endDO
             endDO
          endDO

          CWCT(:,:,:,1) = RAD_QIT !QI
          CWCT(:,:,:,2) = RAD_QLT !QL
          CWCT(:,:,:,3) = 0.0     !QR - this is not generally available in the linear model
          CWCT(:,:,:,4) = 0.0     !QI - this is not generally available in the linear model

          REFFT(:,:,:,1) = RAD_RIT * 1.0e6 !RI
          REFFT(:,:,:,2) = RAD_RLT * 1.0e6 !RL
          REFFT(:,:,:,3) = 100.e-6 * 1.0e6 !RR
          REFFT(:,:,:,4) = 140.e-6 * 1.0e6 !RS

          !Determine the model level seperating high and middle clouds
          LCLDMH = 1
          do L = 1, LM
             if ( PREF(L) >= 40000. ) then
                LCLDMH = L
                exit
             end if
          end do

          !Determine the model level seperating low and middle clouds
          LCLDLM = LM
          do L = 1, LM
             if ( PREF(L) >= 70000.  ) then
                LCLDLM = L
                exit
             end if
          end do

          !Set surface quantities
          EGT = 0.0
          do L = 1, 10
             EGT(:,:,1,L) = EMIS(:,:)
          end do
          FST         = 1.0
          TGT         = 0.0
          TGT(:,:,1)  = TS
          TVT         = 0.0
          TVT(:,:,1)  = TS
          EVT         = 0.0
          RVT         = 0.0

          !Get IRRAD aerosols
          if (NA .eq. 0) then

             TAUA_IRT = 0.0
             SSAA_IRT = 0.0
             ASYA_IRT = 0.0

             TAUA_SOT = 0.0
             SSAA_SOT = 0.0
             ASYA_SOT = 0.0

          elseif (NA .gt. 0) then

             !Place aerosol fields into rank 4 array that can be looped over
             AEROT(:,:,:,1) = DU001T
             AEROT(:,:,:,2) = DU002T
             AEROT(:,:,:,3) = DU003T
             AEROT(:,:,:,4) = DU004T
             AEROT(:,:,:,5) = DU005T

             !Names of aerosols as named in NL model
             AEROSOLS(1) = 'du001'
             AEROSOLS(2) = 'du002'
             AEROSOLS(3) = 'du003'
             AEROSOLS(4) = 'du004'
             AEROSOLS(5) = 'du005'

             !Optical property calculation needs pressure as (1/g)dp/dz
             DO L = 1, LM
                DO J = 1, JM
                   DO I = 1, IM
                      X = ((PLE(I,J,L) - PLE(I,J,L-1))*0.01)*(100./MAPL8_GRAV)
                      DO AI = 1, NA
                         SAEROT(I,J,L,AI) = X*AEROT(I,J,L,AI)
                      END DO
                   END DO
                END DO
             END DO

             !Get aerosol optical properties using Mie Table lookup:
             ! TAUA - Aerosol optical thickness
             ! SSAA - Single scattering albedo
             ! ASYA - Symmetry factor

             ! Note that the aerosol optical properties are dependent on relative humidity, except for dust.
             ! This means that while only considering dust any relative humidity can be used and perturbations
             ! of RH are not required. If adding other aerosols this will need to be changed. For dust SSAAP 
             ! and ASYAP are zero.

             !Note that the TL/AD code calls a modified Get_AeroOptPropAllAero that just returns the table values
             !rather than the optical properties themselves. This saves making two calls, one for the trajectory
             !and one for the perturbation, and makes linearization more straightforward. A more involved linearsation
             !of Get_AeroOptPropAllAero is required if extending to perturbations of RH for non-dust aerosols.

             ! Tables (variables ending in _C) are in single precision

             !Doing dust only so fix RH at some constant between 0 and 1.
             RH = 0.5

             !Initialize aerosol tables. [LooseEnd - this should really be in an initilize to save time]
             call Create_AeroOptProp(CF,STATUS)
             VERIFY_(STATUS)

             !Compute aerosol optical properties for IR
             OFFSET = NBCHOU_SO
             !CALL Get_AeroOptPropAllAero_tlad(NA,AEROSOLS,NBCHOU_IR,OFFSET,RH,TAUA_IR_C,SSAA_IR_C,ASYA_IR_C,STATUS)

             TAUA_IRT = 0.0
             SSAA_IRT = 0.0
             ASYA_IRT = 0.0
             do J = 1,NBCHOU_IR
                do L = 1,NA
                   TAUA_IRT(:,:,:,J) = TAUA_IRT(:,:,:,J) + SAEROT(:,:,:,L)*TAUA_IR_C(:,:,:,J,L)
                   SSAA_IRT(:,:,:,J) = SSAA_IRT(:,:,:,J) + SAEROT(:,:,:,L)*TAUA_IR_C(:,:,:,J,L)*SSAA_IR_C(:,:,:,J,L)
                   ASYA_IRT(:,:,:,J) = ASYA_IRT(:,:,:,J) + SAEROT(:,:,:,L)*TAUA_IR_C(:,:,:,J,L)*SSAA_IR_C(:,:,:,J,L)*ASYA_IR_C(:,:,:,J,L)
                enddo
             enddo

             !Compute aerosol optical properties for SO
             OFFSET = 0
             !CALL Get_AeroOptPropAllAero_tlad(NA,AEROSOLS,NBCHOU_SO,OFFSET,RH,TAUA_SO_C,SSAA_SO_C,ASYA_SO_C,STATUS)

             TAUA_SOT = 0.0
             SSAA_SOT = 0.0
             ASYA_SOT = 0.0
             do J = 1,NBCHOU_SO
                do L = 1,NA
                   TAUA_SOT(:,:,:,J) = TAUA_SOT(:,:,:,J) + SAEROT(:,:,:,L)*TAUA_SO_C(:,:,:,J,L)
                   SSAA_SOT(:,:,:,J) = SSAA_SOT(:,:,:,J) + SAEROT(:,:,:,L)*TAUA_SO_C(:,:,:,J,L)*SSAA_SO_C(:,:,:,J,L)
                   ASYA_SOT(:,:,:,J) = ASYA_SOT(:,:,:,J) + SAEROT(:,:,:,L)*TAUA_SO_C(:,:,:,J,L)*SSAA_SO_C(:,:,:,J,L)*ASYA_SO_C(:,:,:,J,L)
                enddo
             enddo

          endif

          !if (SolCycFileName /= '/dev/null') THEN
          !   call MAPL_SunGetSolarConstant(CLOCK,trim(SolCycFileName),SC,HK=HK_4,rc=STATUS)
          !   VERIFY_(STATUS)
 
          !   HK = dble(HK_4)
          !   HK_UV_TEMP = HK(:5)
       
          !   do L=1,3
          !      HK_IR_TEMP(L,:)=HK_IR_OLD(L,:)*(HK(5+L)/sum(HK_IR_OLD(L,:)))
          !   end do
          !else if (SC<0.0) then
          !   call MAPL_SunGetSolarConstant(CURRENTTIME,SC,HK_4)
          !   HK = dble(HK_4)
          !   HK_UV_TEMP = HK(:5)
          !   do L=1,3
          !      HK_IR_TEMP(L,:)=HK_IR_OLD(L,:)*(HK(5+L)/sum(HK_IR_OLD(L,:)))
          !   end do
          !else
          !   HK_UV_TEMP = HK_UV_OLD
          !   HK_IR_TEMP = HK_IR_OLD
          !end if

          COSZTMP = max(.0001,COSZ)

          !Initialize the fluxes and outputs
          IR_UFLXT = 0.0
          IR_DFLXT = 0.0
          IR_DFDTsT = 0.0

          SO_NFLXT = 0.0

          !Copy rad constants to double precision
          aib_ir_8 = dble(aib_ir)
          awb_ir_8 = dble(awb_ir)
          aiw_ir_8 = dble(aiw_ir)
          aww_ir_8 = dble(aww_ir)
          aig_ir_8 = dble(aig_ir)
          awg_ir_8 = dble(awg_ir)
          xkw_8 = dble(xkw)
          xke_8 = dble(xke)
          aw_8 = dble(aw)
          bw_8 = dble(bw)
          pm_8 = dble(pm)
          fkw_8 = dble(fkw)
          gkw_8 = dble(gkw)
          cb_8 = dble(cb)
          dcb_8 = dble(dcb)
          w11_8 = dble(w11)
          w12_8 = dble(w12)
          w13_8 = dble(w13)
          p11_8 = dble(p11)
          p12_8 = dble(p12)
          p13_8 = dble(p13)
          dwe_8 = dble(dwe)
          dpe_8 = dble(dpe)
          c1_8 = dble(c1)
          c2_8 = dble(c2)
          c3_8 = dble(c3)
          oo1_8 = dble(oo1)
          oo2_8 = dble(oo2)
          oo3_8 = dble(oo3)
          h11_8 = dble(h11)
          h12_8 = dble(h12)
          h13_8 = dble(h13)
          h21_8 = dble(h21)
          h22_8 = dble(h22)
          h23_8 = dble(h23)
          h81_8 = dble(h81)
          h82_8 = dble(h82)
          h83_8 = dble(h83)

          wk_uv_8 = dble(wk_uv)
          zk_uv_8 = dble(zk_uv)
          ry_uv_8 = dble(ry_uv)
          xk_ir_8 = dble(xk_ir)
          ry_ir_8 = dble(ry_ir)
          cah_8 = dble(cah)
          coa_8 = dble(coa)
          aig_uv_8 = dble(aig_uv)
          awg_uv_8 = dble(awg_uv)
          arg_uv_8 = dble(arg_uv)
          aib_uv_8 = dble(aib_uv)
          awb_uv_8 = dble(awb_uv)
          arb_uv_8 = dble(arb_uv)
          aib_nir_8 = dble(aib_nir)
          awb_nir_8 = dble(awb_nir)
          arb_nir_8 = dble(arb_nir)
          aia_nir_8 = dble(aia_nir)
          awa_nir_8 = dble(awa_nir)
          ara_nir_8 = dble(ara_nir)
          aig_nir_8 = dble(aig_nir)
          awg_nir_8 = dble(awg_nir)
          arg_nir_8 = dble(arg_nir)
          caib_8 = dble(caib)
          caif_8 = dble(caif)

! Perturbation calulations down when alarm in ringing
! ---------------------------------------------------

          if (PHASE==TLMphase) then 

! Tangent linear of preparation of inputs for radiation
! -----------------------------------------------------

             !Compute 2m temperature
             T2MP = TP(:,:,LM)*(0.5*(1.0 + PLE(:,:,LM-1)/PLE(:,:,LM)))**(-MAPL8_KAPPA)       

             !Ozone ppmv to kgkg
             O3mmP = O3P / 1.e6

             !Splitting to large scale and convective parts
             QILSP = QIP * ILSF
             QICNP = QIP * ICNF
             QLLSP = QLP * LLSF
             QLCNP = QLP * LCNF

             !Initialize RAD/CLOUD variables
             RAD_CFT = 0.0
             RAD_QLT = 0.0
             RAD_QIT = 0.0
             RAD_RLT = 0.0
             RAD_RIT = 0.0
             RAD_CFP = 0.0
             RAD_QLP = 0.0
             RAD_QIP = 0.0
             RAD_RLP = 0.0
             RAD_RIP = 0.0

             DO I = 1,IM
                DO J = 1,JM
                   DO L = 1,LM
                      call RADCOUPLE_D( TT(I,J,L),      TP(I,J,L),         &
                                        PLOf(I,J,L),                       &
                                        CFLST(I,J,L),   CFLSP(I,J,L),      &
                                        CFCNT(I,J,L),   CFCNP(I,J,L),      &
                                        QLLST(I,J,L),   QLLSP(I,J,L),      & 
                                        QILST(I,J,L),   QILSP(I,J,L),      &
                                        QLCNT(I,J,L),   QLCNP(I,J,L),      & 
                                        QICNT(I,J,L),   QICNP(I,J,L),      & 
                                        RAD_QLT(I,J,L), RAD_QLP(I,J,L),    & 
                                        RAD_QIT(I,J,L), RAD_QIP(I,J,L),    &  
                                        RAD_CFT(I,J,L), RAD_CFP(I,J,L),    & 
                                        RAD_RLT(I,J,L), RAD_RLP(I,J,L),    & 
                                        RAD_RIT(I,J,L), RAD_RIP(I,J,L),    & 
                                        TEMPOR(I,J)                        )
                   endDO
                endDO
             endDO

             CWCP(:,:,:,1) = RAD_QIP
             CWCP(:,:,:,2) = RAD_QLP
             CWCP(:,:,:,3) = 0.0
             CWCP(:,:,:,4) = 0.0

             REFFP(:,:,:,1) = RAD_RIP * 1.0e6
             REFFP(:,:,:,2) = RAD_RLP * 1.0e6
             REFFP(:,:,:,3) = 0.0
             REFFP(:,:,:,4) = 0.0 

             !Aerosols mass to optical
             if (NA .eq. 0) then

                TAUA_IRP = 0.0
                SSAA_IRP = 0.0
                ASYA_IRP = 0.0

                TAUA_SOP = 0.0
                SSAA_SOP = 0.0
                ASYA_SOP = 0.0

             elseif (NA .gt. 0) then

                AEROP = 0.0
                AEROP(:,:,:,1) = DU001P
                AEROP(:,:,:,2) = DU002P
                AEROP(:,:,:,3) = DU003P
                AEROP(:,:,:,4) = DU004P
                AEROP(:,:,:,5) = DU005P

                SAEROP = 0.0
                DO L = 1, LM
                   DO J = 1, JM
                      DO I = 1, IM
                         X = ((PLE(I,J,L) - PLE(I,J,L-1))*0.01)*(100./MAPL8_GRAV)
                         DO AI = 1, NA
                            SAEROP(I,J,L,AI) = X*AEROP(I,J,L,AI)
                         END DO
                      END DO
                   END DO
                END DO

                !Compute perturbation aerosol optical properties for IR
                TAUA_IRP = 0.0
                SSAA_IRP = 0.0
                ASYA_IRP = 0.0
                do J = 1,NBCHOU_IR
                   do L = 1,NA
                      TAUA_IRP(:,:,:,J) = TAUA_IRP(:,:,:,J) &
                                        + SAEROP(:,:,:,L)*TAUA_IR_C(:,:,:,J,L)
                      SSAA_IRP(:,:,:,J) = SSAA_IRP(:,:,:,J) &
                                        + SAEROP(:,:,:,L)*TAUA_IR_C(:,:,:,J,L)*SSAA_IR_C(:,:,:,J,L)
                      ASYA_IRP(:,:,:,J) = ASYA_IRP(:,:,:,J) &
                                        + SAEROP(:,:,:,L)*TAUA_IR_C(:,:,:,J,L)*SSAA_IR_C(:,:,:,J,L)*ASYA_IR_C(:,:,:,J,L)
                   enddo
                enddo

                !Compute perturbation aerosol optical properties for SO
                TAUA_SOP = 0.0
                SSAA_SOP = 0.0
                ASYA_SOP = 0.0
                do J = 1,NBCHOU_SO
                   do L = 1,NA
                      TAUA_SOP(:,:,:,J) = TAUA_SOP(:,:,:,J) &
                                        + SAEROP(:,:,:,L)*TAUA_SO_C(:,:,:,J,L)
                      SSAA_SOP(:,:,:,J) = SSAA_SOP(:,:,:,J) &
                                        + SAEROP(:,:,:,L)*TAUA_SO_C(:,:,:,J,L)*SSAA_SO_C(:,:,:,J,L)
                      ASYA_SOP(:,:,:,J) = ASYA_SOP(:,:,:,J) &
                                        + SAEROP(:,:,:,L)*TAUA_SO_C(:,:,:,J,L)*SSAA_SO_C(:,:,:,J,L)*ASYA_SO_C(:,:,:,J,L)
                   enddo
                enddo

             END IF

! Tangnet linear short wave (SO) radiation scheme
! -----------------------------------------------

             !Initialize the perturbation fluxes
             SO_NFLXP = 0.0

             DO I = 1,IM
                DO J = 1,JM

                   !Call scheme, looped to minimise checkpoint memory demands
                   CALL SORAD_D( 1,                                    &
                                 LM,                                   &
                                 NBCHOU_SO,                            &
                                 COSZTMP(I,J),                         &
                                 PLESO(I,J,:),                         &
                                 TT(I,J,:), TP(I,J,:),                 &
                                 QVT(I,J,:), QVP(I,J,:),               &
                                 O3mmT(I,J,:), O3mmP(I,J,:),           &
                                 CO2,                                  &
                                 CWCT(I,J,:,:), CWCP(I,J,:,:),         &
                                 RAD_CFT(I,J,:), RAD_CFP(I,J,:),       &
                                 LCLDMH, LCLDLM,                       &
                                 REFFT(I,J,:,:), REFFP(I,J,:,:),       &
                                 HK_UV_TEMP, HK_IR_TEMP,               &
                                 TAUA_SOT(I,J,:,:), TAUA_SOP(I,J,:,:), &
                                 SSAA_SOT(I,J,:,:), SSAA_SOP(I,J,:,:), &
                                 ASYA_SOT(I,J,:,:), ASYA_SOP(I,J,:,:), &
                                 RGBUV(I,J),                           &
                                 RGFUV(I,J),                           &
                                 RGBIR(I,J),                           & 
                                 RGFIR(I,J),                           & 
                                 SO_NFLXT(I,J,:), SO_NFLXP(I,J,:),     &
                                 MAPL8_GRAV,                           &
                                 wk_uv_8, zk_uv_8, ry_uv_8,            &
                                 xk_ir_8, ry_ir_8,                     & 
                                 cah_8, coa_8,                         &
                                 aig_uv_8, awg_uv_8, arg_uv_8,         &
                                 aib_uv_8, awb_uv_8, arb_uv_8,         &
                                 aib_nir_8, awb_nir_8, arb_nir_8,      &
                                 aia_nir_8, awa_nir_8, ara_nir_8,      &
                                 aig_nir_8, awg_nir_8, arg_nir_8,      &
                                 caib_8, caif_8                        )

                endDO
             endDO

             !Put perturbation fluxes in internal state
             SO_NFLXP_TL_INT = real(SO_NFLXP,4)

! Tangnet linear long wave (IR) radiation scheme
! ----------------------------------------------

             !Initialize the perturbation fluxes
             IR_UFLXP = 0.0
             IR_DFLXP = 0.0
             IR_DFDTsP = 0.0

             !Call scheme
             DO I = 1,IM
                DO J = 1,JM

                   CALL IRRAD_D( LM,                                     &
                                 PLE(I,J,:),                             &
                                 TT(I,J,:), TP(I,J,:),                   &
                                 QVT(I,J,:), QVP(I,J,:),                 &
                                 O3mmT(I,J,:), O3mmP(I,J,:),             &
                                 T2MT(I,J), T2MP(I,J),                   &
                                 CO2,                                    &
                                 TRACE,                                  &
                                 N2OT(I,J,:),                            &
                                 CH4T(I,J,:),                            &
                                 CFC11T(I,J,:),                          &
                                 CFC12T(I,J,:),                          &
                                 CFC22T(I,J,:),                          &
                                 CWCT(I,J,:, :), CWCP(I,J,:,:),          &
                                 RAD_CFT(I,J,:), RAD_CFP(I,J,:),         &
                                 LCLDMH, LCLDLM,                         &
                                 REFFT(I,J,:,:), REFFP(I,J,:, :),        &
                                 NS,                                     &
                                 FST(I,J,:),                             &
                                 TGT(I,J,:),                             &
                                 EGT(I,J,:,:),                           &
                                 TVT(I,J,:),                             &
                                 EVT(I,J,:,:),                           &
                                 RVT(I,J,:,:),                           &
                                 NA,                                     &
                                 NBCHOU_IR,                              &
                                 TAUA_IRT(I,J,:,:), TAUA_IRP(I,J,:,:),   &
                                 SSAA_IRT(I,J,:,:), SSAA_IRP(I,J,:,:),   &
                                 ASYA_IRT(I,J,:,:), ASYA_IRP(I,J,:,:),   &
                                 IR_UFLXT(I,J,:), IR_UFLXP(I,J,:),       &
                                 IR_DFLXT(I,J,:), IR_DFLXP(I,J,:),       &
                                 IR_DFDTsT(I,J,:), IR_DFDTsP(I,J,:),     &
                                 aib_ir_8, awb_ir_8, aiw_ir_8,           &
                                 aww_ir_8, aig_ir_8, awg_ir_8,           &
                                 xkw_8, xke_8,                           &
                                 mw, aw_8, bw_8, pm_8,                   &
                                 fkw_8, gkw_8, cb_8, dcb_8,              &
                                 w11_8, w12_8, w13_8, p11_8,             &
                                 p12_8, p13_8, dwe_8, dpe_8,             &
                                 c1_8, c2_8, c3_8,                       &
                                 oo1_8, oo2_8, oo3_8,                    &
                                 h11_8, h12_8, h13_8,                    &
                                 h21_8, h22_8, h23_8,                    &
                                 h81_8, h82_8, h83_8                     )

                endDO
             endDO

             !Compute net flux for IR
             IR_NFLXP = IR_UFLXP + IR_DFLXP

             !Place fluxes in internal states
             IR_NFLXP_TL_INT = real(IR_NFLXP,4)
             IR_DFDTsP_TL_INT = real(IR_DFDTsP,4)

       elseif (PHASE==ADJphase) then

! Adjoint long wave (IR) radiation scheme
! ---------------------------------------


             !Initilize perturbation quantities output by adjoint IR radiation
             T2MP = 0.0
             O3mmP = 0.0
             CWCP = 0.0
             RAD_CFP = 0.0
             REFFP = 0.0
             TAUA_IRP = 0.0
             SSAA_IRP = 0.0
             ASYA_IRP = 0.0

             IR_DFLXP = 0.0
             IR_UFLXP = 0.0
             IR_UFLXP = IR_NFLXP
             IR_DFLXP = IR_NFLXP

             DO I = 1,IM
                DO J = 1,JM

                   CALL IRRAD_B( LM,                                     &
                                 PLE(I,J,:),                             &
                                 TT(I,J,:), TP(I,J,:),                   &
                                 QVT(I,J,:), QVP(I,J,:),                 &
                                 O3mmT(I,J,:), O3mmP(I,J,:),             &
                                 T2MT(I,J), T2MP(I,J),                   &
                                 CO2,                                    &
                                 TRACE,                                  &
                                 N2OT(I,J,:),                            &
                                 CH4T(I,J,:),                            &
                                 CFC11T(I,J,:),                          &
                                 CFC12T(I,J,:),                          &
                                 CFC22T(I,J,:),                          &
                                 CWCT(I,J,:, :), CWCP(I,J,:,:),          &
                                 RAD_CFT(I,J,:), RAD_CFP(I,J,:),         &
                                 LCLDMH, LCLDLM,                         &
                                 REFFT(I,J,:,:), REFFP(I,J,:, :),        &
                                 NS,                                     &
                                 FST(I,J,:),                             &
                                 TGT(I,J,:),                             &
                                 EGT(I,J,:,:),                           &
                                 TVT(I,J,:),                             &
                                 EVT(I,J,:,:),                           &
                                 RVT(I,J,:,:),                           &
                                 NA,                                     &
                                 NBCHOU_IR,                              &
                                 TAUA_IRT(I,J,:,:), TAUA_IRP(I,J,:,:),   &
                                 SSAA_IRT(I,J,:,:), SSAA_IRP(I,J,:,:),   &
                                 ASYA_IRT(I,J,:,:), ASYA_IRP(I,J,:,:),   &
                                 IR_UFLXT(I,J,:), IR_UFLXP(I,J,:),       &
                                 IR_DFLXT(I,J,:), IR_DFLXP(I,J,:),       &
                                 IR_DFDTsT(I,J,:), IR_DFDTsP(I,J,:),     &
                                 aib_ir_8, awb_ir_8, aiw_ir_8,           &
                                 aww_ir_8, aig_ir_8, awg_ir_8,           &
                                 xkw_8, xke_8,                           &
                                 mw, aw_8, bw_8, pm_8,                   &
                                 fkw_8, gkw_8, cb_8, dcb_8,              &
                                 w11_8, w12_8, w13_8, p11_8,             &
                                 p12_8, p13_8, dwe_8, dpe_8,             &
                                 c1_8, c2_8, c3_8,                       &
                                 oo1_8, oo2_8, oo3_8,                    &
                                 h11_8, h12_8, h13_8,                    &
                                 h21_8, h22_8, h23_8,                    &
                                 h81_8, h82_8, h83_8                     )

                endDO
             endDO
 
! Adjoint short wave (SO) radiation scheme
! ----------------------------------------

             TAUA_SOP = 0.0
             SSAA_SOP = 0.0
             ASYA_SOP = 0.0

             DO I = 1,IM
                DO J = 1,JM

                   !Call scheme, looped to minimise checkpoint memory demands
                   CALL SORAD_B( 1,                                    &
                                 LM,                                   &
                                 NBCHOU_SO,                            &
                                 COSZTMP(I,J),                         &
                                 PLESO(I,J,:),                         &
                                 TT(I,J,:), TP(I,J,:),                 &
                                 QVT(I,J,:), QVP(I,J,:),               &
                                 O3mmT(I,J,:), O3mmP(I,J,:),           &
                                 CO2,                                  &
                                 CWCT(I,J,:,:), CWCP(I,J,:,:),         &
                                 RAD_CFT(I,J,:), RAD_CFP(I,J,:),       &
                                 LCLDMH, LCLDLM,                       &
                                 REFFT(I,J,:,:), REFFP(I,J,:,:),       &
                                 HK_UV_TEMP, HK_IR_TEMP,               &
                                 TAUA_SOT(I,J,:,:), TAUA_SOP(I,J,:,:), &
                                 SSAA_SOT(I,J,:,:), SSAA_SOP(I,J,:,:), &
                                 ASYA_SOT(I,J,:,:), ASYA_SOP(I,J,:,:), &
                                 RGBUV(I,J),                           &
                                 RGFUV(I,J),                           &
                                 RGBIR(I,J),                           & 
                                 RGFIR(I,J),                           & 
                                 SO_NFLXT(I,J,:), SO_NFLXP(I,J,:),     &
                                 MAPL8_GRAV,                           &
                                 wk_uv_8, zk_uv_8, ry_uv_8,            &
                                 xk_ir_8, ry_ir_8,                     & 
                                 cah_8, coa_8,                         &
                                 aig_uv_8, awg_uv_8, arg_uv_8,         &
                                 aib_uv_8, awb_uv_8, arb_uv_8,         &
                                 aib_nir_8, awb_nir_8, arb_nir_8,      &
                                 aia_nir_8, awa_nir_8, ara_nir_8,      &
                                 aig_nir_8, awg_nir_8, arg_nir_8,      &
                                 caib_8, caif_8                        )

                endDO
             endDO

! Adjoint of preparing inputs for radiation schemes
! -------------------------------------------------

             if (NA .ne. 0) then

                SAEROP = 0.0

                !Solar aerosol optical properties
                DO J = NBCHOU_SO,1,-1
                   DO L = NA,1,-1
                      SAEROP(:,:,:, L) = SAEROP(:,:,:,L)                                               &
                                       + TAUA_SO_C(:,:,:,J,L)*TAUA_SOP(:,:,:,J)                        &
                                       + TAUA_SO_C(:,:,:,J,L)*SSAA_SOP(:,:,:,J)*SSAA_SO_C(:,:,:,J,L)   &
                                       + TAUA_SO_C(:,:,:,J,L)*ASYA_SOP(:,:,:,J)*SSAA_SO_C(:,:,:,J,L)*ASYA_SO_C(:,:,:,J,L)
                   END DO
                END DO

                !IR aerosol optical properties
                DO J = NBCHOU_IR,1,-1
                   DO L = NA,1,-1
                      SAEROP(:,:,:, L) = SAEROP(:,:,:,L)                                               &
                                       + TAUA_IR_C(:,:,:,J,L)*TAUA_IRP(:,:,:,J)                        &
                                       + TAUA_IR_C(:,:,:,J,L)*SSAA_IRP(:,:,:,J)*SSAA_IR_C(:,:,:,J,L)   &
                                       + TAUA_IR_C(:,:,:,J,L)*ASYA_IRP(:,:,:,J)*SSAA_IR_C(:,:,:,J,L)*ASYA_IR_C(:,:,:,J,L)
                   END DO
                END DO

                AEROP = 0.0
                DO L = LM,1,-1
                   DO J = JM,1,-1
                      DO I = IM,1,-1

                         X = ((PLE(I,J,L) - PLE(I,J,L-1))*0.01)*(100./MAPL8_GRAV)

                         DO AI = NA,1,-1

                            AEROP(I,J,L,AI)  = AEROP(I,J,L,AI) + X * SAEROP(I,J,L,AI)
                            SAEROP(I,J,L,AI) = 0.0

                      END DO
                    END DO
                  END DO
                END DO

                DU001P = DU001P + AEROP(:,:,:,1)
                DU002P = DU002P + AEROP(:,:,:,2)
                DU003P = DU003P + AEROP(:,:,:,3)
                DU004P = DU004P + AEROP(:,:,:,4)
                DU005P = DU005P + AEROP(:,:,:,5)
                AEROP = 0.0

             endif

             RAD_RIP = 0.0
             RAD_RIP = 1.0e6*REFFP(:,:,:, 1)
             REFFP(:,:,:,1) = 0.0
             RAD_RLP = 0.0
             RAD_RLP = 1.0e6*REFFP(:,:,:, 2)
             REFFP(:,:,:,2) = 0.0

             REFFP(:,:,:,3) = 0.0
             REFFP(:,:,:,4) = 0.0

             RAD_QIP = 0.0
             RAD_QIP = CWCP(:,:,:, 1)
             CWCP(:,:,:,1) = 0.0
             RAD_QLP = 0.0
             RAD_QLP = CWCP(:,:,:, 2)
             CWCP(:,:,:,2) = 0.0

             CWCP(:,:,:,3) = 0.0
             CWCP(:,:,:,4) = 0.0

             QLCNP = 0.0
             QLLSP = 0.0
             QICNP = 0.0
             QILSP = 0.0

             RAD_CFT = 0.0
             RAD_QLT = 0.0
             RAD_QIT = 0.0
             RAD_RLT = 0.0
             RAD_RIT = 0.0

             DO I = IM,1,-1
               DO J = JM,1,-1
                 DO L = LM,1,-1
                      call RADCOUPLE_b( TT(I,J,L),      TP(I,J,L),         &
                                        PLOf(I,J,L),                       &
                                        CFLST(I,J,L),   CFLSP(I,J,L),      &
                                        CFCNT(I,J,L),   CFCNP(I,J,L),      &
                                        QLLST(I,J,L),   QLLSP(I,J,L),      & 
                                        QILST(I,J,L),   QILSP(I,J,L),      &
                                        QLCNT(I,J,L),   QLCNP(I,J,L),      & 
                                        QICNT(I,J,L),   QICNP(I,J,L),      & 
             	                        RAD_QLT(I,J,L), RAD_QLP(I,J,L),    & 
                                        RAD_QIT(I,J,L), RAD_QIP(I,J,L),    &  
                                        RAD_CFT(I,J,L), RAD_CFP(I,J,L),    & 
                                        RAD_RLT(I,J,L), RAD_RLP(I,J,L),    & 
                                        RAD_RIT(I,J,L), RAD_RIP(I,J,L),    & 
                                        TEMPOR(I,J)                        )
                      RAD_QLP(I,J,L) = 0.0
                      RAD_QIP(I,J,L) = 0.0
                      RAD_CFP(I,J,L) = 0.0
                      RAD_RLP(I,J,L) = 0.0
                      RAD_RIP(I,J,L) = 0.0
                   END DO
                END DO
             END DO

             QIP = QIP + QILSP * ILSF + QICNP * ICNF
             QLP = QLP + QLLSP * LLSF + QLCNP * LCNF

             O3P = O3P + O3mmP / 1.e6
             TP(:,:,LM) = TP(:,:,LM) + T2MP*(0.5*(1.0 + PLE(:,:,LM-1)/PLE(:,:,LM)))**(-MAPL8_KAPPA)

             !Now that we are done with internal fluxes for this phase, reset them to zero.
             IR_NFLXP_AD_INT = 0.0
             SO_NFLXP_AD_INT = 0.0
             IR_DFDTsP_AD_INT = 0.0

             !Go back to single precision
             QVP4    = real(QVP,4)
             O3P4    = real(O3P,4)
             QIP4    = real(QIP,4)
             QLP4    = real(QLP,4)
             CFLSP4  = real(CFLSP,4)
             CFCNP4  = real(CFCNP,4)
             if (NA > 0) then !Doing aerosols
                DU001P4 = real(DU001P,4)
                DU002P4 = real(DU002P,4)
                DU003P4 = real(DU003P,4)
                DU004P4 = real(DU004P,4)
                DU005P4 = real(DU005P,4)
             endif

          endif 

! Deallocate memory used in the flux recalculation
! ------------------------------------------------

          deallocate(UT,QVT,O3T,CFLST,CFCNT,QVP,O3P,CFLSP,CFCNP)
          deallocate(QLST,QCNT,QIP,QLP,TS)
          deallocate(EMIS,COSZ,RGBUV,RGFUV,RGBIR,RGFIR)
          deallocate(T2MT,T2MP,PLESO,COSZTMP)
          deallocate(O3mmT,O3mmP,N2OT,CH4T,CFC11T,CFC12T,CFC22T)
          deallocate(QILST,QLLST,QICNT,QLCNT,QILSP,QLLSP,QICNP,QLCNP)
          deallocate(QIT,QLT,fQi,ILSF,ICNF,LLSF,LCNF,CWCT,CWCP)
          deallocate(REFFT,REFFP,TEMPOR)
          deallocate(RAD_QLT,RAD_QLP,RAD_QIT,RAD_QIP,RAD_CFT)
          deallocate(RAD_CFP,RAD_RLT,RAD_RLP,RAD_RIT,RAD_RIP)
          deallocate(FST,TGT,TVT,EGT,EVT,RVT)
          deallocate(aib_ir_8,awb_ir_8,aiw_ir_8,aww_ir_8,aig_ir_8,awg_ir_8)
          deallocate(xkw_8,xke_8,aw_8,bw_8,pm_8,fkw_8,gkw_8,cb_8,dcb_8)
          deallocate(c1_8,c2_8,c3_8,oo1_8,oo2_8,oo3_8)
          deallocate(h11_8,h12_8,h13_8,h21_8,h22_8,h23_8,h81_8,h82_8,h83_8)
          deallocate(wk_uv_8,zk_uv_8,ry_uv_8,xk_ir_8,ry_ir_8)
          deallocate(cah_8,coa_8,aig_uv_8,awg_uv_8,arg_uv_8)
          deallocate(awb_uv_8,arb_uv_8,awb_nir_8,arb_nir_8,aia_nir_8,awa_nir_8)
          deallocate(ara_nir_8,aig_nir_8,awg_nir_8,arg_nir_8,caib_8,caif_8)
          if (NA > 0) then
             deallocate(DU001T,DU002T,DU003T,DU004T,DU005T)
             deallocate(DU001P,DU002P,DU003P,DU004P,DU005P)
             deallocate(AEROSOLS,AEROT,AEROP,SAEROT,SAEROP,RH)
             deallocate(TAUA_IR_C,SSAA_IR_C,ASYA_IR_C)
             deallocate(TAUA_IRT,TAUA_IRP,SSAA_IRT,SSAA_IRP,ASYA_IRT,ASYA_IRP)
             deallocate(TAUA_SO_C,SSAA_SO_C,ASYA_SO_C)
             deallocate(TAUA_SOT,TAUA_SOP,SSAA_SOT,SSAA_SOP,ASYA_SOT,ASYA_SOP)
          endif

       endif

! Post flux update perturbation caluclulations done at every time step
! --------------------------------------------------------------------

       if (PHASE==TLMphase) then 

          if (RUNALARM == 1) then
             !Use just computed fluxes, saves loss of dot product accuracy
             IR_NFLXP = IR_NFLXP
             SO_NFLXP = SO_NFLXP
             IR_DFDTsP = IR_DFDTsP
          else
             !Grab most recent fluxes from internal state
             IR_NFLXP = dble(IR_NFLXP_TL_INT)
             SO_NFLXP = dble(SO_NFLXP_TL_INT)
             IR_DFDTsP = dble(IR_DFDTsP_TL_INT)
          endif

          IR_FLXP = 0.0
          SO_FLXP = 0.0
          DO L = 0,LM
             IR_FLXP(:,:,L) = IR_NFLXP(:,:,L) + DELT*IR_DFDTsP(:,:,L)
             SO_FLXP(:,:,L) = SLR*SO_NFLXP(:,:,L)
          END DO

          !Weighted temperature tendency
          DTDTP = (MAPL8_GRAV/MAPL8_CP) * (  IR_FLXP(:,:,0:LM-1)-IR_FLXP(:,:,1:LM) &
                                           + SO_FLXP(:,:,0:LM-1)-SO_FLXP(:,:,1:LM) )

          !Unweighted new temperature
          TP = TP + DT * DTDtP / DELTAP

          !FV scaled potential temperature
          PTP = TP/PKT

          !Go back to single precision
          PTP4 = real(PTP,4)

       elseif (PHASE==ADJphase) then

          PTP = PKT * TP

          !Go back to single precision
          PTP4 = real(PTP,4)

       endif

! Deallocate variables used every time step
! -----------------------------------------

       deallocate(PTT,PTP,PS,SLR,DELT)
       deallocate(TT,TP,PLE,PLOf,DELTAP,VERT_AK,VERT_BK,PREF,DTDtP)
       deallocate(IR_UFLXT,IR_UFLXP,IR_DFLXT,IR_DFLXP,IR_NFLXP,IR_DFDTsT,IR_DFDTsP,IR_FLXP)
       deallocate(SO_NFLXT,SO_NFLXP,SO_FLXP)
       deallocate(PKT,PET,PKE,BX,CX)

! Turn off TLM/ADM timers
! -----------------------

       if (PHASE==ADJphase) then         
          call timing_off('RADIA_ADM')             
       elseif (PHASE==TLMphase) then         
          call timing_off('RADIA_TLM')            
       endif

    ENDIF !DO_RAD_PHYS

! Turn off MAPL timers
! --------------------

    call MAPL_TimerOff(MAPL,PHASE_NAME) 
    call MAPL_TimerOff(MAPL,"-RADIATION")
    call MAPL_TimerOff(MAPL,"TOTAL")

! Until next time
! ---------------

    RETURN_(ESMF_SUCCESS)

 contains

 REAL*4 FUNCTION GetCO2Pert(Year, DayOfYear)

! !DESCRIPTION:
!
!  Given the year and day-of-year, this function returns the RCP45 CO2 concentration in 
!  mole fraction (volume mixing ratio).  If Year is less than 1765, the value for 1765
!  is returned.  If Year is greater than 2150, the value for 2150 is returned.  In the
!  original dataset, the value for 2150 is used for all years through 2500.  We choose
!  to truncate the list at 2151 for this application.
!
!  DayOfYear is expected to have a value of 1.00 at 0:00 UTC Jan 1.
!
!  In-line documentation from the source dataset is reproduced below:
! 
! RCP45 Midyear Atmospheric CO2 Concentrations (ppmv)
!
! NEEDS TO BE KEPT UP TO DATE WITH NONLINEAR MODEL
!
! ---------------------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Year
  INTEGER, INTENT(IN) :: DayOfYear

  REAL*4 :: f,i,n
  INTEGER :: previous,current,next
  
  INTEGER, PARAMETER :: firstYear = 1764
  INTEGER, PARAMETER :: finalYear = 2151
  INTEGER, PARAMETER :: tableLength = finalYear-firstYear+1

  REAL*4, SAVE :: CO2ppmv(tableLength) = (/                                278.052, &
   278.052,  278.106,  278.220,  278.343,  278.471,  278.600,  278.733,  278.869, &
   279.009,  279.153,  279.302,  279.457,  279.618,  279.782,  279.943,  280.097, &
   280.243,  280.382,  280.518,  280.657,  280.803,  280.957,  281.118,  281.282, &
   281.443,  281.598,  281.747,  281.891,  282.031,  282.167,  282.299,  282.427, &
   282.551,  282.671,  282.787,  282.899,  283.007,  283.111,  283.211,  283.307, &
   283.400,  283.490,  283.578,  283.661,  283.735,  283.797,  283.847,  283.889, &
   283.926,  283.963,  284.001,  284.043,  284.086,  284.129,  284.167,  284.198, &
   284.223,  284.244,  284.263,  284.281,  284.300,  284.320,  284.340,  284.360, &
   284.380,  284.400,  284.385,  284.280,  284.125,  283.975,  283.825,  283.675, &
   283.525,  283.425,  283.400,  283.400,  283.425,  283.500,  283.600,  283.725, &
   283.900,  284.075,  284.225,  284.400,  284.575,  284.725,  284.875,  285.000, &
   285.125,  285.275,  285.425,  285.575,  285.725,  285.900,  286.075,  286.225, &
   286.375,  286.500,  286.625,  286.775,  286.900,  287.000,  287.100,  287.225, &
   287.375,  287.525,  287.700,  287.900,  288.125,  288.400,  288.700,  289.025, &
   289.400,  289.800,  290.225,  290.700,  291.200,  291.675,  292.125,  292.575, &
   292.975,  293.300,  293.575,  293.800,  294.000,  294.175,  294.325,  294.475, &
   294.600,  294.700,  294.800,  294.900,  295.025,  295.225,  295.500,  295.800, &
   296.125,  296.475,  296.825,  297.200,  297.625,  298.075,  298.500,  298.900, &
   299.300,  299.700,  300.075,  300.425,  300.775,  301.100,  301.400,  301.725, &
   302.075,  302.400,  302.700,  303.025,  303.400,  303.775,  304.125,  304.525, &
   304.975,  305.400,  305.825,  306.300,  306.775,  307.225,  307.700,  308.175, &
   308.600,  309.000,  309.400,  309.750,  310.000,  310.175,  310.300,  310.375, &
   310.375,  310.300,  310.200,  310.125,  310.100,  310.125,  310.200,  310.325, &
   310.500,  310.750,  311.100,  311.500,  311.925,  312.425,  313.000,  313.600, &
   314.225,  314.848,  315.500,  316.272,  317.075,  317.795,  318.397,  318.925, &
   319.647,  320.647,  321.605,  322.635,  323.902,  324.985,  325.855,  327.140, &
   328.677,  329.742,  330.585,  331.747,  333.272,  334.848,  336.525,  338.360, &
   339.728,  340.793,  342.198,  343.783,  345.283,  346.797,  348.645,  350.737, &
   352.487,  353.855,  355.017,  355.885,  356.777,  358.128,  359.837,  361.462, &
   363.155,  365.323,  367.348,  368.865,  370.467,  372.522,  374.760,  376.812, &
   378.812,  380.828,  382.777,  384.800,  386.952,  389.128,  391.274,  393.421, &
   395.583,  397.764,  399.966,  402.184,  404.411,  406.643,  408.882,  411.129, &
   413.378,  415.639,  417.936,  420.274,  422.656,  425.080,  427.538,  430.021, &
   432.523,  435.046,  437.589,  440.131,  442.664,  445.207,  447.770,  450.355, &
   452.963,  455.586,  458.215,  460.845,  463.475,  466.093,  468.678,  471.234, &
   473.780,  476.328,  478.881,  481.438,  483.993,  486.535,  489.060,  491.536, &
   493.932,  496.244,  498.474,  500.645,  502.768,  504.847,  506.884,  508.871, &
   510.799,  512.647,  514.401,  516.065,  517.629,  519.096,  520.488,  521.818, &
   523.089,  524.302,  525.451,  526.509,  527.457,  528.296,  529.027,  529.643, &
   530.144,  530.553,  530.883,  531.138,  531.319,  531.490,  531.702,  531.942, &
   532.205,  532.487,  532.776,  533.070,  533.388,  533.741,  534.131,  534.558, &
   535.011,  535.480,  535.955,  536.435,  536.920,  537.399,  537.871,  538.358, &
   538.872,  539.388,  539.884,  540.352,  540.782,  541.168,  541.510,  541.808, &
   542.053,  542.246,  542.408,  542.559,  542.712,  542.866,  543.013,  543.139, &
   543.239,  543.311,  543.355,  543.360,  543.327,  543.288,  543.266,  543.264, &
   543.282,  543.310,  543.337,  543.354,  543.361,  543.356,  543.330,  543.283, &
   543.238,  543.208,  543.197,  543.205,  543.223,  543.239,  543.246,  543.242, &
   543.227,  543.190,  543.131,  543.072,  543.025,  542.993,  542.979,  542.973, &
   542.963,  542.955,  542.955 /)

! Establish location in table for current, previous and next year.
! ----------------------------------------------------------------
  current  = Year-firstYear+1
  current  = MAX(current,1)
  current  = MIN(current,tableLength)

  previous = current-1
  previous = MAX(previous,1)
  previous = MIN(previous,tableLength)
  IF(Year > finalYear) previous = tableLength

  next     = current+1
  next     = MAX(next,1)
  next     = MIN(next,tableLength)
  IF(Year < firstYear) next = 1

! Divide the year into halves.
! ----------------------------
  IF(dayOfYear <= 182) THEN
   i = CO2ppmv(previous)
   f = CO2ppmv(current)
   n = dayOfYear+183
  ELSE
   i = CO2ppmv(current)
   f = CO2ppmv(next)
   n = dayOfYear-183
  END IF

! Linear interpolation to the given day-of-year.
! ----------------------------------------------
  GetCO2Pert = i + (f-i)*n/365.00

! Convert to mole fraction (volume mixing ratio).
! ----------------------Pert-------------------------
  GetCO2Pert = GetCO2Pert*1.00E-06

 END FUNCTION GetCO2Pert

 end subroutine Run

end module GEOS_RadiationPertGridCompMod
