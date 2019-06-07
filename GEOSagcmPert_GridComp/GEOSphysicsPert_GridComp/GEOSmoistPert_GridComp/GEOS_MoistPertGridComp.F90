! $Id$

! VERIFY_ and RETURN_ macros for error handling.

#include "MAPL_Generic.h"

#define MAPL_FieldBundleGetPointer ESMFL_BundleGetPointerToData

!=============================================================================
!BOP

! !MODULE: GEOS_MoistPert -- A Module to compute linearised moist processes, 
! including convection, large-scale condensation and precipitation and cloud.
! Although the model works with real*4 variables it is found that the moist 
! physics schemes do not satisfy the dot product test of equivalence between 
! the tangent and adjoint in real*4. Precision issues lead to differences 
! between the tlm and adjoint and the tlm at real*4 and the tlm at real*8. 
! By working in real*8 within this module avoids these issues from occuring
! and is found to have little impact on efficiency.

! !INTERFACE:

module GEOS_MoistPertGridCompMod

! GEOS Modules
  use ESMF
  use MAPL_MOD
  use GEOS_Mod
  use GEOS_PertSharedMod

! Saturation table
  use qsat_util

! Moist Schemes
  use CONVECTION
  use CONVECTION_AD
  use CONVECTION_TL
  use CLOUD
  use CLOUD_AD
  use CLOUD_TL

! Pert timers
  use fv_timing_mod,  only: timing_on, timing_off

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION:
! 
!   {\tt GEOS_MoistPertGridCompMod} implements linearised moist processes in GEOS-5. 
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
    call MAPL_TimerAdd(GC,   name="-MOIST"       ,RC=STATUS)
    VERIFY_(STATUS)

! The same run method is registered twice
!  It queries the phase to do TLM or ADJ
! ---------------------------------------

    call MAPL_GridCompSetEntryPoint (gc, ESMF_METHOD_RUN,  Run,        rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint (gc, ESMF_METHOD_RUN,  Run,        rc=status)
    VERIFY_(STATUS)


! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)
     
    RETURN_(ESMF_SUCCESS)
     
  end subroutine SetServices


! Main run component
! ------------------

subroutine Run(gc, import, export, clock, rc)

  !Run arguments
  type(ESMF_GridComp), intent(inout) :: gc
  type (ESMF_State),   intent(inout) :: import
  type (ESMF_State),   intent(inout) :: export
  type (ESMF_Clock),   intent(inout) :: clock
  integer,  optional,  intent(  out) :: rc 

  !Mapl variables
  type (MAPL_MetaComp), pointer      :: MAPL
  integer                            :: STATUS
  character(len=ESMF_MAXSTR)         :: IAm
  character(len=ESMF_MAXSTR)         :: COMP_NAME
  character(len=ESMF_MAXSTR)         :: PHASE_NAME

  !Grid lookup for params
  character(len=ESMF_MAXSTR)         :: GRIDNAME
  character(len=4)                   :: imchar
  character(len=2)                   :: dateline
  integer                            :: imsize, nn
  
  !Control variables
  real(8)                            :: DT
  integer                            :: IM, JM, LM
  integer                            :: PHASE
  integer                            :: I, J, L
  type (ESMF_FieldBundle)            :: Traj3, Traj2, Pert

  !Reference and perturbation trajectories (io single precision)
  real(4), pointer, dimension(:,:,:) :: UT4, VT4, PTT4, QVT4            , QLST4, QCNT4, CFCNT4
  real(4), pointer, dimension(:,:,:) :: UP4, VP4, PTP4, QVP4, QIP4, QLP4              , CFCNP4
  real(4), pointer, dimension(:,:)   :: PS4, KCBL4, TS4, FRLAND4, KHu4, KHl4
  real(4), pointer, dimension(:,:,:) :: TrT4, TrP4

  !Reference and perturbation trajectories (working double precision)
  real(8), allocatable, dimension(:,:,:) :: UT, VT, PTT, QVT
  real(8), allocatable, dimension(:,:,:) :: UP, VP, PTP, QVP            
  real(8), allocatable, dimension(:,:)   :: PS, TS, FRLAND
  integer, allocatable, dimension(:,:)   :: KCBL, KHu, KHl

  !Pressure and temperature variables
  real(8), pointer                       :: ak(:), bk(:)
  real(8), allocatable, dimension(:)     :: pref
  real(8), allocatable, dimension(:,:,:) :: PLE, CNV_PLE, PLO, PK
  real(8), allocatable, dimension(:,:,:) :: TEMP

  !Convection scheme variables
  integer :: ICMIN, MAXCONDEP
  real(8), parameter :: PMIN_DET = 3000.0, AUTOC_CN_OCN  = 2.5e-3, AUTOC_CN_LAND = AUTOC_CN_OCN
  real(8), allocatable, dimension(:,:,:) :: CNV_DQLDTT, CNV_MFDT, CNV_PRC3T, CNV_UPDFT
  real(8), allocatable, dimension(:,:,:) :: CNV_DQLDTP, CNV_MFDP, CNV_PRC3P, CNV_UPDFP
  real(8), allocatable, dimension(:,:,:) :: PTT_C, QVT_C
  real(8), allocatable, dimension(:,:,:) :: CNV_DQLDTT_C, CNV_MFDT_C, CNV_PRC3T_C, CNV_UPDFT_C
  real(8), allocatable, dimension(:,:,:) :: PTT_F, QVT_F
  real(8), allocatable, dimension(:,:,:) :: PTT_L, QVT_L
  integer, allocatable, dimension(:,:)   :: SEEDRAS
  real(8), allocatable, dimension(:,:)   :: CO_AUTO
  real(8), allocatable, dimension(:)     :: SIGE
  real(8), allocatable, dimension(:,:,:) :: WGT0, WGT1
  integer, parameter :: n_rasparams = 25
  real(8) :: RASPARAMS(n_rasparams)
  !Convection filtering
  real(8), allocatable, dimension(:,:,:) :: HEAT
  integer, allocatable, dimension(:,:)   :: DOCONVEC, CTOP
  real(8), allocatable, dimension(:,:)   :: sumHEAT
  real(8), allocatable, dimension(:,:) :: JACOBIAN 
  real(8), allocatable, dimension(:)   :: PT_pert, QV_pert
  real(8), allocatable, dimension(:)   :: PT_pert_in, QV_pert_in
  real(8), allocatable, dimension(:)   :: H_pert, M_pert

  !Cloud scheme variables
  integer, parameter :: N_CLOUD_PARAMS = 57
  real(8) :: CLOUDPARAMS (N_CLOUD_PARAMS)
  real(8), allocatable, dimension(:,:,:) :: QILST, QLLST, QICNT, QLCNT
  real(8), allocatable, dimension(:,:,:) :: QILSP, QLLSP, QICNP, QLCNP
  real(8), allocatable, dimension(:,:,:) :: CFLST, CFCNT
  real(8), allocatable, dimension(:,:,:) :: CFLSP, CFCNP
  real(8), allocatable, dimension(:,:,:) :: fQi
  real(8), allocatable, dimension(:,:,:) :: ILSF, ICNF, LLSF, LCNF

  !R8 versions of the required MAPL_Constants
  real(8) :: MAPL8_CP, MAPL8_ALHL, MAPL8_GRAV, MAPL8_P00, MAPL8_KAPPA
  real(8) :: MAPL8_RGAS, MAPL8_H2OMW, MAPL8_AIRMW, MAPL8_VIREPS
  real(8) :: MAPL8_RUNIV, MAPL8_ALHF, MAPL8_PI, MAPL8_ALHS
  real(8) :: MAPL8_TICE, MAPL8_RVAP

  !Saturation vapour pressure lookup table
  integer, parameter :: DEGSUBS   = 100
  real(8), parameter :: TMINTBL   = 150.0, TMAXTBL = 333.0
  integer, parameter :: TABLESIZE = nint(TMAXTBL-TMINTBL)*DEGSUBS + 1
  real(8)            :: ESTBLX(TABLESIZE)

  !Tracers
  integer                                  :: ITRCR       !Number of tracers
  real(8), allocatable, dimension(:,:,:,:) :: XHOT, XHOP  !Tracer arrays
  real(8), allocatable, dimension(:)       :: FSCAV       !Scavening coefficients

! Begin
!------
    Iam = "Run"
    call ESMF_GridCompGet( GC, name=COMP_NAME, currentPhase=Phase, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the generic state
! -----------------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

! Turn on MAPL timers
! -------------------
    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"-MOIST")
    call MAPL_TimerOn(MAPL,PHASE_NAME)

! Get time step of the linear model
! ---------------------------------
    call MAPL_GetResource(MAPL, DT, Label="RUN_DT:", RC=STATUS)
    VERIFY_(STATUS)

! Set phase
! ---------
    select case(phase)
    case(TLMphase)
       PHASE_NAME = "RUNTLM"
    case(ADJphase)
       PHASE_NAME = "RUNADJ"
    case default
       ASSERT_(.false.)
    end select

! Begin linearised moist physics calculations if asked for in AGCM_apert
! ----------------------------------------------------------------------
    IF (DO_MOIST_PHYS == 1 .or. DO_MOIST_PHYS == 2 .or. DO_MOIST_PHYS == 3 ) THEN

! Get the phase
! -------------
       if (PHASE==ADJphase) then         
          call timing_on('MOIST_ADM')
       elseif (PHASE==TLMphase) then         
          call timing_on('MOIST_TLM')
       endif
 
! Read in reference and perturbation trajectories
! -----------------------------------------------
       call ESMF_StateGet(Import, "TRAJWRK3", Traj3, rc=status)
       VERIFY_(STATUS)
       call ESMF_StateGet(Import, 'TRAJWRK2', Traj2, rc=STATUS)
       VERIFY_(STATUS)
       call ESMF_StateGet(Import, "PERTWRK" , Pert , rc=STATUS)
       VERIFY_(STATUS)

       !U wind speed
       call MAPL_FieldBundleGetPointer(Traj3, "U"     , UT4    , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Pert , "U"     , UP4    , rc=STATUS)
       VERIFY_(STATUS)

       !V wind speed
       call MAPL_FieldBundleGetPointer(Traj3, "V"     , VT4    , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Pert , "V"     , VP4    , rc=STATUS)
       VERIFY_(STATUS)
       
       !Potential temperature
       call MAPL_FieldBundleGetPointer(Traj3, "PT"    , PTT4   , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Pert , "PT"    , PTP4   , rc=STATUS)
       VERIFY_(STATUS)
       
       !Specific humidity
       call MAPL_FieldBundleGetPointer(Traj3, "QV"    , QVT4   , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Pert , "QV"    , QVP4   , rc=STATUS)
       VERIFY_(STATUS)
 
       if (DO_MOIST_PHYS .ne. 3) then
 
          !Total cloud liquid ice and water (pert only)
          call MAPL_FieldBundleGetPointer(Pert , "QI"    , QIP4   , rc=STATUS)
          VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Pert , "QL"    , QLP4   , rc=STATUS)
          VERIFY_(STATUS)
   
          !Cloud totals for large scale and convective (traj only)
          call MAPL_FieldBundleGetPointer(Traj3, "QLS"   , QLST4  , rc=STATUS)
          VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Traj3, "QCN"   , QCNT4  , rc=STATUS)
          VERIFY_(STATUS)
   
          !Anvil fraction
          call MAPL_FieldBundleGetPointer(Traj3, "CFCN"  , CFCNT4 , rc=STATUS)
          VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Pert , "CFCN"  , CFCNP4 , rc=STATUS)
          VERIFY_(STATUS)

       endif

       !Rank 2 variables (traj only)
       call MAPL_FieldBundleGetPointer(Traj2, "PS"     , PS4    , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Traj2, "KCBL"   , KCBL4  , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Traj2, "TS"     , TS4    , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Traj2, "FRLAND" , FRLAND4, rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Traj2, "KHL"    , KHL4   , rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Traj2, "KHU"    , KHU4   , rc=STATUS)
       VERIFY_(STATUS)

! Find grid dimensions for processor
! ----------------------------------
       IM = size(UP4,1)
       JM = size(UP4,2)
       LM = size(UP4,3)

! Allocate memory for all the local variables
! -------------------------------------------
       !Double precision trajectories
       allocate(UT(im,jm,lm),           stat=status); VERIFY_(STATUS)
       allocate(VT(im,jm,lm),           stat=status); VERIFY_(STATUS)
       allocate(PTT(im,jm,lm),          stat=status); VERIFY_(STATUS)
       allocate(QVT(im,jm,lm),          stat=status); VERIFY_(STATUS)
       allocate(UP(im,jm,lm),           stat=status); VERIFY_(STATUS)
       allocate(VP(im,jm,lm),           stat=status); VERIFY_(STATUS)
       allocate(PTP(im,jm,lm),          stat=status); VERIFY_(STATUS)
       allocate(QVP(im,jm,lm),          stat=status); VERIFY_(STATUS)
       allocate(PS(im,jm),              stat=status); VERIFY_(STATUS)
       allocate(TS(im,jm),              stat=status); VERIFY_(STATUS)
       allocate(FRLAND(im,jm),          stat=status); VERIFY_(STATUS)
       allocate(KCBL(im,jm),            stat=status); VERIFY_(STATUS)
       allocate(KHl(IM,JM),             stat=status); VERIFY_(STATUS)
       allocate(KHu(IM,JM),             stat=status); VERIFY_(STATUS)

       !Pressure and temperature variables
       allocate(ak(size(pert_ak)),      stat=status); VERIFY_(STATUS)
       allocate(bk(size(pert_ak)),      stat=status); VERIFY_(STATUS)
       allocate(pref(0:size(AK)-1),     stat=status); VERIFY_(STATUS)
       allocate(PLE(im,jm,0:lm),        stat=status); VERIFY_(STATUS)
       allocate(CNV_PLE(im,jm,0:lm),    stat=status); VERIFY_(STATUS)
       allocate(PLO(im,jm,lm),          stat=status); VERIFY_(STATUS)
       allocate(PK(im,jm,lm),           stat=status); VERIFY_(STATUS)
       allocate(TEMP(IM,JM,LM),         stat=status); VERIFY_(STATUS)

       !Convection varaibles
       allocate(PTT_F(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QVT_F(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(PTT_L(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QVT_L(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(PTT_C(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QVT_C(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(CNV_DQLDTT_C(im,jm,lm), stat=status); VERIFY_(STATUS)
       allocate(CNV_MFDT_C(im,jm,lm),   stat=status); VERIFY_(STATUS)
       allocate(CNV_PRC3T_C(im,jm,lm),  stat=status); VERIFY_(STATUS)
       allocate(CNV_UPDFT_C(im,jm,lm),  stat=status); VERIFY_(STATUS)
       allocate(CNV_DQLDTT(im,jm,lm),   stat=status); VERIFY_(STATUS)
       allocate(CNV_DQLDTP(im,jm,lm),   stat=status); VERIFY_(STATUS)
       allocate(CNV_MFDT(im,jm,lm),     stat=status); VERIFY_(STATUS)
       allocate(CNV_MFDP(im,jm,lm),     stat=status); VERIFY_(STATUS)
       allocate(CNV_PRC3T(im,jm,lm),    stat=status); VERIFY_(STATUS)
       allocate(CNV_PRC3P(im,jm,lm),    stat=status); VERIFY_(STATUS)
       allocate(CNV_UPDFT(im,jm,lm),    stat=status); VERIFY_(STATUS)
       allocate(CNV_UPDFP(im,jm,lm),    stat=status); VERIFY_(STATUS)
       allocate(SEEDRAS(im,jm),         stat=status); VERIFY_(STATUS)
       allocate(SIGE(0:lm),             stat=status); VERIFY_(STATUS)
       allocate(WGT0(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(WGT1(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(CO_AUTO(im,jm),         stat=status); VERIFY_(STATUS)
       allocate(DOCONVEC(im,jm),        stat=status); VERIFY_(STATUS)
       allocate(HEAT(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(CTOP(im,jm),            stat=status); VERIFY_(STATUS)
       allocate(sumHEAT(im,jm),         stat=status); VERIFY_(STATUS)
       allocate(JACOBIAN(2*lm,2),       stat=status); VERIFY_(STATUS)
       allocate(H_pert(lm),             stat=status); VERIFY_(STATUS)
       allocate(M_pert(lm),             stat=status); VERIFY_(STATUS) 
       allocate(PT_pert(lm),            stat=status); VERIFY_(STATUS)
       allocate(QV_pert(lm),            stat=status); VERIFY_(STATUS)
       allocate(PT_pert_in(lm),         stat=status); VERIFY_(STATUS)
       allocate(QV_pert_in(lm),         stat=status); VERIFY_(STATUS)

       !Cloud variables
       allocate(QILST(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QLLST(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QICNT(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QLCNT(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QILSP(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QLLSP(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QICNP(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(QLCNP(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(CFLST(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(CFCNT(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(CFLSP(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(CFCNP(im,jm,lm),        stat=status); VERIFY_(STATUS)
       allocate(fQi(IM,JM,LM),          stat=status); VERIFY_(STATUS)
       allocate(ILSF(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(ICNF(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(LLSF(im,jm,lm),         stat=status); VERIFY_(STATUS)
       allocate(LCNF(im,jm,lm),         stat=status); VERIFY_(STATUS)

! Create double precision versions of required constants and varibles
! -------------------------------------------------------------------
       !MAPL_Constants
       MAPL8_CP     = dble(MAPL_CP)
       MAPL8_ALHL   = dble(MAPL_ALHL)
       MAPL8_GRAV   = dble(MAPL_GRAV)
       MAPL8_P00    = dble(MAPL_P00)
       MAPL8_KAPPA  = dble(MAPL_KAPPA)
       MAPL8_RGAS   = dble(MAPL_RGAS)
       MAPL8_H2OMW  = dble(MAPL_H2OMW)
       MAPL8_AIRMW  = dble(MAPL_AIRMW)
       MAPL8_VIREPS = dble(MAPL_VIREPS)
       MAPL8_RUNIV  = dble(MAPL_RUNIV)
       MAPL8_ALHF   = dble(MAPL_ALHF)
       MAPL8_PI     = dble(MAPL_PI)
       MAPL8_ALHS   = dble(MAPL_ALHS)
       MAPL8_TICE   = dble(MAPL_TICE)
       MAPL8_RVAP   = dble(MAPL_RVAP)

       !Copy reference trajectory to double precision
       UT     = dble(UT4)
       VT     = dble(VT4)
       PTT    = dble(PTT4)*(MAPL8_P00**MAPL8_KAPPA)
       QVT    = dble(QVT4)
       CFLST  = 0.0
       if (DO_MOIST_PHYS .ne. 3) CFCNT  = dble(CFCNT4)
       PS     = dble(PS4)
       TS     = dble(TS4)
       FRLAND = dble(FRLAND4)
       KCBL   = nint(KCBL4)
       KHL    = nint(KHL4)
       KHU    = nint(KHU4)

       !Copy perturbation trajectory to double precision
       UP = dble(UP4)
       VP = dble(VP4)
       if (PHASE==ADJphase) then
          PTP = dble(PTP4)/(MAPL8_P00**MAPL8_KAPPA)
       elseif (PHASE==TLMphase) then 
          PTP = dble(PTP4)*(MAPL8_P00**MAPL8_KAPPA)
       endif
       QVP = dble(QVP4)
       CFLSP = 0.0
       if (DO_MOIST_PHYS .ne. 3) CFCNP = dble(CFCNP4)

! Create additional versions of the trajectory to be overwritten
! --------------------------------------------------------------
       PTT_C = PTT
       QVT_C = QVT
       CNV_DQLDTT_C  = 0.0
       CNV_MFDT_C    = 0.0
       CNV_PRC3T_C   = 0.0
       CNV_UPDFT_C   = 0.0

       PTT_F = PTT
       QVT_F = QVT

       PTT_L = PTT
       QVT_L = QVT

! Define up the constants used by the convection scheme
! -----------------------------------------------------
       !RASPARAMS 23 (MAXDALLOWED) and CLOUDPARAMS 41 (MINRHCRIT) Set based on resolution.
       !The following reads PHYSICS_PERT_GRIDNAME: from AGCM_APERT. E.G. "PHYSICS_PERT_GRIDNAME: PE90x540-CF"
       !90X540 means face of cube is 90 BY 90 (1 degree), 540/90 = 6 sides of the cube)
       !Equivalnet format for lat-long is "GRIDNAME: PC288x181-DC"
       call MAPL_GetResource(MAPL,GRIDNAME,'PHYSICS_PERT_GRIDNAME:', RC=STATUS)
       VERIFY_(STATUS)
       GRIDNAME =  AdjustL(GRIDNAME)
             nn = len_trim(GRIDNAME)
       dateline = GRIDNAME(nn-1:nn)
       imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
       read(imchar,*) imsize
       if(dateline.eq.'CF') imsize = imsize*4

       RASPARAMS( 1) = 1.000
       RASPARAMS( 2) = 0.05
       RASPARAMS( 3) = 0.0   ! NOW IN CO_AUTO (CONTAINED WITHIN SUBROUTINE)
       RASPARAMS( 4) = 8.0e-4
       RASPARAMS( 5) = 1800.
       RASPARAMS( 6) = 43200.0
       RASPARAMS( 7) = -300. !RASNCL, CONTROLS FINDDTLS, USE OF RANDOM NUMBER
       RASPARAMS( 8) = 4.0 
       RASPARAMS( 9) = 0.0 
       RASPARAMS(10) = 200.
       RASPARAMS(11) = 7.5e-4
       RASPARAMS(12) = 1.0 
       RASPARAMS(13) =-1.0 
       RASPARAMS(14) = 1.3 
       RASPARAMS(15) = 1.3 
       RASPARAMS(16) = 263.
       RASPARAMS(17) = 0.5 
       RASPARAMS(18) = 1.0 
       RASPARAMS(19) = 0.0 
       RASPARAMS(20) = 0.1 
       RASPARAMS(21) = 0.8 
       RASPARAMS(22) = 1.0
       if( imsize .le. 200                      ) RASPARAMS(23) = 4000.0
       if( imsize .gt. 200 .and. imsize.le.400  ) RASPARAMS(23) = 2000.0
       if( imsize .gt. 400 .and. imsize.le.800  ) RASPARAMS(23) = 700.0
       if( imsize .gt. 800 .and. imsize.le.1600 ) RASPARAMS(23) = 450.0
       if( imsize .gt. 1600                     ) RASPARAMS(23) = 450.0
       RASPARAMS(24) = 0.5
       RASPARAMS(25) = 0.65

! Define constants that are used by the cloud scheme
! --------------------------------------------------
       !SET SBAC PARAMETERS
       CLOUDPARAMS( 1) = 10.0
       CLOUDPARAMS( 2) = 4.0 
       CLOUDPARAMS( 3) = 4.0
       CLOUDPARAMS( 4) = 1.0
       CLOUDPARAMS( 5) = 2.0e-3
       CLOUDPARAMS( 6) = 8.0e-4
       CLOUDPARAMS( 7) = 2.0   
       CLOUDPARAMS( 8) = 1.0  
       CLOUDPARAMS( 9) = -1.0  
       CLOUDPARAMS(10) = 0.0   
       CLOUDPARAMS(11) = 1.3   
       CLOUDPARAMS(12) = 1.0e-9
       CLOUDPARAMS(13) = 3.3e-4
       CLOUDPARAMS(14) = 20.   
       CLOUDPARAMS(15) = 4.8   
       CLOUDPARAMS(16) = 4.8   
       CLOUDPARAMS(17) = 230.  
       CLOUDPARAMS(18) = 1.0   
       CLOUDPARAMS(19) = 1.0   
       CLOUDPARAMS(20) = 230.  
       CLOUDPARAMS(21) = 14400.
       CLOUDPARAMS(22) = 50.   
       CLOUDPARAMS(23) = 0.01  
       CLOUDPARAMS(24) = 0.1   
       CLOUDPARAMS(25) = 200.  
       CLOUDPARAMS(26) = 0.    
       CLOUDPARAMS(27) = 0.    
       CLOUDPARAMS(28) = 0.5   
       CLOUDPARAMS(29) = 0.5   
       CLOUDPARAMS(30) = 2000. 
       CLOUDPARAMS(31) = 0.8   
       CLOUDPARAMS(32) = 0.5   
       CLOUDPARAMS(33) = -40.0 
       CLOUDPARAMS(34) = 1.0 
       CLOUDPARAMS(35) = 4.0 
       CLOUDPARAMS(36) = 0.0 
       CLOUDPARAMS(37) = 0.0 
       CLOUDPARAMS(38) = 0.0 
       CLOUDPARAMS(39) = 1.0e-3
       CLOUDPARAMS(40) = 8.0e-4
       CLOUDPARAMS(41) = 1.0   
       if( imsize .le. 200                      ) CLOUDPARAMS(42) = 0.80
       if( imsize .gt. 200 .and. imsize.le.400  ) CLOUDPARAMS(42) = 0.90
       if( imsize .gt. 400 .and. imsize.le.800  ) CLOUDPARAMS(42) = 0.93
       if( imsize .gt. 800 .and. imsize.le.1600 ) CLOUDPARAMS(42) = 0.95
       if( imsize .gt. 1600                     ) CLOUDPARAMS(42) = 0.97
       CLOUDPARAMS(43) = 1.0
       CLOUDPARAMS(44) = 0.0
       CLOUDPARAMS(45) = 750.0 
       CLOUDPARAMS(46) = CLOUDPARAMS(42)+0.01 
       CLOUDPARAMS(47) = 1.0
       CLOUDPARAMS(48) = 1.0
       CLOUDPARAMS(49) = 0.0
       CLOUDPARAMS(50) = 0.0
       CLOUDPARAMS(51) = 10.e-6
       CLOUDPARAMS(52) = 20.e-6
       CLOUDPARAMS(53) = 21.e-6
       CLOUDPARAMS(54) = 40.e-6
       CLOUDPARAMS(55) = 30.e-6
       CLOUDPARAMS(56) = 1.0  
       CLOUDPARAMS(57) = 1.0  

! Generatre saturation vapour pressure loopup table
! -------------------------------------------------
       call ESINIT(ESTBLX)

! Compute reference pressure, pressure, Exner pressure and temperature
! -------------------------------------------------------------------
       !Compute pref from ak and bk
       ak = pert_ak
       bk = pert_bk
       DO L = 0,lm
          pref(L) = AK(L+1) + BK(L+1)*MAPL8_P00
       enddo

       !Pressure at the half levels from Ps
       DO L = 0,lm
          PLE(:,:,L) = AK(L+1) + BK(L+1)*PS(:,:)
       enddo

       !Pressure in hPa, Exner pressure and temperature
       CNV_PLE  = 0.01*PLE
       PLO      = 0.5*(CNV_PLE(:,:,0:LM-1) +  CNV_PLE(:,:,1:LM  ) )
       PK       = (PLO/1000.)**(MAPL8_RGAS/MAPL8_CP)
       TEMP     = PTT*PK

! Prepare convection inputs and call nonlinear convection scheme 
! --------------------------------------------------------------
       ICMIN = max(1,count(PREF < PMIN_DET))
       SIGE = PREF/PREF(LM)

       !Not linearised as could produce unpredictable behaviour (and very sensitive).
       SEEDRAS(:,:) = 1000000 * ( 100*TEMP(:,:,LM) - INT( 100*TEMP(:,:,LM) ) )

       !Strapping levels
       DO I = 1,IM
          DO J = 1,JM
             WGT0(I,J,:)            = 0.0
             WGT0(I,J,KCBL(I,J):LM) = 1.0 
             WGT1(I,J,:)            = 0.0
             WGT1(I,J,KCBL(I,J):LM) = 1.0 
          ENDDO
       ENDDO

       where (FRLAND<0.1) 
          CO_AUTO = AUTOC_CN_OCN   ! ocean value
       elsewhere
          CO_AUTO = AUTOC_CN_LAND  ! land value
       end where

       if (DO_MOIST_PHYS == 3) then

          !This just performs the tracer scavening.

          !call MAPL_FieldBundleGetPointer(Traj3, "TR" , TrT4 , rc=STATUS)
          !VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Pert , "O3" , TrP4 , rc=STATUS)
          VERIFY_(STATUS)

          !Set the number of tracers
          ITRCR = 1

          allocate(XHOT(im,jm,lm,ITRCR),  stat=status); VERIFY_(STATUS)
          allocate(XHOP(im,jm,lm,ITRCR),  stat=status); VERIFY_(STATUS)
          allocate(FSCAV(ITRCR),          stat=status); VERIFY_(STATUS)

          !Copy tracers into double precision holders
          XHOT(:,:,:,1) = 0.0 !dble(TrT4)
          XHOP(:,:,:,1) = dble(TrP4)
             
          !CO2 Value
          FSCAV(1) = 0.0

          !Call TLM/ADM for the tracer only version of RAS
          if (PHASE==ADJphase) then

             DO I = 1,IM
                DO J = 1,JM 

                   !ADJOINT GOES HERE
                   call RASE_TRACER_B( 1, 1, LM, ICMIN, DT,                              &
                                       MAPL8_CP, MAPL8_ALHL, MAPL8_GRAV, MAPL8_RGAS,     &
                                       MAPL8_H2OMW, MAPL8_AIRMW, MAPL8_VIREPS,           &
                                       SEEDRAS(I,J), SIGE,                               &
                                       KCBL(I,J),                                        &
                                       WGT0(I,J,:), WGT1(I,J,:),                         &
                                       FRLAND(I,J), Ts(I,J),                             &
                                       PTT(I,J,:),                                       &
                                       QVT(I,J,:),                                       &
                                       UT(I,J,:),                                        &
                                       VT(I,J,:),                                        &
                                       CO_AUTO(I,J), CNV_PLE(I,J,:),                     &
                                       RASPARAMS, ESTBLX,                                &
                                       ITRCR,                                            &
                                       XHOT(I,J,:,:), XHOP(I,J,:,:),                     &
                                       FSCAV                                             )


                endDO
             endDO

          elseif (PHASE==TLMphase) then

             DO I = 1,IM
                DO J = 1,JM 

                   call RASE_TRACER_D( 1, 1, LM, ICMIN, DT,                              &
                                       MAPL8_CP, MAPL8_ALHL, MAPL8_GRAV, MAPL8_RGAS,     &
                                       MAPL8_H2OMW, MAPL8_AIRMW, MAPL8_VIREPS,           &
                                       SEEDRAS(I,J), SIGE,                               &
                                       KCBL(I,J),                                        &
                                       WGT0(I,J,:), WGT1(I,J,:),                         &
                                       FRLAND(I,J), Ts(I,J),                             &
                                       PTT(I,J,:),                                       &
                                       QVT(I,J,:),                                       &
                                       UT(I,J,:),                                        &
                                       VT(I,J,:),                                        &
                                       CO_AUTO(I,J), CNV_PLE(I,J,:),                     &
                                       RASPARAMS, ESTBLX,                                &
                                       ITRCR,                                            &
                                       XHOT(I,J,:,:), XHOP(I,J,:,:),                     &
                                       FSCAV                                             )

                endDO
             endDO


          endif

          !Move double precision tracer perturbation back to single precision holder (not trajectory)
          TrP4 = real(XHOP(:,:,:,1))

          !Deallocate memory
          deallocate(XHOT,XHOP,FSCAV)

       else

          !Call nonlinear convection scheme. We need to do this because the
          !cloud adjoint scheme needs the outputs from the convection.
          !We will also only call the convective linearizations for profiles
          !where convection is occuring, which is a big time saver.
          CALL RASE0(IM*JM, IM*JM, LM, ICMIN, DT,                   &
                     MAPL8_CP, MAPL8_ALHL, MAPL8_GRAV, MAPL8_RGAS,  &
                     MAPL8_H2OMW, MAPL8_AIRMW, MAPL8_VIREPS,        &
                     SEEDRAS, SIGE,                                 &
                     KCBL, WGT0, WGT1, FRLAND, TS,                  &
                     PTT_C, QVT_C,                                  &
                     CO_AUTO, CNV_PLE,                              &
                     CNV_DQLDTT_C, CNV_MFDT_C,                      &
                     CNV_PRC3T_C, CNV_UPDFT_C,                      &
                     RASPARAMS, ESTBLX                              )
   
! Do the filtering to determine whether linear convection should be called
! ------------------------------------------------------------------------
          !Figure out whether or not each convective profile should be linearized.
          !Conditions  - Convection is happening (heating rate is nonzero)
          !            - Convection is deep enough (>= 10 levels)
          !            - The heating rate profile is not just a spike at one level
          !            - Gradients are not too steep (Jacobian filtering).
          DOCONVEC = 0
          HEAT = 0.0
          CTOP = LM
          sumHEAT = 0.0
   
          !Be more lenient on profiles let in for 4DVAR, less likely to encounter problems
          !at shorter lead times and gives more realistic low level cloud perturbations.
          if (DO_MOIST_PHYS == 1) then
             MAXCONDEP = 1
          elseif (DO_MOIST_PHYS == 2) then
             MAXCONDEP = 10
          endif
   
          DO I = 1,IM
             DO J = 1,JM
   
                !Compute the heating rate.
                HEAT(I,J,:) = (PTT_C(I,J,:) - PTT(I,J,:))/DT
   
                !Starting at the top scan downwards to look for nonzero heating rate,
                !record index of highest level convection reaches (ignoring v.small heating rate)
                DO L = 1,LM
                   IF ( abs(HEAT(I,J,L)) .gt. 0.01*maxval(abs(HEAT(I,J,:)),1) ) THEN
                      CTOP(I,J) = L
                      exit
                   ENDIF
                ENDDO
                
                !Compute sort of integral of the heating rate.
                if ( (CTOP(I,J) .ne. LM) .and. (KCBL(I,J) - CTOP(I,J) > 0) ) then
                   sumHEAT(I,J) = (   sum(abs(HEAT(I,J,CTOP(I,J):KCBL(I,J)-1))) &
                                          - maxval(abs(HEAT(I,J,CTOP(I,J):KCBL(I,J)-1)),1) ) &
                                          / ( KCBL(I,J) - CTOP(I,J)  )
                endif
   
                !Compare `integral` to maximum absolute heating rate
                IF ( KCBL(I,J) - CTOP(I,J) >= MAXCONDEP ) THEN
                   IF (sumHEAT(I,J) / maxval(abs(HEAT(I,J,1:KCBL(I,J)-1)),1) > 0.125 ) then
                      DOCONVEC(I,J) = 1
                   endif
                endif
   
                !Compute two columns of the Jacobian to check for steep gradients
                !This prevents instability and floating point issues that can cause failure of the TLM/ADJ dot prod test
                IF ( DOCONVEC(I,J) == 1 ) THEN
                   call jacobian_filter_tlm
                ENDIF
   
             enddo
          enddo
   
! Prepare the inputs for the cloud scheme 
! ---------------------------------------
          !Compute the ice fraction for each grid box
          DO i = 1,IM
             DO j = 1,JM
                DO l = 1,LM
                   call IceFraction( TEMP(i,j,l), fQi(i,j,l) )
                enddo
             enddo
          enddo
   
          !Split the input large scale and convective cloud into ice and liquid parts.
          QILST = QLST4 * fQi
          QLLST = QLST4 * (1-fQi)
          QICNT = QCNT4 * fQi
          QLCNT = QCNT4 * (1-fQi)
   
          !Split the perturbations for total cloud water and ice into the convective and large scale parts.
          !Spilitting is based on fraction of total cloud ice/water that is attributed to large scale and anvil
          !in the trajectory at this time step. Note that we don't really need to protect for small values of
          !total cloud sum, if the denominator is small then so is the numerator, though compiler may complain.
          ILSF = 0.0
          ICNF = 0.0
          LLSF = 0.0
          LCNF = 0.0
          DO I = 1,IM
             DO J = 1,JM
                DO L = 1,LM
    
                   if ( QILST(i,j,l) + QICNT(i,j,l) .gt. 0.0 ) then
                      ILSF(i,j,l) = QILST(i,j,l) / ( QILST(i,j,l) + QICNT(i,j,l) )
                      ICNF(i,j,l) = QICNT(i,j,l) / ( QILST(i,j,l) + QICNT(i,j,l) )
                   endif
                   if ( QLLST(i,j,l) + QLCNT(i,j,l) .gt. 0.0 ) then
                      LLSF(i,j,l) = QLLST(i,j,l) / ( QLLST(i,j,l) + QLCNT(i,j,l) )
                      LCNF(i,j,l) = QLCNT(i,j,l) / ( QLLST(i,j,l) + QLCNT(i,j,l) )
                   endif
   
                enddo
             enddo
          enddo
   
          !Use fraction to split the perturbation total liquid ice and water.
          if (PHASE==ADJphase) then
             !Adjoint of summing to QIT/QLT
             QILSP = QIP4
             QICNP = QIP4
             QLLSP = QLP4
             QLCNP = QLP4
          elseif (PHASE==TLMphase) then
             !Splitting to large scale and convective parts
             QILSP = QIP4 * ILSF
             QICNP = QIP4 * ICNF
             QLLSP = QLP4 * LLSF
             QLCNP = QLP4 * LCNF
          endif
   
! Reset/initialise the RAS outputs to zero
! -------------------------------------------
          CNV_DQLDTT  = 0.0
          CNV_DQLDTP  = 0.0
          CNV_MFDT    = 0.0
          CNV_MFDP    = 0.0
          CNV_PRC3T   = 0.0
          CNV_PRC3P   = 0.0
          CNV_UPDFT   = 0.0
          CNV_UPDFP   = 0.0
   
! Perform tangent linear and adjoint computations of the main convection and cloud schemes
! ----------------------------------------------------------------------------------------
          if (PHASE==ADJphase) then ! Do adjoint
   
             !Call the adjoint of the cloud scheme
             CALL CLOUD_DRIVER_B ( DT, IM, JM, LM, PTT_C, PTP, QVT_C, QVP, PLE,        &
                                   CNV_DQLDTT_C, CNV_DQLDTP, CNV_MFDT_C, CNV_MFDP,     &
                                   CNV_PRC3T_C, CNV_PRC3P, CNV_UPDFT_C, CNV_UPDFP,     &
                                   QILST, QILSP, QLLST, QLLSP,                         &
                                   QICNT, QICNP, QLCNT, QLCNP,                         &
                                   CFLST, CFLSP, CFCNT, CFCNP,                         &
                                   FRLAND, CLOUDPARAMS, ESTBLX, KHu, KHl,              &
                                   MAPL8_RUNIV, MAPL8_KAPPA, MAPL8_airmw, MAPL8_h2omw, &
                                   MAPL8_GRAV, MAPL8_ALHL, MAPL8_ALHF, MAPL8_PI,       &
                                   MAPL8_RGAS, MAPL8_CP, MAPL8_VIREPS, MAPL8_ALHS,     &
                                   MAPL8_TICE, MAPL8_RVAP, MAPL8_P00, DO_MOIST_PHYS    )
   
             !Call the adjoint of the convection scheme
             DO I = 1,IM
                DO J = 1,JM
                   if ( DOCONVEC(I,J) == 1) then
   
                      call RASE_B (1, 1, LM, ICMIN, DT,                              &
                                   MAPL8_CP, MAPL8_ALHL, MAPL8_GRAV, MAPL8_RGAS,     &
                                   MAPL8_H2OMW, MAPL8_AIRMW, MAPL8_VIREPS,           &
                                   SEEDRAS(I,J), SIGE,                               &
                                   KCBL(I,J),                                        &
                                   WGT0(I,J,:), WGT1(I,J,:),                         &
                                   FRLAND(I,J), Ts(I,J),                             &
                                   PTT_L(I,J,:), PTP(I,J,:),                         &
                                   QVT_L(I,J,:), QVP(I,J,:),                         &
                                   UT(I,J,:), UP(I,J,:),                             &
                                   VT(I,J,:), VP(I,J,:),                             &
                                   CO_AUTO(I,J), CNV_PLE(I,J,:),                     &
                                   CNV_DQLDTT(I,J,:), CNV_DQLDTP(I,J,:),             &
                                   CNV_MFDT(I,J,:),   CNV_MFDP(I,J,:),               &
                                   CNV_PRC3T(I,J,:),  CNV_PRC3P(I,J,:),              &
                                   CNV_UPDFT(I,J,:),  CNV_UPDFP(I,J,:),              &
                                   RASPARAMS, ESTBLX                                 )
   
                   endif
                enddo
             enddo
          
          else
   
             !Call the tangent linear convection scheme
             DO I = 1,IM
                DO J = 1,JM
                   if ( DOCONVEC(I,J) == 1) then
   
                      call RASE_D (1, 1, LM, ICMIN, DT,                              &
                                   MAPL8_CP, MAPL8_ALHL, MAPL8_GRAV, MAPL8_RGAS,     &
                                   MAPL8_H2OMW, MAPL8_AIRMW, MAPL8_VIREPS,           &
                                   SEEDRAS(I,J), SIGE,                               &
                                   KCBL(I,J),                                        &
                                   WGT0(I,J,:), WGT1(I,J,:),                         &
                                   FRLAND(I,J), Ts(I,J),                             &
                                   PTT(I,J,:), PTP(I,J,:),                           &
                                   QVT(I,J,:), QVP(I,J,:),                           &
                                   UT(I,J,:), UP(I,J,:),                             &
                                   VT(I,J,:), VP(I,J,:),                             &
                                   CO_AUTO(I,J), CNV_PLE(I,J,:),                     &
                                   CNV_DQLDTT(I,J,:), CNV_DQLDTP(I,J,:),             &
                                   CNV_MFDT(I,J,:),   CNV_MFDP(I,J,:),               &
                                   CNV_PRC3T(I,J,:),  CNV_PRC3P(I,J,:),              &
                                   CNV_UPDFT(I,J,:),  CNV_UPDFP(I,J,:),              &
                                   RASPARAMS, ESTBLX                                 )
   
                   endif
                enddo
             enddo
   
             !Call the tangent linear cloud scheme.
             CALL CLOUD_DRIVER_D ( DT, IM, JM, LM, PTT_C, PTP, QVT_C, QVP, PLE,        &
                                   CNV_DQLDTT_C, CNV_DQLDTP, CNV_MFDT_C, CNV_MFDP,     &
                                   CNV_PRC3T_C, CNV_PRC3P, CNV_UPDFT_C, CNV_UPDFP,     &
                                   QILST, QILSP, QLLST, QLLSP,                         &
                                   QICNT, QICNP, QLCNT, QLCNP,                         &
                                   CFLST, CFLSP, CFCNT, CFCNP,                         &
                                   FRLAND, CLOUDPARAMS, ESTBLX, KHu, KHl,              &
                                   MAPL8_RUNIV, MAPL8_KAPPA, MAPL8_airmw, MAPL8_h2omw, &
                                   MAPL8_GRAV, MAPL8_ALHL, MAPL8_ALHF, MAPL8_PI,       &
                                   MAPL8_RGAS, MAPL8_CP, MAPL8_VIREPS, MAPL8_ALHS,     &
                                   MAPL8_TICE, MAPL8_RVAP, MAPL8_P00, DO_MOIST_PHYS    )
   
          endif 
   
! Move perturbation tractory back into the original single precision containers
! -----------------------------------------------------------------------------
          UP4 = real(UP,4)
          VP4 = real(VP,4)
          if (PHASE==ADJphase) then
             PTP4 = real(PTP*(MAPL8_P00**MAPL8_KAPPA),4)
          elseif (PHASE==TLMphase) then 
             PTP4 = real(PTP/(MAPL8_P00**MAPL8_KAPPA),4)
          endif
          QVP4 = real(QVP,4)
          if (PHASE==ADJphase) then
             QIP4 = real(QILSP*ILSF + QICNP*ICNF,4)
             QLP4 = real(QLLSP*LLSF + QLCNP*LCNF,4)
          elseif (PHASE==TLMphase) then
             QIP4 = real(QILSP + QICNP,4)
             QLP4 = real(QLLSP + QLCNP,4)
          endif
          CFCNP4 = real(CFCNP,4)

       endif !Tracer only or full moist physics package.

! Deallocate memory used to store all the local arrays
! ----------------------------------------------------
       deallocate(UT,VT,PTT,QVT,UP,VP,PTP,QVP,PS,TS,FRLAND,KCBL,KHl,KHu)
       deallocate(ak,bk,pref,PLE,CNV_PLE,PLO,PK,TEMP)
       deallocate(PTT_F,QVT_F,PTT_L,QVT_L,PTT_C,QVT_C)
       deallocate(CNV_DQLDTT_C,CNV_MFDT_C,CNV_PRC3T_C,CNV_UPDFT_C)
       deallocate(CNV_DQLDTT,CNV_DQLDTP,CNV_MFDT,CNV_MFDP,CNV_PRC3T,CNV_PRC3P,CNV_UPDFT,CNV_UPDFP)
       deallocate(SEEDRAS,SIGE,WGT0,WGT1,CO_AUTO,DOCONVEC,HEAT,CTOP,sumHEAT)
       deallocate(JACOBIAN,H_pert,M_pert,PT_pert,QV_pert,PT_pert_in,QV_pert_in)
       deallocate(QILST,QLLST,QICNT,QLCNT,QILSP,QLLSP,QICNP,QLCNP)
       deallocate(CFLST,CFCNT,CFLSP,CFCNP)
       deallocate(fQi,ILSF,ICNF,LLSF,LCNF)

! Turn of the timers
! ------------------
       if (PHASE==ADJphase) then         
          call timing_off('MOIST_ADM')
       elseif (PHASE==TLMphase) then         
          call timing_off('MOIST_TLM')
       endif

    ENDIF !DO_MOIST_PHYS

! All done
! --------
    call MAPL_TimerOff(MAPL,PHASE_NAME) 
    call MAPL_TimerOff(MAPL,"TOTAL")
    call MAPL_TimerOff(MAPL,"-MOIST")

    RETURN_(ESMF_SUCCESS)

 contains


! Subroutines
! -----------

 subroutine jacobian_filter_tlm

 !Compute two specific columns of the Jacobian for a single profile. Then filter if values in those
 !columns are over a certain value, implying too steep gradients.

  JACOBIAN = 0.0
  
  DO L = 1,1 !Perturb level by level
 
     PT_pert = 0.0
     QV_pert = 0.0
 
     if (L == 1) then
 
        PT_pert(KCBL(I,J)) = 1.0 
 
     elseif (L == 2) then
 
        if (KCBL(I,J) == lm) then
          QV_pert(KCBL(I,J)) = 1.0 
        else
           QV_pert(KCBL(I,J) + 1) = 1.0 
        endif
 
     endif
 
     !SAVE PRECALL PROGNOSTIC VARIABLES
     PT_pert_in = PT_pert
     QV_pert_in = QV_pert
         
     !DO NONLINEAR CONVECTION TO CREATE INPUTS TRAJECTORY FOR LARGE SCALE ADJOINT
     call RASE0_d ( 1, 1, LM, ICMIN, DT,                           &
                    MAPL8_CP, MAPL8_ALHL, MAPL8_GRAV, MAPL8_RGAS,  &
                    MAPL8_H2OMW, MAPL8_AIRMW, MAPL8_VIREPS,        &
                    SEEDRAS(I,J), SIGE,                      &
                    KCBL(I,J),                               &
                    WGT0(I,J,:), WGT1(I,J,:),          &
                    FRLAND(I,J), Ts(I,J),              &
                    PTT_F(I,J,:), PT_PERT,                   &
                    QVT_F(I,J,:), QV_PERT,                   &
                    CO_AUTO(I,J), CNV_PLE(I,J,:),      &
                    RASPARAMS, ESTBLX                              )
 
     !COMPUTE PERTURBATION HEATING AND MOISTENING RATES
     H_pert = (PT_pert - PT_pert_in)/DT
     M_pert = (QV_pert - QV_pert_in)/DT
       
     !Uncomment here if just doing two columns of the Jacobian
     if (L == 1) then
        JACOBIAN(0*LM+1:1*LM,1) = H_pert
        JACOBIAN(1*LM+1:2*LM,1) = M_pert
     elseif (L == 2) then
        JACOBIAN(0*LM+1:1*LM,2) = H_pert
        JACOBIAN(1*LM+1:2*LM,2) = M_pert
     endif
 
  endDO

  !Constants here determined so as to remove as many of the problematic points as possible.
  !The constants used in this if loop are tuned from looking at many Jacobians for many time steps. Values are choosen
  !so as to balance between keeping the natural behahiour for as many points as possible without 
  IF ( (maxval(abs(Jacobian(1:lm     ,1))) .gt. 0.00010  ) .or. &
       (maxval(abs(Jacobian(1:lm     ,2))) .gt. 0.25     ) .or. & 
       (maxval(abs(Jacobian(lm+1:2*lm,1))) .gt. 1.0e-07  ) .or. &
       (maxval(abs(Jacobian(lm+1:2*lm,2))) .gt. 0.000250 ) ) then

     DOCONVEC(I,J) = 0

  else

     DOCONVEC(I,J) = 1

  endIF 

 endsubroutine jacobian_filter_tlm

 end subroutine Run

end module GEOS_MoistPertGridCompMod


