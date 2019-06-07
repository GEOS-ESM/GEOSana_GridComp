!   $Id$

#include "MAPL_Generic.h"

#define MAPL_FieldBundleGetPointer ESMFL_BundleGetPointerToData
!=============================================================================

module GEOS_TurbulencePertGridCompMod

!BOP

!  !MODULE: GEOS_Turbulence --- An GEOS generic atmospheric turbulence component

! !USES:

  use ESMF
  use GEOS_Mod
  use GEOS_PertSharedMod

  use BLDRIVER

  use fv_timing_mod,  only: timing_on, timing_off

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

! !DESCRIPTION:
! 
! Linearised version of the TurbulenceGridComp.F90 module. Two options exist,
! triggered by DO_BL_PHYS in AGCM_apert.rc. Under option 0 nothing is done 
! here, otherwise a linearisation of the full nonlinear GEOS-5
! Louis/Lock/Belajaars scheme is enacted, where different Ks are used for wind,
! temperature and tracers. Perturbations of the turbulence 
! parameters Km and Kh are neglected. In doing so the turbulence problem
! is solved exactly as it is solved in the full nonlinear model, i.e.
! by solving the tri-diagonal backward implicit finite difference.

! REVISION HISTORY:
! Nov2013 D.Holdaway - Added the BL_OPTION: 2 full boundary layer linearisation.
! Mar2016 D.Holdaway - Updated solver to work for tracers, removed option 1 as obsolete

contains

!=============================================================================
!=============================================================================
!=============================================================================
!=============================================================================

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !DESCRIPTION: This version uses the {\tt GEOS\_GenericSetServices}, which sets
!               the Initialize and Finalize services to generic versions. It also
!               allocates our instance of a generic state and puts it in the 
!               gridded component (GC). Here we only set the two-stage run method and
!               declare the data services.
! 
! !REVISION HISTORY: 
!   ??Jul2006 E.Novak./Todling - Added output defining TLM/ADM trajectory

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code
!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! The same run method is registered twice
!  It queries the phase to do TLM or ADJ
! ---------------------------------------

    call MAPL_GridCompSetEntryPoint (gc, ESMF_METHOD_RUN, Run, rc=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint (gc, ESMF_METHOD_RUN, Run, rc=status)
    VERIFY_(STATUS)

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,   name="RUNTLM"       ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="RUNADJ"       ,RC=STATUS)
    VERIFY_(STATUS)
    
! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices


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

  type(ESMF_GridComp), intent(inout)  :: gc
  type (ESMF_State),   intent(inout)  :: import
  type (ESMF_State),   intent(inout)  :: export
  type (ESMF_Clock),   intent(inout)  :: clock
  integer,  optional,  intent(  out)  :: rc 

!EOP

! !Local Variables:
  
  type (MAPL_MetaComp),      pointer  :: MAPL 
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: IAm
  character(len=ESMF_MAXSTR)          :: COMP_NAME
  character(len=ESMF_MAXSTR)          :: PHASE_NAME
  integer                             :: PHASE
  integer                             :: IM, JM, LM
  integer                             :: I, J, L

  type (ESMF_FieldBundle)             :: Traj3, Traj2, Pert
  real, pointer                       :: UDT(:,:,:), UDP(:,:,:)
  real, pointer                       :: VDT(:,:,:), VDP(:,:,:)
  real, pointer                       :: PTT(:,:,:), PTP(:,:,:)
  real, pointer                       :: QVT(:,:,:), QVP(:,:,:)
  real, pointer                       :: QLT(:,:,:), QLP(:,:,:)
  real, pointer                       :: QIT(:,:,:), QIP(:,:,:)
  real, pointer                       ::             O3P(:,:,:)
!  real, pointer                       ::             TrP(:,:,:)

  !Winds used by the physics
  real, allocatable, dimension(:,:,:) :: UT, VT
  real, allocatable, dimension(:,:,:) :: UP, VP

  !Cloud variables - if doing moist physics
  real, pointer                       :: QLST(:,:,:), QCNT(:,:,:)
  real, allocatable, dimension(:,:,:) :: Ph, PIh, TEMP, fQi

  real, allocatable                   :: QIT1(:,:,:), QLT1(:,:,:)
  real, allocatable                   :: PrT(:,:,:), DPT_PS(:,:,:)

  real, pointer                       :: ak(:), bk(:)
  real, allocatable, dimension(:)     :: pref
  real, pointer, dimension(:,:)       :: PS, FROCEAN, FRLAND, VARFLT
  real, pointer, dimension(:,:)       :: USTAR, BSTAR, ZPBL, CM, CT, CQ

  real, dimension(22)   :: TURBPARAMS
  integer, dimension(4) :: TURBPARAMSI

  real                                :: DT

  real, allocatable, dimension(:,:,:) :: AKQ, BKQ, CKQ, TKQ
  real, allocatable, dimension(:,:,:) :: AKS, BKS, CKS, TKS
  real, allocatable, dimension(:,:,:) :: AKV, BKV, CKV, TKV, EKV, FKV

  real, allocatable, dimension(:,:,:) :: pkt, pet, pke, bx, cx 

  real, allocatable, dimension(:,:,:) :: PTT1
  real, allocatable, dimension(:,:)   :: ZPBL1, CT1


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

    !Get physics time step
    call MAPL_GetResource(MAPL, DT, Label="RUN_DT:", RC=STATUS)
    VERIFY_(STATUS)

    if (DO_BL_PHYS .ne. 0) then

       !Get TLM/ADM Phase
       select case(phase)
       case(TLMphase)
          PHASE_NAME = "RUNTLM"
       case(ADJphase)
          PHASE_NAME = "RUNADJ"
       case default
          ASSERT_(.false.)
       end select

       !Timers on
       call MAPL_TimerOn(MAPL,"TOTAL")
       call MAPL_TimerOn(MAPL,PHASE_NAME)
       if (PHASE==ADJphase) then  
          call timing_on('BDLAY_ADM')
       elseif (PHASE==TLMphase) then
          call timing_on('BDLAY_TLM')
       endif
      
       !Pointer to trajectory
       call ESMF_StateGet(import, "TRAJWRK3", Traj3, rc=status)
       VERIFY_(STATUS)
       call ESMF_StateGet(Import, "PERTWRK",  Pert,  rc=STATUS)
       VERIFY_(STATUS)
       call ESMF_StateGet(Import, 'TRAJWRK2', Traj2, rc=STATUS)
       VERIFY_(STATUS)
   
       !Fill rank 3 pointers
       call MAPL_FieldBundleGetPointer(Traj3, "U" , UDT, rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Pert,  "U" , UDP, rc=STATUS)
       VERIFY_(STATUS)
   
       call MAPL_FieldBundleGetPointer(Traj3, "V" , VDT, rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Pert,  "V" , VDP, rc=STATUS)
       VERIFY_(STATUS)
   
       call MAPL_FieldBundleGetPointer(Traj3, "PT", PTT, rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Pert,  "PT", PTP, rc=STATUS)
       VERIFY_(STATUS)
   
       call MAPL_FieldBundleGetPointer(Pert,  "O3", O3P, rc=STATUS)
       VERIFY_(STATUS)
   
 !      call MAPL_FieldBundleGetPointer(Pert,  "TR", TrP, rc=STATUS)
 !      VERIFY_(STATUS)

       call MAPL_FieldBundleGetPointer(Traj3, "QV", QVT, rc=STATUS)
       VERIFY_(STATUS)
       call MAPL_FieldBundleGetPointer(Pert,  "QV", QVP, rc=STATUS)
       VERIFY_(STATUS)
   
       call MAPL_FieldBundleGetPointer(Pert,  "QL", QLP, rc=STATUS)
       VERIFY_(STATUS)
       if (DO_MOIST_PHYS .eq. 0) then
          call MAPL_FieldBundleGetPointer(Traj3,  "QL", QLT, rc=STATUS)
          VERIFY_(STATUS)
       endif
   
       call MAPL_FieldBundleGetPointer(Pert,  "QI", QIP, rc=STATUS)
       VERIFY_(STATUS)
       if (DO_MOIST_PHYS .eq. 0) then
          call MAPL_FieldBundleGetPointer(Traj3,  "QI", QIT, rc=STATUS)
          VERIFY_(STATUS)
       endif
   
       if (DO_MOIST_PHYS .ne. 0) then 
          !Cloud Trajectory Large scale and Convective
          call MAPL_FieldBundleGetPointer(Traj3,  "QLS", QLST, rc=STATUS)
          VERIFY_(STATUS)
          call MAPL_FieldBundleGetPointer(Traj3,  "QCN", QCNT, rc=STATUS)
          VERIFY_(STATUS)
       endif
   
       !Fill rank 2 pointers
       call MAPL_FieldBundleGetPointer(Traj2, "PS"     , PS     , rc=STATUS)
       VERIFY_(STATUS)
   
       call MAPL_FieldBundleGetPointer(Traj2, "FROCEAN", FROCEAN, rc=STATUS)
       VERIFY_(STATUS)
   
       call MAPL_FieldBundleGetPointer(Traj2, "FRLAND" , FRLAND , rc=STATUS)
       VERIFY_(STATUS)
   
       call MAPL_FieldBundleGetPointer(Traj2, "VARFLT" , VARFLT , rc=STATUS)
       VERIFY_(STATUS)
   
       call MAPL_FieldBundleGetPointer(Traj2, "USTAR"  , USTAR  , rc=STATUS)
       VERIFY_(STATUS)
   
       call MAPL_FieldBundleGetPointer(Traj2, "BSTAR"  , BSTAR  , rc=STATUS)
       VERIFY_(STATUS)
   
       call MAPL_FieldBundleGetPointer(Traj2, "ZPBL"   , ZPBL   , rc=STATUS)
       VERIFY_(STATUS)
   
       call MAPL_FieldBundleGetPointer(Traj2, "CM"     , CM     , rc=STATUS)
       VERIFY_(STATUS)
   
       call MAPL_FieldBundleGetPointer(Traj2, "CT"     , CT     , rc=STATUS)
       VERIFY_(STATUS)
   
       call MAPL_FieldBundleGetPointer(Traj2, "CQ"     , CQ     , rc=STATUS)
       VERIFY_(STATUS)
   
       !Grid dimensions for allocations
       im = size(UDP,1)
       jm = size(UDP,2)
       lm = size(UDP,3)
  
       !Winds
       allocate(UT(im,jm,lm)      , stat=status); VERIFY_(STATUS)
       allocate(VT(im,jm,lm)      , stat=status); VERIFY_(STATUS)
       allocate(UP(im,jm,lm)      , stat=status); VERIFY_(STATUS)
       allocate(VP(im,jm,lm)      , stat=status); VERIFY_(STATUS)
 
       !Cloud trajectory arrays
       allocate(Ph(im,jm,lm)      , stat=status); VERIFY_(STATUS)
       allocate(PIh(im,jm,lm)     , stat=status); VERIFY_(STATUS)
       allocate(TEMP(im,jm,lm)    , stat=status); VERIFY_(STATUS)
       allocate(fQi(im,jm,lm)     , stat=status); VERIFY_(STATUS)
   
       !Diagonals of arrays for solving diffusion.
       allocate(AKV(im,jm,lm)     , stat=status); VERIFY_(STATUS)
       allocate(BKV(im,jm,lm)     , stat=status); VERIFY_(STATUS)
       allocate(CKV(im,jm,lm)     , stat=status); VERIFY_(STATUS)
       allocate(TKV(im,jm,lm)     , stat=status); VERIFY_(STATUS)
       allocate(AKS(im,jm,lm)     , stat=status); VERIFY_(STATUS)
       allocate(BKS(im,jm,lm)     , stat=status); VERIFY_(STATUS)
       allocate(CKS(im,jm,lm)     , stat=status); VERIFY_(STATUS)
       allocate(TKS(im,jm,lm)     , stat=status); VERIFY_(STATUS)
       allocate(AKQ(im,jm,lm)     , stat=status); VERIFY_(STATUS)
       allocate(BKQ(im,jm,lm)     , stat=status); VERIFY_(STATUS)
       allocate(CKQ(im,jm,lm)     , stat=status); VERIFY_(STATUS)
       allocate(TKQ(im,jm,lm)     , stat=status); VERIFY_(STATUS)
   
       !Momentum mixing factor and topographic roughness factor
       allocate(EKV(im,jm,lm)     , stat=status); VERIFY_(STATUS)
       allocate(FKV(im,jm,lm)     , stat=status); VERIFY_(STATUS)
   
       !Pressure variables
       allocate(PrT(im,jm,0:lm)   , stat=status); VERIFY_(STATUS)    
       allocate(DPT_PS(im,jm,lm)  , stat=status); VERIFY_(STATUS)
       allocate(pkt(im,jm,  lm)   , stat=status); VERIFY_(STATUS)
       allocate(pet(im,jm,0:lm)   , stat=status); VERIFY_(STATUS)
       allocate(bx (im,jm,  lm)   , stat=status); VERIFY_(STATUS)
       allocate(cx (im,jm,  lm)   , stat=status); VERIFY_(STATUS)
       allocate(pke(im,jm,0:lm)   , stat=status); VERIFY_(STATUS)
   
       allocate(QIT1(im,jm,lm)    , stat=status); VERIFY_(STATUS)
       allocate(QLT1(im,jm,lm)    , stat=status); VERIFY_(STATUS)
   
       allocate(PTT1(im,jm,lm)    , stat=status); VERIFY_(STATUS)
       allocate(ZPBL1(im,jm)      , stat=status); VERIFY_(STATUS)
       allocate(CT1(im,jm)        , stat=status); VERIFY_(STATUS)
   
       !AK, BK, Pref
       allocate(AK(size(pert_ak)) , stat=status); VERIFY_(STATUS)
       allocate(BK(size(pert_bk)) , stat=status); VERIFY_(STATUS)
       allocate(Pref(0:size(AK)-1), stat=status); VERIFY_(STATUS)
   
       !Compute Pref from ak and bk
       AK = pert_ak
       BK = pert_bk
       DO i = 0,lm
          Pref(i) = AK(i+1) + BK(i+1)*MAPL_P00
       enddo
   
       !Transform winds to the A-grid
       !if (DO_PHYS_D2A == 0) then
          UT = UDT
          VT = VDT
          UP = UDP
          VP = VDP
       !else
          !call timing_on('PHYS_D2A')       
          !call ReStaggerWindsPhys(UDT, VDT, UT, VT, im, jm, lm, .true.)
          !call timing_off('PHYS_D2A')
          !UP = UDP !Dont restagger perturbation winds, will damp too much
          !VP = VDP
       !endif

       !Use local copied to avoid overwrite
       ZPBL1 = ZPBL
       CT1 = CT

       !Convert potential temperature to P0=10e5
       PTT1 = PTT*(MAPL_P00**MAPL_KAPPA)   

       !Turbulence Parameters
       TurbParams(1)  = 5.0          !LOUIS
       TurbParams(2)  = 160.0        !LAMBDAM
       TurbParams(3)  = 1.0          !LAMBDAM2
       TurbParams(4)  = 160.0        !LAMBDAH
       TurbParams(5)  = 1.0          !LAMBDAH2
       TurbParams(6)  = 3000.        !ZKMENV
       TurbParams(7)  = 3000.        !ZKHENV
       TurbParams(8)  = 0.1          !MINTHICK
       TurbParams(9)  = 0.0030       !MINSHEAR
       TurbParams(10) = 2.5101471e-8 !C_B
       TurbParams(11) = 1500.        !LAMBDA_B 
       TurbParams(12) = 500.         !AKHMMAX
       TurbParams(13) = 1.0          !PRANDTLSFC
       TurbParams(14) = 0.75         !PRANDTLRAD
       TurbParams(15) = 0.50         !BETA_RAD
       TurbParams(16) = 0.25         !BETA_SURF
       TurbParams(17) = 0.85         !KHRADFAC
       TurbParams(18) = 0.45         !KHSFCFAC
       TurbParams(19) = 20.0         !TPFAC_SURF
       TurbParams(20) = 1.5e-3       !ENTRATE_SURF
       TurbParams(21) = 0.5          !PCEFF_SURF
       TurbParams(22) = -999.        !LOUIS_MEMORY
       TurbParamsI(1) = count(PREF < 50000.)   !KPBLMIN
       TurbParamsI(2) = 1                      !LOCK_ON
       TurbParamsI(3) = 1                      !PBLHT_OPTION
       TurbParamsI(4) = 0                      !RADLW_DEP
   
       !Pressure at the half levels from Ps
       do i = 0,lm
          PrT(:,:,i) = AK(i+1) + BK(i+1)*PS(:,:)
       enddo
       !Pressure thickness
       DPT_PS(:,:,1:LM) = PrT(:,:,1:LM) - PrT(:,:,0:LM-1)
   
       call GetPressVarsTraj(DPT_PS, PET, PKE, PKT, BX, CX)
   
       !Calculate total cloud ice and liquid trajectory
       if (DO_MOIST_PHYS == 0) then
   
          QIT1 = QIT
          QLT1 = QLT
   
       else
   
          Ph = 0.01 * 0.5 * (PrT(:,:,0:LM-1) +  PrT(:,:,1:LM  ) ) 
          PIh = (Ph/1000.)**(MAPL_RGAS/MAPL_CP)
          TEMP = PTT1*PIh
          DO i = 1,IM
             DO j = 1,JM
                DO l = 1,LM
   
                   call IceFraction( TEMP(i,j,l), fQi(i,j,l) )
   
                enddo
             enddo
          enddo
   
          QIT1 = (QLST + QCNT) * fQi
          QLT1 = (QLST + QCNT) * (1 - fQi)
   
       endif
   
       !Initialize tri-diagonal arrays
       AKV = 0.0
       BKV = 0.0
       CKV = 0.0
       AKS = 0.0
       BKS = 0.0
       CKS = 0.0
       AKQ = 0.0
       BKQ = 0.0
       CKQ = 0.0
   
       !Call boundary layer routine. These routines will return the lower (AK*), main (BK*),
       !and upper (CK*) diagonals.
      
       call BL_DRIVER( IM           , &
                       JM           , &
                       LM           , &
                       DT           , &
                       UT          , &
                       VT          , &
                       PTT1         , &
                       QVT          , &
                       PET          , &
                       QIT1         , &
                       QLT1         , &
                       FRLAND       , &
                       FROCEAN      , &
                       VARFLT       , &
                       ZPBL1        , &
                       CM           , &
                       CT1          , &
                       CQ           , &
                       TURBPARAMS   , &
                       TURBPARAMSI  , &
                       USTAR        , &
                       BSTAR        , &
                       AKS, BKS, CKS, &
                       AKQ, BKQ, CKQ, &
                       AKV, BKV, CKV, &
                       EKV, FKV       &
                     )
       
       !Ks for the full BL are defined using P0=10^5
       if (PHASE==ADJphase) then          
          PTP = PTP/(MAPL_P00**MAPL_KAPPA)             
       elseif (PHASE==TLMphase) then
          PTP = PTP*(MAPL_P00**MAPL_KAPPA)
       endif
   
       !SOLVE DIFFUSION PROBLEM
   
       !Solver part 1, perform LU decomposition
       call VTRILUPERT(IM,JM,LM,AKV,BKV,CKV)
       call VTRILUPERT(IM,JM,LM,AKS,BKS,CKS)
       call VTRILUPERT(IM,JM,LM,AKQ,BKQ,CKQ)
      
       !Solver part 2, solve tri-diagonal in LU form
       !The last argument is either 1 or 0 and tells the solver whether YG would be present.
       !For winds, temperature and pressure it should be 1, for tracers its 0. Although
       !YGpert is 0 at the moment this can mess up the tracers if not taken into account.
       !Normally fluxes of tracers at the surface are not handled by the surface/turb component.
       call VTRISOLVEPERT(IM,JM,LM,AKV,BKV,CKV,UP ,PHASE,1)
       call VTRISOLVEPERT(IM,JM,LM,AKV,BKV,CKV,VP ,PHASE,1)
       call VTRISOLVEPERT(IM,JM,LM,AKS,BKS,CKS,PTP,PHASE,1)
       call VTRISOLVEPERT(IM,JM,LM,AKQ,BKQ,CKQ,QVP,PHASE,1)
       call VTRISOLVEPERT(IM,JM,LM,AKQ,BKQ,CKQ,QIP,PHASE,0)
       call VTRISOLVEPERT(IM,JM,LM,AKQ,BKQ,CKQ,QLP,PHASE,0)
       call VTRISOLVEPERT(IM,JM,LM,AKQ,BKQ,CKQ,O3P,PHASE,0)
       !call VTRISOLVEPERT(IM,JM,LM,AKQ,BKQ,CKQ,TrP,PHASE,0)
   
       !Convert potential temperature back to P0 = 1
       if (PHASE==ADJphase) then 
          PTP = PTP*(MAPL_P00**MAPL_KAPPA)
       elseif (PHASE==TLMphase) then
          PTP = PTP/(MAPL_P00**MAPL_KAPPA)
       endif
   
       UDP = UP
       VDP = VP

       deallocate(UT,VT,UP,VP)
       deallocate(Ph,PIh,TEMP,fQi)
       deallocate(AKV,BKV,CKV,TKV,AKS,BKS,CKS,TKS,AKQ,BKQ,CKQ,TKQ)
       deallocate(EKV,FKV)
       deallocate(PrT,DPT_PS,PKT,PET,BX,CX,PKE)
       deallocate(QIT1,QLT1)
       deallocate(AK,BK,PREF)
       deallocate(PTT1,ZPBL1,CT1)
      
       !Timers off
       call MAPL_TimerOff(MAPL,PHASE_NAME)
       call MAPL_TimerOff(MAPL,"TOTAL")
       if (PHASE==ADJphase) then  
          call timing_off('BDLAY_ADM')
       elseif (PHASE==TLMphase) then
          call timing_off('BDLAY_TLM')
       endif

    endif !DO_BL_PHYS

   ! All done
   ! --------

    RETURN_(ESMF_SUCCESS)

  end subroutine Run

  subroutine VTRILUPERT(IM,JM,LM,A,B,C)

   !Perform LU decomposition of tridiagonal system

   implicit none

   integer,                   intent(IN   ) :: IM, JM, LM
   real, dimension(IM,JM,LM), intent(IN   ) :: C
   real, dimension(IM,JM,LM), intent(INOUT) :: A, B

   integer :: L

    B(:,:,1) = 1. / B(:,:,1)
    do L = 2,LM
       A(:,:,L) = A(:,:,L) * B(:,:,L-1)
       B(:,:,L) = 1. / ( B(:,:,L) - C(:,:,L-1) * A(:,:,L) )
    end do

  end subroutine VTRILUPERT

  subroutine VTRISOLVEPERT(IM,JM,LM,A,B,C,Y,PHASE,YGSWITCH)

   !Solve LU decomposed tridiagonal system

   implicit none

   !Arguments
   integer,                   intent(IN   ) :: IM, JM, LM, PHASE, YGSWITCH
   real, dimension(IM,JM,LM), intent(IN   ) :: A, B, C
   real, dimension(IM,JM,LM), intent(INOUT) :: Y

   !Locals
   integer :: L

    if (PHASE == TLMPhase) then

       !Solve (LU)ynew = yold

       !Sweep down, modifying rhs with multiplier A
       do l=2,LM
          Y(:,:,l) = Y(:,:,l) - A(:,:,l)*Y(:,:,l-1)
       enddo

       !Sweep up, solving for updated value. Note B has the inverse of the main diagonal
       if (YGSWITCH == 1) then !Winds, temperature and q
          Y(:,:,LM) = Y(:,:,LM)*B(:,:,LM) !YGprime = 0
       else !Tracers
          Y(:,:,LM) = Y(:,:,LM)*B(:,:,LM-1)/(B(:,:,LM-1) - A(:,:,LM)*(1.0+C(:,:,LM-1)*B(:,:,LM-1) ))
       endif

       do l=LM-1,1,-1
          Y(:,:,l) = B(:,:,l)*(Y(:,:,l)-C(:,:,l)*Y(:,:,l+1) )
       enddo

    elseif (PHASE == ADJPhase) then

       if (YGSWITCH == 1) then !Can solve (LU)'ynew = U'L'ynew = yold

          !Sweep down but with U' instead of L
          Y(:,:,1) = Y(:,:,1)*B(:,:,1)
          do l=2,LM
             Y(:,:,l) = B(:,:,l) * (Y(:,:,l) - C(:,:,l-1)*Y(:,:,l-1))
          enddo

          !Sweep up but with L' instead of U
          do l=LM-1,1,-1
             Y(:,:,l) = Y(:,:,l) - A(:,:,l+1)*Y(:,:,l+1)
          enddo

       else !Change for surface means transpose doesnt work so use line-by-line adjoint

          !Adjoint of sweep up, solving for updated value. Note B has the inverse of the main diagonal
          DO l=1,lm-1,1
             Y(:,:,l+1) = Y(:,:,l+1) - C(:,:,l)*B(:,:,l)*Y(:,:,l)
             Y(:,:,l) = B(:,:,l)*Y(:,:,l)
          END DO

          !Adjoint of surface fix
          Y(:,:,lm) = B(:,:,lm-1)*Y(:,:,lm)/(B(:,:,lm-1)-A(:,:,lm)*(C(:,:,lm-1)*B(:,:,lm-1)+1.0))

          !Adjoint of sweep down, modifying rhs with multiplier A
          DO l=lm,2,-1
             Y(:,:,l-1) = Y(:,:,l-1) - A(:,:,l)*Y(:,:,l)
          END DO

       endif

    endif

  end subroutine VTRISOLVEPERT

  subroutine solve_tridiag(a,b,c,v,n)

   !d.holdaway - this solver does not do the correct thing as it doenst account for YG.
   !For tracers this is problematic since surface fluxes are not handled by
   !this part of the code, this solver assumes they are. This would also be problematic
   !if YGprime is computed in the future and not set to zero.

   implicit none

    !  a - sub-diagonal (means it is the diagonal below the main diagonal)
    !  b - the main diagonal
    !  c - sup-diagonal (means it is the diagonal above the main diagonal)
    !  v - right part
    !  x - the answer
    !  n - number of equations
 
   integer :: n
   real,dimension(n) :: a,b,c,v
   real,dimension(n) :: x
   real,dimension(n) :: bp,vp
   real :: m
   integer i

   !Make copies of the b and v variables so that they are unaltered by this sub
   bp(1) = b(1)
   vp(1) = v(1)
 
   !The first pass (setting coefficients):
   firstpass: do i = 2,n
      m = a(i)/bp(i-1)
      bp(i) = b(i) - m*c(i-1)
      vp(i) = v(i) - m*vp(i-1)
   end do firstpass
 
   x(n) = vp(n)/bp(n)

   !The second pass (back-substition)
   backsub:do i = n-1, 1, -1
      x(i) = (vp(i) - c(i)*x(i+1))/bp(i)
   end do backsub
 
   v(:) = x(:)

 end subroutine solve_tridiag

end module GEOS_TurbulencePertGridCompMod

