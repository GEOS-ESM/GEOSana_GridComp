
!  $Id$
#include "MAPL_Generic.h"
#define MAPL_FieldBundleGetPointer ESMFL_BundleGetPointerToData

!=============================================================================
!BOP

! !MODULE: GEOS_Gwd -- A Module to compute the forcing due to parameterized
!          gravity wave drag

! !INTERFACE:

module GEOS_GwdPertGridCompMod

! !USES:

  use ESMF
  use MAPL_Mod
  use GEOS_Mod
  use GEOS_PertSharedMod

  use GW_DRAG_D
  use GW_DRAG_B

  use fv_timing_mod,  only: timing_on, timing_off

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION:
! 
!   {\tt GWD} is a light-weight gridded component to compute the forcing
! due to gravity wave drags. It operates on the ESMF grid that appears in the
! gridded component passed to its {\tt Initialize} method. Unlike
! heavier gridded components, it does not enforce its own grid.
! The only restrictions are that it be a 3-dimensional grid
! in which one dimension is aligned with the vertical coordinate and
! only the horizontal dimensions are decomposed.
!
! The gravity wave drag scheme is based on NCAR WACCM1b gw_drag routine.
! The scheme includes parameterizations for orographic (stationary) gravity
! waves (Kiehl et al. 1996), and for a spectrum of traveling gravity waves 
!(Sassi et al. 2003; http://acd.ucar.edu/models/WACCM). Both parameteriz-
! ations are based on Lindzen's [1981] formulation. The interested reader 
! is referred to those publications for details of the mathematical
! derivations.
!

!EOP

contains

!BOP

! ! IROUTINE: SetServices -- Sets ESMF services for this component

! ! INTERFACE:

  subroutine SetServices ( GC, RC )

! ! ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! ! DESCRIPTION: This version uses the MAPL_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF_State INTERNAL, which is in the MAPL_MetaComp.

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

! Set the Run entry point
! -----------------------


! Set the state variable specs.
! -----------------------------
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

  subroutine Initialize ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: gc       ! composite gridded component 
  type(ESMF_State),    intent(inout) :: import   ! import state
  type(ESMF_State),    intent(inout) :: export   ! export state
  type(ESMF_Clock),    intent(inout) :: clock    ! the clock
  integer, optional,   intent(  out) :: rc       ! Error code:
                                                 ! = 0 all is well
                                                 ! otherwise, error

! Locals
!-------

  type (MAPL_MetaComp),      pointer :: MAPL 
  
  integer                            :: STATUS
  character(len=ESMF_MAXSTR)         :: IAm
  character(len=ESMF_MAXSTR)         :: COMP_NAME

  type (ESMF_Field)                  :: Field
  type (ESMF_Grid)                   :: ESMFGRID

  integer                            :: Num, I, CommID

  type (ESMF_State)                  :: INTERNAL


! Begin
! -----

    Iam = "Initialize"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize

  subroutine Run(gc, import, export, clock, rc)

! !ARGUMENTS:
  type(ESMF_GridComp), intent(inout) :: gc
  type (ESMF_State),   intent(inout) :: import
  type (ESMF_State),   intent(inout) :: export
  type (ESMF_Clock),   intent(inout) :: clock
  integer,  optional,  intent(  out) :: rc 

! !Local Variables:  
  type (MAPL_MetaComp),      pointer :: MAPL 
  type (ESMF_Alarm)                  :: ALARM

  integer                            :: STATUS
  character(len=ESMF_MAXSTR)         :: IAm
  character(len=ESMF_MAXSTR)         :: COMP_NAME
  character(len=ESMF_MAXSTR)         :: PHASE_NAME
  integer                            :: PHASE
  integer                            :: IM, JM, LM
  integer                            :: I, J, K, L

!EOP
  type (ESMF_FieldBundle)            :: Traj3, Traj2, Pert
  type (ESMF_State)                  :: INTERNAL

  real, pointer                      :: UT (:,:,:), UP (:,:,:)
  real, pointer                      :: VT (:,:,:), VP (:,:,:)
  real, pointer                      :: PTT(:,:,:)
  real, pointer                      :: QVT(:,:,:)
  real, pointer                      :: SGH(:,:)
  real, pointer                      :: Ps(:,:)
  real, pointer                      :: LATS(:,:) 
  real, pointer                      :: LONS(:,:) 
  real                               :: DT

  real, allocatable                  :: UT8 (:,:,:)
  real, allocatable                  :: VT8 (:,:,:)

  real, allocatable                  :: ak(:), bk(:)
  real, allocatable                  :: PLE(:,:,:) 
  real, allocatable                  :: PTT1(:,:,:) 
  real, allocatable                  :: PMID(:,:,:) ! pressure at the layers
  real, allocatable                  :: PDEL(:,:,:) ! pressure thickness at the layers
  real, allocatable                  :: RPDEL(:,:,:)
  real, allocatable                  :: PILN(:,:,:)
  real, allocatable                  :: PMLN(:,:,:)
  real, allocatable                  :: EFFGWO(:,:), EFFGWB(:,:)
  real, allocatable                  :: ZI(:,:,:), ZM(:,:,:)
  real, allocatable                  :: pref(:)
  real, allocatable                  :: ATT(:,:,:)
  real, allocatable                  :: UP_gwd(:,:,:)
  real, allocatable                  :: VP_gwd(:,:,:)

  integer, parameter                 :: pgwv    = 4
  real, parameter                    :: effgworo= 0.125, effgwbkg= 0.25
  real, parameter                    :: efforonh= 2.5  , efforosh= 1.0

  real, parameter                    :: PI_GWD  = 4.0*atan(1.0)  
  real, parameter                    :: ROG     = MAPL_RGAS/MAPL_GRAV
  real, parameter                    :: P00     = MAPL_P00

  real                               :: HVSD

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

    call MAPL_Get ( MAPL, lats = LATS, lons = LONS,  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL, DT, Label="RUN_DT:", RC=STATUS)
    VERIFY_(STATUS)

    select case(phase)
    case(TLMphase)
       PHASE_NAME = "RUNTLM"
    case(ADJphase)
       PHASE_NAME = "RUNADJ"
    case default
       ASSERT_(.false.)
    end select

    IF (DO_GWD_PHYS .ne. 0) THEN

    if (PHASE==ADJphase) then         
       call timing_on('GWAVD_ADM')             
    elseif (PHASE==TLMphase) then         
       call timing_on('GWAVD_TLM')            
    endif

    call ESMF_StateGet(import, "TRAJWRK3", Traj3, rc=status)
    VERIFY_(STATUS)
    call ESMF_StateGet(Import, 'TRAJWRK2', Traj2, rc=STATUS)
    VERIFY_(STATUS)
    call ESMF_StateGet(Import, "PERTWRK", Pert, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_FieldBundleGetPointer(Traj3, "U", UT, rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_FieldBundleGetPointer(Pert,  "U", UP, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_FieldBundleGetPointer(Traj3, "V", VT, rc=STATUS)
    VERIFY_(STATUS)
    call MAPL_FieldBundleGetPointer(Pert,  "V", VP, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_FieldBundleGetPointer(Traj3, "PT", PTT, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_FieldBundleGetPointer(Traj3, "QV", QVT, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_FieldBundleGetPointer(Traj2, "HS_STDV", SGH, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_FieldBundleGetPointer(Traj2, "PS"     , PS,  rc=STATUS)
    VERIFY_(STATUS)

    IM = size(PTT,1)
    JM = size(PTT,2)
    LM = size(PTT,3)

    allocate(PLE(IM,JM,LM+1))
    allocate(PTT1(IM,JM,LM))
    allocate(PMID(IM,JM,LM))
    allocate(PDEL(IM,JM,LM))
    allocate(RPDEL(IM,JM,LM))
    allocate(PILN(IM,JM,LM+1))
    allocate(PMLN(IM,JM,LM))
    allocate(EFFGWO(IM,JM), EFFGWB(IM,JM))
    allocate(ZI(IM,JM,LM+1))
    allocate(ZM(IM,JM,LM))
    allocate(ak(LM+1))
    allocate(bk(LM+1))
    allocate(PREF(LM+1))
    allocate(ATT(IM,JM,LM))
    allocate(UP_gwd(IM,JM,LM))
    allocate(VP_gwd(IM,JM,LM))

    allocate(  UT8(IM,JM,LM))
    allocate(  VT8(IM,JM,LM))

    DO L = 1, LM
    DO J=1,JM
      DO I=1,IM
         UT8(I,J,L)  = UT(I,J,L)
         VT8(I,J,L)  = VT(I,J,L)
       ENDDO
    ENDDO
    ENDDO

    ak = pert_ak
    bk = pert_bk

    PREF = AK + BK * P00

    DO L = 1,LM+1
       PLE(:,:,L) = AK(L) + BK(L)*PS(:,:)
    ENDDO

    PTT1 = PTT

    IF (PHASE_NAME .eq. "RUNTLM") then
       DO L = 1, LM
          DO J=1,JM
             DO I=1,IM
                UP_gwd(I,J,L)  = UP(I,J,L)
                VP_gwd(I,J,L)  = VP(I,J,L)
             ENDDO
          ENDDO
       END DO
    ELSE ! ADJ PHASE
       DO L = 1, LM
          DO J=1,JM
             DO I=1,IM
                UP_gwd(I,J,L)  = UP(I,J,L)
                VP_gwd(I,J,L)  = VP(I,J,L)
             ENDDO
          ENDDO
       END DO
    ENDIF

    DO J = 1, JM
    DO I = 1, IM
       DO K = 1, LM
           PMID(I,J,K) = 0.5*(  PLE(I,J,K  ) + PLE(I,J,K+1) )
           PDEL(I,J,K) =        PLE(I,J,K+1) - PLE(I,J,K  )
           RPDEL(I,J,K)= 1.0 / PDEL(I,J,K)
           PILN(I,J,K) = log(   PLE(I,J,K) )
           PMLN(I,J,K) = log(  PMID(I,J,K) ) !
       END DO
       PILN(I,J,LM+1)  = log( PLE(I,J,LM+1)  )
       IF (LATS(I,J)*180./PI_GWD < -20.) THEN
          HVSD = efforosh - 0.75*exp(-abs(LATS(I,J)*180./PI_GWD+20.)/5.)
       ELSE
          HVSD = efforonh + 0.75*exp(-abs(LATS(I,J)*180./PI_GWD+20.)/5.)
       END IF
       EFFGWO(I,J) = EFFGWORO*HVSD
       EFFGWB(I,J) = EFFGWBKG
    ENDDO
    ENDDO

    IF (PHASE_NAME .eq. "RUNTLM") then
    CALL GW_main_d(IM*JM,    LM,      DT,                       &
                   PGWV,     EFFGWO,  EFFGWB,                    &
                   PLE,      PTT1,    UT8,       up_gwd,         &
                   VT8,      vp_gwd,  SGH,       PREF,           &
                   PMID,     PDEL,    RPDEL,     PILN,  ZM,  QVT, rog, mapl_vireps, LATS)
    ELSE ! ADM PHASE 
    CALL GW_main_b(IM*JM,    LM,      DT,                       &
                   PGWV,     EFFGWO,  EFFGWB,                    &
                   PLE,      PTT1,    UT8,       up_gwd,         &
                   VT8,      vp_gwd,  SGH,       PREF,           &
                   PMID,     PDEL,    RPDEL,     PILN,  ZM,  QVT, rog, mapl_vireps, LATS)
    ENDIF

    IF (PHASE_NAME .eq. "RUNTLM") then
       DO L = 1, LM
          DO J=1,JM
             DO I=1,IM
                UP(I,J,L)  = UP_gwd(I,J,L) 
                VP(I,J,L)  = VP_gwd(I,J,L) 
             ENDDO
          ENDDO
       END DO
    ELSE ! ADM PHASE
       DO L = 1, LM
          DO J=1,JM
             DO I=1,IM
                UP(I,J,L) = UP_gwd(I,J,L) 
                VP(I,J,L) = VP_gwd(I,J,L) 
             ENDDO
          ENDDO
       END DO
    ENDIF

    deallocate(PLE)
    deallocate(PTT1)
    deallocate(PMID)
    deallocate(PDEL)
    deallocate(RPDEL)
    deallocate(PILN)
    deallocate(PMLN)
    deallocate(EFFGWO, EFFGWB)
    deallocate(ZI, ZM)
    deallocate(ak)
    deallocate(bk)
    deallocate(PREF)
    deallocate(ATT)
    deallocate(UP_gwd)
    deallocate(VP_gwd)

    deallocate(UT8)
    deallocate(VT8)

    if (PHASE==ADJphase) then         
       call timing_off('GWAVD_ADM')             
    elseif (PHASE==TLMphase) then         
       call timing_off('GWAVD_TLM')            
    endif


    ENDIF !DO_GWD_PHYS

    RETURN_(ESMF_SUCCESS)
  end subroutine Run
end module GEOS_GwdPertGridCompMod
