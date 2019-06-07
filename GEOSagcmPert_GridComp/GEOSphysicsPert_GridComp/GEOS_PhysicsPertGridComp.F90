! $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: 

! GEOS\_PhysicsGridCompMod -- A Module to combine Short-Wave, Long-Wave Radiation,
!                             Moist-Physics and Turbulence Gridded Components

! !INTERFACE:

module GEOS_PhysicsPertGridCompMod

! !USES:

  use ESMF
  use MAPL_Mod

  use GEOS_SurfacePertGridCompMod,    only : SurfSetServices      => SetServices
  use GEOS_MoistPertGridCompMod,      only : MoistSetServices     => SetServices
  use GEOS_TurbulencePertGridCompMod, only : TurblSetServices     => SetServices
  use GEOS_RadiationPertGridCompMod,  only : RadiationSetServices => SetServices
  use GEOS_ChemPertGridCompMod,       only : AChemSetServices     => SetServices
  use GEOS_GwdPertGridCompMod,        only : GwdSetServices       => SetServices
  use GEOS_PertSharedMod,             only : TLMPhase, ADJPhase
!  integer, parameter :: TLMPhase    = 1
!  integer, parameter :: ADJPhase    = 2

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION: This gridded component (GC) combines the Radiation (Short-Wave and Long-Wave), 
!   Moist-Physics, Chem, Surface and Turbulence GCs into a new composite Physics GC.

!EOP

  integer ::        GWD
  integer ::        SURF
  integer ::        CHEM
  integer ::        MOIST
  integer ::        TURBL
  integer ::        RAD

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer,             intent(  OUT) :: RC  ! return code

! !DESCRIPTION:  The SetServices for the Physics GC needs to register its
!   Initialize and Run.  It uses the MAPL_Generic construct for defining 
!   state specs and couplings among its children.  In addition, it creates the   
!   children GCs (SURF, CHEM, RADIATION, MOIST, TURBULENCE) and runs their
!   respective SetServices.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Locals

    CHARACTER(LEN=ESMF_MAXSTR)              :: RATsProviderName
    integer                                 :: I
    type (ESMF_Config)                      :: CF

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "::" // Iam

! Register services for this component
! ------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run,        RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run,        RC=STATUS )
    VERIFY_(STATUS)

! Create children`s gridded components and invoke their SetServices
! -----------------------------------------------------------------

    GWD   = MAPL_AddChild(GC, NAME='GWD_PERT',        SS=GwdSetServices,       RC=STATUS)
    VERIFY_(STATUS)
    MOIST = MAPL_AddChild(GC, NAME='MOIST_PERT',      SS=MoistSetServices,     RC=STATUS)
    VERIFY_(STATUS)
    SURF  = MAPL_AddChild(GC, NAME='SURFACE_PERT',    SS=SurfSetServices,      RC=STATUS)
    VERIFY_(STATUS)
    TURBL = MAPL_AddChild(GC, NAME='TURBULENCE_PERT', SS=TurblSetServices,     RC=STATUS)
    VERIFY_(STATUS)
    CHEM  = MAPL_AddChild(GC, NAME='CHEMISTRY_PERT',  SS=AChemSetServices,     RC=STATUS)
    VERIFY_(STATUS)
    RAD   = MAPL_AddChild(GC, NAME='RADIATION_PERT',  SS=RadiationSetServices, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerAdd(GC, name="RUNTLM"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="RUNADJ"    ,RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GenericSetServices ( GC, RC=STATUS )
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

  type (ESMF_State),        pointer  :: GIM(:)
  type (ESMF_FieldBundle)            :: Traj3, Traj2, Pert
  integer                            :: nc
 
! Begin
! -----

    Iam = "Initialize"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Retrieve the pointer to the MAPL object
! ---------------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"INITIALIZE")
    call MAPL_TimerOn(MAPL,"TOTAL")

! Get my children GCs and IMPORT States
!--------------------------------------

    call MAPL_Get ( MAPL, GIM=GIM, RC=STATUS )
    VERIFY_(STATUS)

! Get working trajectory and pertubation from IMPORT  
!---------------------------------------------------

    call ESMF_StateGet(Import, 'TRAJWRK3', Traj3, rc=STATUS)
    VERIFY_(STATUS)
    call ESMF_StateGet(Import, 'TRAJWRK2', Traj2, rc=STATUS)
    VERIFY_(STATUS)
    call ESMF_StateGet(Import, 'PERTWRK', Pert, rc=STATUS)
    VERIFY_(STATUS)

! Put the trajectory and perturbation in the childrens imports
!-------------------------------------------------------------

    do nc=1,size(GIM)
       call MAPL_StateAdd(GIM(nc), Traj3, rc=status)
       VERIFY_(STATUS)
       call MAPL_StateAdd(GIM(nc), Traj2, rc=status)
       VERIFY_(STATUS)
       call MAPL_StateAdd(GIM(nc), Pert, rc=status)
       VERIFY_(STATUS)
    end do

    call MAPL_TimerOff(MAPL,"TOTAL")

! Call Generic Initialize
! -----------------------
    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

! All done
! --------

    call MAPL_TimerOff(MAPL,"INITIALIZE")

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize

!---------------------------------------------------------------------

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
  integer,  optional,  intent(  out) :: rc 

!EOP

! !Local Variables:
  
    type (MAPL_MetaComp),      pointer :: MAPL 
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: IAm
    character(len=ESMF_MAXSTR)         :: COMP_NAME

    character(len=ESMF_MAXSTR)         :: PHASE_NAME
    integer                            :: nn, nc, phase, idmodel
    logical,save                       :: first=.true.
    type (ESMF_GridComp),      pointer :: GCS(:)
    type (ESMF_State),         pointer :: GIM(:)
    type (ESMF_State),         pointer :: GEX(:)

! Begin
!------

    Iam = "Run"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

! Retrieve the pointer to the generic state
! -----------------------------------------

    call MAPL_GetObjectFromGC (GC, MAPL,  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, idmodel, 'AGCM_PERT_IDMODEL:', default=0, RC=STATUS )
    VERIFY_(STATUS)
    if (idmodel/=0) then
       if (first) call write_parallel("Skip physics")
       first = .false.
       RETURN_(ESMF_SUCCESS)
    endif

    call MAPL_TimerOn(MAPL,"TOTAL")

    call ESMF_GridCompGet(GC, currentPhase=Phase, rc=STATUS)
    VERIFY_(STATUS)

    select case(phase)
    case(TLMphase)
       PHASE_NAME = "RUNTLM"
    case(ADJphase)
       PHASE_NAME = "RUNADJ"
    case default
       ASSERT_(.false.)
    end select

    call MAPL_TimerOn(MAPL,PHASE_NAME)

! Get my children GCs and IM/EX States
!--------------------------------------

    call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX, RC=STATUS )
    VERIFY_(STATUS)

! Call children
! -------------

    do nc=1,size(GCS)
       if(phase==ADJphase) then
          nn = size(GCS) + 1 - nc
       else
          nn = nc
       end if

       call ESMF_GridCompRun (GCS(nn), importState=GIM(nn), exportState=GEX(nn), &
            clock=CLOCK, PHASE=PHASE, userRC=STATUS )
       VERIFY_(STATUS)
    end do

! All done
! --------
    call MAPL_TimerOff(MAPL,PHASE_NAME)
    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)

  end subroutine Run

end module GEOS_PhysicsPertGridCompMod
