#include "MAPL_Exceptions.h"

!=============================================================================
!BOP

! !MODULE: GEOS\_ChemGridCompMod -- Parent Aerosol/Chemistry Component

! !INTERFACE:

module GEOS_ChemPertGridCompMod

! !USES:

  use ESMF
  use MAPL_Mod
!  use Chem_Mod

!  use GEOS_ChemEnvPertGridCompMod, only : EChemSetServices  => SetServices
!  use GOCART_PertGridCompMod,      only : AChemSetServices  => SetServices
!  use StratChem_PertGridCompMod,   only : SChemSetServices  => SetServices
!  use GMIchem_PertGridCompMod,     only : GChemSetServices  => SetServices
  use GEOS_PChemPertGridCompMod,   only : PChemSetServices  => SetServices

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION: This gridded component (GC) combines 
 
!EOP

! IMPORTANT: If adding a new component, make sure to update private function GetProvider_()
! -----------------------------------------------------------------------------------------
  integer ::      CHEMENV
  integer ::        PCHEM
  integer ::       GOCART
  integer ::    STRATCHEM
  integer ::      GMICHEM

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer,             intent(  OUT) :: RC  ! return code

! !DESCRIPTION:  The SetServices for the Chemistry GC needs to register its
!   Initialize and Run.  It uses the MAPL_Generic construct for defining 
!   state specs and couplings among its children.  In addition, it creates the   
!   children GCs and runs their respective SetServices.

!EOP

!=============================================================================
!
! ErrLog Variables

    __Iam__('SetServices')
    character(len=ESMF_MAXSTR) :: COMP_NAME

! Locals

   type (ESMF_GridComp), pointer :: GCS(:)

    integer                    :: I, RATS_PROVIDER, AERO_PROVIDER
    type (ESMF_Config)         :: CF

    integer                    :: n, id
    
    INTEGER :: numRATs
    CHARACTER(LEN=ESMF_MAXSTR) :: specieName(8)
    CHARACTER(LEN=ESMF_MAXSTR) :: providerName
    CHARACTER(LEN=ESMF_MAXSTR) :: shortName

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )
    Iam = trim(COMP_NAME) // '::' // trim(Iam)

! Register services for this component
! ------------------------------------
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Init, __RC__ )
!    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run,  __RC__ )

! Create children's gridded components and invoke their SetServices
! -----------------------------------------------------------------
        PCHEM = -1
       GOCART = -1
    STRATCHEM = -1
      GMICHEM = -1

      PCHEM = MAPL_AddChild(GC, NAME='PCHEM_PERT',  SS=PChemSetServices, __RC__)

! A container for the friendly tracers
! ------------------------------------
    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'CHEM_TRACERS',                              &
         LONG_NAME  = 'chemistry_friendly_tracers',                &
         UNITS      = 'X',                                         &
         DATATYPE   = MAPL_BundleItem,                             &
         __RC__  )

! Radiatively Active Tracers (RATs)
   specieName(1) = "OX"
   specieName(2) = "O3"
   specieName(3) = "CH4"
   specieName(4) = "N2O"
   specieName(5) = "CFC11"
   specieName(6) = "CFC12"
   specieName(7) = "CFC22"
   specieName(8) = "O3PPMV"
   numRATs = 8
   RATS_PROVIDER = PCHEM

  IF(MAPL_AM_I_ROOT()) THEN
   PRINT *," "
   PRINT *,TRIM(Iam)//": Radiatively active tracers (RATs):"
  END IF

! Add export specs for RATs.
! NOTE: See O3PPMV below.
! --------------------------
  DO i = 1, numRATs
  
   IF(RATS_PROVIDER == PCHEM) THEN
    shortName = TRIM(specieName(i))
   ELSE
    shortName = TRIM(providerName)//"::"//TRIM(specieName(i))
   END IF

   CALL MAPL_AddExportSpec( GC, SHORT_NAME = TRIM(shortName), &
  			    CHILD_ID = RATS_PROVIDER, __RC__ )
   IF(MAPL_AM_I_ROOT()) PRINT *," "//TRIM(shortName)
  END DO

! O3PPMV is treated as a special case, since it is not
! a member of the species bundle in STRATCHEM or GMICHEM
! ------------------------------------------------------
  IF(RATS_PROVIDER /= PCHEM) THEN
   CALL MAPL_AddExportSpec( GC, SHORT_NAME = "O3PPMV", &
  			    CHILD_ID = RATS_PROVIDER, __RC__ )
   IF(MAPL_AM_I_ROOT()) PRINT *," "//TRIM(providerName)//"::O3PPMV"
  END IF

! Aerosol for radiation
! ---------------------
  call ESMF_ConfigGetAttribute(CF, providerName, Default="PCHEM", &
                               Label="AERO_PROVIDER:", __RC__ )

  SELECT CASE (TRIM(providerName))
  CASE ("GOCART")
   AERO_PROVIDER = GOCART
  CASE DEFAULT
   AERO_PROVIDER = PCHEM
  END SELECT
 
  call MAPL_AddExportSpec ( GC, SHORT_NAME = 'AERO', &
                            CHILD_ID = AERO_PROVIDER, __RC__  )

  IF(MAPL_AM_I_ROOT()) THEN
   PRINT *," "
   PRINT *,"AERO Provider: ", TRIM(providerName)
  END IF

! H2O/O3 Tendencies
! -----------------
  call MAPL_AddExportSpec ( GC, SHORT_NAME = 'H2O_TEND', &
                            CHILD_ID = RATS_PROVIDER, __RC__  )

  call MAPL_AddExportSpec ( GC, SHORT_NAME = 'OX_TEND', &
                            CHILD_ID = RATS_PROVIDER, __RC__  )


! Finally, set the services
! -------------------------
  call MAPL_GenericSetServices ( GC, __RC__ )

  RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: Initialize -- Initialized method for composite Aero-Chemistry

! !INTERFACE:

  subroutine Init ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Initialize method of the Chemistry Composite Gridded 
!  Component. It acts as a driver for the initializtion of the children.

!EOP

! ErrLog Variables

  __Iam__('Init')
  character(len=ESMF_MAXSTR)           :: COMP_NAME

! Local derived type aliases

   type (MAPL_MetaComp),       pointer :: MAPL
   type (ESMF_GridComp),       pointer :: GCS(:)
   type (ESMF_State)                   :: INTERNAL
   type (ESMF_State),          pointer :: GEX(:)
   type (ESMF_FieldBundle)                  :: fBUNDLE, xBUNDLE, cBUNDLE
   type (ESMF_Field)                   :: FIELD
   type (ESMF_Config)                  :: CF

   integer                             :: AERO_PROVIDER
   integer                             :: I, xNA, cNA

!=============================================================================
 
! Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet ( GC, name=COMP_NAME, CONFIG=CF, __RC__ )
    Iam = trim(COMP_NAME) // "::" // trim(Iam)

!   Call GenericInitialize for every Child
!   --------------------------------------
    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK, __RC__ )

!   Get my MAPL_Generic state
!   --------------------------
    call MAPL_GetObjectFromGC ( GC, MAPL, __RC__ )

!   Get children and their im/ex states from my generic state.
!   ----------------------------------------------------------
    call MAPL_Get ( MAPL, GCS=GCS, GEX=GEX,              &
                    INTERNAL_ESMF_STATE=INTERNAL, __RC__ )

!   Fill in INTERNAL friendly bundle: CHEM_TRACERS
!   VERY IMPORTANT: Only the RATS provider can make OX friendly
!                   to ANALYSIS.
!   -----------------------------------------------------------
    call ESMF_StateGet (INTERNAL, 'CHEM_TRACERS', fBUNDLE, __RC__ )
    call MAPL_GridCompGetFriendlies(GCS,                     &
                                         (/ "ANALYSIS  ",    &
                                            "TURBULENCE",    &
                                            "DYNAMICS  ",    &
                                            "MOIST     " /), &
                                           fBUNDLE, __RC__ )

!   In early versions of MAPL the AERO bundle was not filled in
!   automatically with the contents of he AERO bundle from the
!   AERO_PROVIDER child. Here it checks whether this is the case, 
!   filling in the AERO bundle if necessary
!   -----------------------------------------------------------
    AERO_PROVIDER = GetProvider_(CF,'AERO_PROVIDER:',              __RC__ )
    call ESMF_StateGet(EXPORT, 'AERO',             xBUNDLE, __RC__ )
    call ESMF_StateGet(GEX(AERO_PROVIDER), 'AERO', cBUNDLE, __RC__ )

!   Let's see how many we've got
!   ----------------------------
    call ESMF_FieldBundleGet(xBUNDLE, fieldCount=xNA, __RC__ )
    call ESMF_FieldBundleGet(cBUNDLE, fieldCount=cNA, __RC__ )
    
!   If the AERO bundle is empty and the provider has an AERO bundle
!   which is not empty, fill the one in our EXPORT with its contents
!   ----------------------------------------------------------------
    if ( xNA==0 .AND. cNA > xNA ) then 
       DO i = 1, cNA
          call ESMF_FieldBundleGet(cBUNDLE, i, FIELD, __RC__ ) ! from child
          call MAPL_FieldBundleAdd(xBUNDLE,    FIELD, __RC__ ) ! to Chem export
       end do
    else
       if (xNA/=cNA) then
          IF(MAPL_AM_I_ROOT()) PRINT *, TRIM(Iam)//": Inconsistent number of fields in AERO bundles from export and child states."
       end if
    end if


#ifdef PRINT_STATES

!   Print what my states are
!   ------------------------
    if ( MAPL_am_I_root() ) then

       print *,  trim(Iam)//": IMPORT State" 
                                             call ESMF_StatePrint ( IMPORT)
       print *,  trim(Iam)//": INTERNAL State" 
                                             call ESMF_StatePrint ( INTERNAL )
       print *,  trim(Iam)//": EXPORT State" 
                                             call ESMF_StatePrint ( EXPORT )

       print *,  trim(Iam)//": AERO Bundle (EXPORT)"
                                             call ESMF_FieldBundlePrint ( xBUNDLE )

       print *,  trim(Iam)//": AERO Bundle (PROVIDER)",  AERO_PROVIDER  
                                             call ESMF_FieldBundlePrint ( cBUNDLE )

       print *,  trim(Iam)//": Friendly Tracers (INTERNAL)" 
                                             call ESMF_FieldBundlePrint ( fBUNDLE )
    end if

#endif

!   All Done
!   --------
    RETURN_(ESMF_SUCCESS)

  end subroutine Init

!-----------------------------------------------------------------------

     integer function GetProvider_ ( CF, Label, RC )
!
!    Returns provider name as per resource file.
!
     type (ESMF_Config), intent(inout) :: CF
     character(len=*),   intent(in)    :: Label
     integer, intent(out)              :: RC

! ErrLog Variables

     __Iam__('GetProvider_')

     character(len=ESMF_MAXSTR)     :: provider


     call ESMF_ConfigGetAttribute(CF, provider, Default='PCHEM', &
                                  Label=Label, __RC__ )

     select case ( trim(provider) )

           case ('PCHEM')
                                    GetProvider_ = PCHEM
           case ('GOCART')
                                    GetProvider_ = GOCART
           case ('STRATCHEM')
                                    GetProvider_ = STRATCHEM
           case ('GMICHEM')
                                    GetProvider_ = GMICHEM
           case DEFAULT

                __raise__(MAPL_RC_ERROR,'Unknown provider '//TRIM(provider))

     end select
     
     RC = ESMF_SUCCESS

   end function GetProvider_

 end module GEOS_ChemPertGridCompMod


