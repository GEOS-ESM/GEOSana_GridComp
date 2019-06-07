#include "MAPL_Generic.h"

#ifndef _CONNECT_ANA_TO_GCS_

module GEOS_AnaGridCompMod
   use ESMF
   use MAPL_Mod,  only:  ProxySetServices => MAPL_GenericSetServices
   private
   public SetServices

contains
   subroutine SetServices ( GC, RC )
      type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
      integer, optional                  :: RC  ! return code
      call ProxySetServices ( GC, RC=RC )
   end subroutine SetServices

end module GEOS_AnaGridCompMod

#else /* _CONNECT_ANA_TO_GCS_ */

!=============================================================================
!BOP

! !MODULE: GEOS\_AnaGridCompMod -- Parent Aerosol/Chemistry Component

! !INTERFACE:

module GEOS_AnaGridCompMod

! !USES:

  use ESMF
  use MAPL_Mod

  use GEOSAana_GridCompMod, only : AanaSetServices => SetServices
! use GEOSOana_GridCompMod, only : OanaSetServices => SetServices

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION: This gridded component (GC) combines 
 
!EOP

  integer :: AANA ! Atmospheric Analysis
  integer :: OANA ! Oceanic Analysis

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION:  The SetServices for the Chemistry GC needs to register its
!   Initialize and Run.  It uses the MAPL_Generic construct for defining 
!   state specs and couplings among its children.  In addition, it creates the   
!   children GCs and runs their respective SetServices.
!
! !REVISION HISTORY:
!
!   13Aug2007 Todling/daSilva Add this component; still doesn't work as meant
!   08Jul2008 Todling         Merge fdda-b1p3 with das-215/MAPL update
!   04Aug2008 Todling         Somehow dependencies in compilation don't work w/ fake ocean

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Locals

    integer                    :: I

    integer                    :: n, id

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------
    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // '::' // Iam

! No services to register for this component
! ------------------------------------------

! Create childrens gridded components and invoke their SetServices
! ----------------------------------------------------------------
    AANA = MAPL_AddChild(GC, NAME='AANA', SS=AanaSetServices, RC=STATUS)
    VERIFY_(STATUS)
!   OANA = MAPL_AddChild(GC, NAME='OANA', SS=OanaSetServices, RC=STATUS)
!   VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)

! Create children's gridded components and invoke their SetServices
! -----------------------------------------------------------------
   call MAPL_GenericSetServices  ( GC, RC=STATUS )
   VERIFY_(STATUS)

   RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

!--------------------------------------------------------------
  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Initialize method of the Composite Analysis Gridded Component.

!EOP

! ErrLog Variables

    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: IAm

! Local vars
    type(ESMF_Grid)       :: grid
    type(ESMF_FieldBundle):: bundle

!=============================================================================

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

! Create grid for this component
!-------------------------------
    call MAPL_GridCreate(GC, rc=status)
    VERIFY_(STATUS)

! Call Initialize for every Child
!--------------------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize
end module GEOS_AnaGridCompMod

#endif /* _CONNECT_ANA_TO_GCS_ */
