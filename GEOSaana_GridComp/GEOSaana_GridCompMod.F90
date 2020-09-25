#include "MAPL_Generic.h"
!#define PRINT_STATES

   module GEOSaana_GridCompMod
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GEOSaana_GridCompMod --- GEOS Atmospheric Analysis Component 
! !DESCRIPTION: This gridded component provides a wrapper around the
!               GSI Grid Component in order to make it fully compatible
!  with the GEOS-5/MAPL system.
!
! This component is necessary because GEOS-5 and GSI have different grids 
! and since MAPL does not yet provide a coupling mechanism between components
! with different grids the coupling must be performed "manually". Hence the 
! need for this component.
!
! !REVISION HISTORY:
!
! 01Feb2007  Cruz  First implementation.

! !USES:

   use ESMF
   use MAPL
   use GSI_GridCompMod, only: GSI_SetServices => GSI_GridCompSetServices
   use GSI_GridCompMod, only: AANA_SetupSpecs => GSI_GridCompSetupSpecs
   use GSI_GridCompMod, only: AANA_SetAlarms  => GSI_GridCompSetAlarms
   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:

   public SetServices

! !REVISION HISTORY:
!
!  08Jul2008  Todling  Merge fdda-b1p3 w/ das-215 (MAPL update)
!
!
!EOP
!-------------------------------------------------------------------------


! Internal state
! --------------
   type TGEOSaana
     private
     type(ESMF_GridComp) :: gcGSI      ! GSI gridded component
     type(ESMF_State)    :: impGSI     ! import for GSI
     type(ESMF_State)    :: expGSI     ! export for GSI
     type(ESMF_Clock)    :: clock          
   end type TGEOSaana

! Wrapper for extracting internal state
! -------------------------------------
   type Twrap
     type(TGEOSaana), pointer :: intSt
   end type Twrap

   integer, save        :: BKGfreq_sc
   logical, save        :: analysis_is_done=.false.

   CONTAINS

!BOPI

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

   subroutine SetServices ( gcAANA, RC )

! !ARGUMENTS:

   type(ESMF_GridComp)     :: gcAANA  ! this gridded component
   integer, intent(out)    :: RC      ! return code

! !DESCRIPTION: This function sets the IRF services,and creates the GSI 
!       sub-components (g.c., grid, and coupling states).
! 
! !REVISION HISTORY:
!
!  11Sep2007  Todling  Add logics to skip Init/Run/Final
!  11Dec2011  Todling  Do not finalize until destroy bug is resolved
!
!EOPI
!-------------------------------------------------------------------------

! Locals

   character(len=ESMF_MAXSTR) :: IAm
   integer                    :: STATUS
   character(len=ESMF_MAXSTR) :: COMP_NAME
   type(TGEOSaana), pointer   :: ThisIntSt
   type(Twrap)                :: wrap
   type(ESMF_Config)          :: cf 
   character(len=ESMF_MAXSTR) :: observer
   logical                    :: do_observer
   integer                    :: BKGfreq, BKGfreq_hr, BKGfreq_mn

! Begin...

   rc = 0
   call ESMF_GridCompGet( gcAANA, NAME=COMP_NAME, config=cf, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // 'SetServices'

   call ESMF_ConfigGetAttribute( cf, BKGfreq, label ='BKG_FREQUENCY:', rc = STATUS )
   VERIFY_(STATUS)

   BKGfreq_hr = BKGfreq/10000
   BKGfreq_mn = mod(BKGfreq,10000)/100
   BKGfreq_sc = mod(BKGfreq,100)
   BKGfreq_sc = BKGfreq_hr*3600 + BKGFreq_mn*60 + BKGfreq_sc

   do_observer = .false.
   if ( BKGfreq_sc>0 ) do_observer = .true.

   if ( do_observer ) then

!       Register services for this component
!       ------------------------------------

        call MAPL_GridCompSetEntryPoint ( gcAANA, ESMF_METHOD_INITIALIZE, Initialize,   &
                               STATUS )
        VERIFY_(STATUS)
        call MAPL_GridCompSetEntryPoint ( gcAANA, ESMF_METHOD_RUN,  Run       ,   &
                               STATUS )
        VERIFY_(STATUS)
!       call MAPL_GridCompSetEntryPoint ( gcAANA, ESMF_METHOD_FINALIZE,  Finalize,   &
!                              STATUS )
!       VERIFY_(STATUS)

!       Create the GSI sub-components
!       -----------------------------
        allocate(ThisIntSt, stat=status)
        VERIFY_(STATUS)
        call GEOSaana_Create_( gcAANA, ThisIntSt, rc=status)
        VERIFY_(STATUS)

!       Save the "internal" state
!       -------------------------
        wrap%intSt => ThisIntSt
        call ESMF_GridCompSetInternalState(gcAANA, wrap, status)
        VERIFY_(STATUS)
   
!       Register the GSI gridded component
!       ----------------------------------
        call ESMF_GridCompSetServices ( ThisIntSt%gcGSI, GSI_SetServices, RC=status )
        VERIFY_(STATUS)

!       Set Import/Export Coupling SPECS
!       --------------------------------
        call  AANA_SetupSpecs (gcAANA, opthw=0, rc=STATUS)
        VERIFY_(STATUS)

        if (MAPL_AM_I_ROOT()) &
        print *, trim(Iam)//': ACTIVE '

   else

        if (MAPL_AM_I_ROOT()) &
        print *, trim(Iam)//': NOT ACTIVE, defaulting to Generic No-op stubs'

   endif

!ALT: GSI and AANA imports are not compatible!
!     they have different halo width
! Although this is not an issue with the current version of GSI
! the next statement will not hurt, and if GSI ever becomes a MAPL child
! the next statement is REQUIRED!

   call MAPL_TerminateImport    ( gcAANA, ALL=.true., RC=STATUS  )
   VERIFY_(STATUS)

!  Generic SetServices
!  --------------------
   call MAPL_GenericSetServices (gcAANA, RC=STATUS)
   VERIFY_(STATUS)

   rc = STATUS
   RETURN_(ESMF_SUCCESS)
  
   end subroutine SetServices

!-------------------------------------------------------------------
!BOPI

! !IROUTINE: Initialize -- Initialize method

! !INTERFACE:

   subroutine Initialize ( gcAANA, impAANA, expAANA, CLOCK, RC )

! !ARGUMENTS:

   type(ESMF_GridComp), intent(inout) :: gcAANA  ! this gridded component 
   type(ESMF_State),    intent(inout) :: impAANA ! Import state
   type(ESMF_State),    intent(inout) :: expAANA ! Export state
   type(ESMF_Clock),    intent(inout) :: CLOCK   ! The clock
   integer, optional,   intent(  out) :: RC      ! Error code

! !DESCRIPTION: This function initializes the AANA gridded component.
!       It is a wrapper to the initialization of the GSI g.c.
!
! !REVISION HISTORY:
!
!   14Apr2007 Todling removed reference to numbkgs
! 
!EOPI
!-------------------------------------------------------------------------

   type (MAPL_MetaComp),  pointer  :: STATE
   character(len=ESMF_MAXSTR) :: IAm
   integer                    :: STATUS
   character(len=ESMF_MAXSTR) :: COMP_NAME
   type(TGEOSaana), pointer   :: ThisIntSt
   type(Twrap)                :: wrap
   type(ESMF_Config)          :: cf
   type(ESMF_Alarm)           :: Final_End_Time_Alarm

! Begin...

   call ESMF_GridCompGet ( gcAANA, name=COMP_NAME, config=cf, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // "Initialize"

! Retrieve internal state

   call ESMF_GridCompGetInternalState ( gcAANA, wrap, status )
   VERIFY_(STATUS)
   ThisIntSt => wrap%intST

! GEOS generic initialize
   
   call MAPL_GenericInitialize(gcAANA, impAANA, expAANA, clock, rc=STATUS)
   VERIFY_(STATUS)

!  Initialize GSI
!  --------------
   if ( MAPL_AM_I_ROOT() ) print *, 'Initialize GSI'      
   call ESMF_GridCompInitialize ( ThisIntSt%gcGSI, &
                                  importState=ThisIntSt%impGSI, &
                                  exportState=ThisIntSt%expGSI, &
                                  clock=clock, &
                                  rc=status )
   VERIFY_(STATUS) 

   call MAPL_GetObjectFromGC ( gcAANA, STATE, RC=STATUS)
   VERIFY_(STATUS)

   call AANA_SetAlarms(STATE, gcAANA, cf, clock)

#ifdef PRINT_STATES
    call WRITE_PARALLEL ( trim(Iam)//": IMPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( impAANA, rc=STATUS )
    call WRITE_PARALLEL ( trim(Iam)//": EXPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( expAANA, rc=STATUS )
#endif

   RETURN_(ESMF_SUCCESS)

   end subroutine Initialize

!-------------------------------------------------------------------
!BOPI

! !IROUTINE: Run -- Run method

! !INTERFACE:

   subroutine Run ( gcAANA, impAANA, expAANA, CLOCK, RC )

! !ARGUMENTS:

   type(ESMF_GridComp), intent(inout) :: gcAANA  ! This gridded component 
   type(ESMF_State),    intent(inout) :: impAANA ! Import state
   type(ESMF_State),    intent(inout) :: expAANA ! Export state
   type(ESMF_Clock),    intent(inout) :: CLOCK   ! The clock
   integer, optional,   intent(  out) :: RC      ! Error code

! !DESCRIPTION: This method sets up the import and export fields
!               for regridding and runs the GSI g.c.  
!
! !REVISION HISTORY:
!
!   05May2007 Todling Renamed background clock
!   14Sep2007 Todling Return when freq of background not satisfied
!   04Dec2009 Todling Allow gcm to run beyond last-background/analysis 
!
!EOPI
!-------------------------------------------------------------------------

! Locals

   type (MAPL_MetaComp), pointer  :: STATE
   character(len=ESMF_MAXSTR) :: IAm 
   integer                    :: STATUS, yy, mm, dd, hh, min, sec, mysecs
   character(len=ESMF_MAXSTR) :: COMP_NAME
   logical                    :: get_background
   type(ESMF_Config)          :: cf
   type(ESMF_Time)            :: currT   
   type(TGEOSaana), pointer   :: ThisIntSt
   type(Twrap)                :: wrap
   type(ESMF_Alarm)           :: Alarm
   type(ESMF_Time)            :: AlaTime

! Begin...

   if(analysis_is_done) then
      RETURN_(ESMF_SUCCESS)
   endif

   Iam = "Run"
   call ESMF_GridCompGet ( gcAANA, name=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // trim(Iam) 
  
   call ESMF_GridCompGetInternalState ( gcAANA, wrap, STATUS )
   VERIFY_(STATUS)
   ThisIntSt => wrap%intST

   call MAPL_GetObjectFromGC ( gcAANA, STATE, RC=STATUS)
   VERIFY_(STATUS)

   call ESMF_ClockGet ( clock, currTime=CurrT,  rc=STATUS ) ; VERIFY_(STATUS)
   call ESMF_TimeGet(currT, yy=YY, mm=MM, dd=DD, h=HH, m=MIN, s=SEC, rc=STATUS); VERIFY_(STATUS)
   if(MAPL_AM_I_ROOT()) then
      write(6,'(a,1x,i4.4,5(a,i2.2))') trim(Iam)//" The current TIME is ", &
                                       YY, "/", MM, "/", DD, " ", HH, ":", MIN, ":", SEC
   endif
   
   mysecs =  HH*3600 + MIN*60
   get_background = mod(mysecs,BKGfreq_sc).eq.0
   if ( .not. get_background ) then
        RETURN_(ESMF_SUCCESS)
   endif

   ! Regrid impAANA to   impGSI
   ! (halowidth=0)  -> (halowidth=1)
   ! -------------------------------  
   call ESMFL_Regrid (impAANA, ThisIntSt%impGSI, rc=status)
   VERIFY_(status)
   
   !   RUN GSI
   !   -------
   call ESMF_GridCompRun ( ThisIntSt%gcGSI, &
        importState=ThisIntSt%impGSI, &
        exportState=ThisIntSt%expGSI, &
        clock=clock, rc=status )
   VERIFY_(STATUS)

   call MAPL_StateAlarmGet(STATE, ALARM, NAME='last-bkg', RC=STATUS)
   VERIFY_(STATUS)
   analysis_is_done = ESMF_AlarmIsRinging(ALARM, rc=status)

   if (analysis_is_done) then   
     ! Regrid impGSI  to    expAANA
     ! (halowidth=1)  -> (halowidth=0)
     ! -------------------------------  
      call ESMFL_Regrid (ThisIntSt%expGSI, expAANA, rc=status)
      VERIFY_(status)  
   end if

   RETURN_(ESMF_SUCCESS)

   end subroutine Run

!-------------------------------------------------------------------
!BOPI

! !IROUTINE: Finalize -- Finalize method

! !INTERFACE:

   subroutine Finalize ( gcAANA, impAANA, expAANA, CLOCK, RC )

! !ARGUMENTS:

   type(ESMF_GridComp), intent(inout) :: gcAANA  ! This gridded component 
   type(ESMF_State),    intent(inout) :: impAANA ! Import state
   type(ESMF_State),    intent(inout) :: expAANA ! Export state
   type(ESMF_Clock),    intent(inout) :: CLOCK   ! The clock
   integer, optional,   intent(  out) :: RC      ! Error code

! !DESCRIPTION: 
 
!EOPI
!-------------------------------------------------------------------------
! Locals

   character(len=ESMF_MAXSTR)          :: IAm 
   integer                             :: STATUS
   character(len=ESMF_MAXSTR)          :: COMP_NAME
   type (MAPL_MetaComp),  pointer  :: STATE
   type(TGEOSaana), pointer            :: ThisIntSt
   type(Twrap)                         :: wrap

! Begin...

   call ESMF_GridCompGet ( gcAANA, name=COMP_NAME,  RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // "Finalize"
 
   call ESMF_GridCompGetInternalState ( gcAANA, wrap, rc )
   ThisIntSt => wrap%intST

!  Finalize GSI
!  ------------
   call ESMF_GridCompFinalize ( ThisIntSt%gcGSI, &
        importState=ThisIntSt%impGSI, &
        exportState=ThisIntSt%expGSI, &
        clock=clock, &
        rc=status )
   VERIFY_(STATUS) 

   call GEOSaana_Destroy_ ( ThisIntSt, rc=status)
   VERIFY_(STATUS)

   deallocate(ThisIntSt, stat=status)
   VERIFY_(STATUS)
   if(present(rc)) rc = STATUS

   RETURN_(ESMF_SUCCESS)

   end subroutine Finalize

!-------------------------------------------------------------------
!BOPI
!
! !IROUTINE: GEOSaana_Create -- 
 
! !INTERFACE:

   subroutine GEOSaana_Create_ ( gc, this, rc ) 

! !ARGUMENTS:

   type(ESMF_GridComp), intent(inout) :: gc
   type(TGEOSaana), pointer           :: this 
   integer, intent(out), optional     :: rc

! !DESCRIPTION:
!
!   27Feb2009 Todling - declare gc inout
!         
!EOPI
!----------------------------------------------------------------------

   type(ESMF_Config)  :: cf
   type(ESMF_Grid)    :: geosGrid
   integer            :: status
   character(len=*), parameter :: IAm = 'GEOSaana_Create'

! start

   call ESMF_GridCompGet ( gc,  config=cf, RC=STATUS )
   VERIFY_(STATUS)

! Create GSI coupling states
! --------------------------
   this%expGSI  = ESMF_StateCreate (name="GSI_Export",  &   
        stateIntent = ESMF_STATEINTENT_EXPORT, & 
        RC=STATUS ) 
   VERIFY_(STATUS)
   this%impGSI  = ESMF_StateCreate (name="GSI_Import",  &
        stateIntent = ESMF_STATEINTENT_IMPORT, &
        RC=STATUS )
   VERIFY_(STATUS)

!  Create GSI gridded component
!  ----------------------------
   this%gcGSI = ESMF_GridCompCreate(    &
        name='Atmos_Analysis', &
        configfile="GSI_GridComp.rc", &
!        gridcomptype=ESMF_ATM, &
        rc=status)
   VERIFY_(status)
   if(present(rc)) rc = STATUS
  
   RETURN_(ESMF_SUCCESS)

   end subroutine GEOSaana_Create_

!BOPI
!
! !IROUTINE: GEOSaana_Destroy_ --
 
! !INTERFACE:

   subroutine GEOSaana_Destroy_ ( this, rc )

! !ARGUMENTS:

   type(TGEOSaana), intent(inout) :: this
   integer, intent(out), optional :: rc

! !DESCRIPTION:
!
!EOPI
!----------------------------------------------------------------------

   integer :: status
   character(len=*), parameter :: IAm = 'GEOSaana_Destroy'

   call ESMF_GridCompDestroy  ( this%gcGSI, rc=status)
   VERIFY_(STATUS)
   call ESMF_StateDestroy  ( this%impGSI, rc=status )
   VERIFY_(STATUS)
   call ESMF_StateDestroy  ( this%expGSI, rc=status)
   VERIFY_(STATUS)
!   call ESMF_ClockDestroy  ( this%clock, rc=status)
!   VERIFY_(STATUS)
   if(present(rc)) rc = STATUS

   RETURN_(ESMF_SUCCESS)

   end subroutine GEOSaana_Destroy_

end module GEOSaana_GridCompMod
