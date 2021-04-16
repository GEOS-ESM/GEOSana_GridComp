#include "GEOS_Generic.h"
#define USING_GEOS_FILE_MODEL

!-------------------------------------------------------------------------
!BOP

! !MODULE: 

! GEOSaana_GridCompMod -- A Module to provide coupling between GEOS and GSI

! !INTERFACE:

   module GEOSaana_GridCompMod

   use ESMF
   use GEOS_Mod
   use ESMF_CfioMod
   use GEOS_CfioMod 
   use m_StrTemplate         ! grads style templates
   use GSI_GridCompMod, only: GSI_SetServices => GSI_GridCompSetServices

   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:

   public SetServices

! Internal state
! --------------
   type TGEOSaana
     private
     type(ESMF_GridComp) :: gsiGC              ! GSI gridded component
     type(ESMF_State)    :: GSIimp, GSIexp     ! import/export for GSI
     type(ESMF_Clock)    :: clock          

! include this ifdef temporarily to signify that the GEOS g.c. is a file model
! The entire block should be removed once we can run a single executable.

#ifdef USING_GEOS_FILE_MODEL
     type(ESMF_Bundle)   :: GEOSdynBUN1   ! bundle to hold GEOS export UPA fields
     type(ESMF_Bundle)   :: GEOSdynBUN2   !  "          "       import UPA   "

#endif

   end type TGEOSaana

! Wrapper for extracting internal state
! -------------------------------------
   type Twrap
     type(TGEOSaana), pointer :: intSt
   end type Twrap

! !DESCRIPTION: This gridded component provides a layer for coupling
!               between GEOS-5 and GSI. 
!
! This component is necessary because GEOS-5 and GSI have different grids 
! and since MAPL does not provide a coupling mechanism between components
! with different grids the coupling must be performed "manually". Hence the 
! need for this component.

! !REVISION HISTORY:
!
!  23Oct2006  Cruz  Initial code.
!
!EOP
!-------------------------------------------------------------------------

   CONTAINS

!BOP

!-------------------------------------------------------------------------
! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

   subroutine SetServices ( GC, RC )

! !ARGUMENTS:

   type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
   integer, optional                  :: RC  ! return code

! ! DESCRIPTION: This function sets the IRF services,and creates the GSI 
!       sub-components (g.c., grid, and coupling states).

! !REVISION HISTORY:
!
!  23Oct2006  Cruz  Initial code.
!
!EOP
!-------------------------------------------------------------------------

! Locals

   character(len=ESMF_MAXSTR) :: IAm
   integer                    :: STATUS
   character(len=ESMF_MAXSTR) :: COMP_NAME
   type(TGEOSaana), pointer   :: ThisIntSt
   type(Twrap)                :: wrap

! Begin...

   call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // 'SetServices'

! Register services for this component
! ------------------------------------

   call ESMF_GridCompSetEntryPoint ( GC, ESMF_SETINIT, Initialize,   &
        ESMF_SINGLEPHASE, STATUS )
   VERIFY_(STATUS)
   call ESMF_GridCompSetEntryPoint ( GC, ESMF_SETRUN,  Run       ,   &
        ESMF_SINGLEPHASE, STATUS )
   VERIFY_(STATUS)
   call ESMF_GridCompSetEntryPoint ( GC, ESMF_SETFINAL,  Finalize       ,   &
        ESMF_SINGLEPHASE, STATUS )
   VERIFY_(STATUS)

! Create the GSI sub-components
! -----------------------------
   allocate(ThisIntSt, stat=status)
   VERIFY_(STATUS)
   call GEOSaana_Create_( gc, ThisIntSt, rc=status)
   VERIFY_(STATUS)

! Save the "internal" state
! -------------------------
   wrap%intSt => ThisIntSt
   call ESMF_GridCompSetInternalState(GC, wrap, rc)
   
!  Register the GSI gridded component
!  ----------------------------------
   call ESMF_GridCompSetServices ( ThisIntSt%gsiGC, GSI_SetServices, status )
   VERIFY_(STATUS)

   RETURN_(ESMF_SUCCESS)
  
   end subroutine SetServices

!-------------------------------------------------------------------------
! !IROUTINE: Initialize -- Initialize method

! !INTERFACE:

   subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

   type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
   type(ESMF_State),    intent(inout) :: IMPORT ! Import state
   type(ESMF_State),    intent(inout) :: EXPORT ! Export state
   type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
   integer, optional,   intent(  out) :: RC     ! Error code

! ! DESCRIPTION: This function initializes the AANA gridded component.
!       It is a wrapper to the initialization of the GSI g.c.
 

! !REVISION HISTORY:
!
!  23Oct2006  Cruz  Initial code.
!
!EOP
!-------------------------------------------------------------------------

   type (GEOS_GenericState),  pointer  :: STATE
   character(len=ESMF_MAXSTR) :: IAm
   integer                    :: STATUS
   character(len=ESMF_MAXSTR) :: COMP_NAME
   type(TGEOSaana), pointer   :: ThisIntSt
   type(Twrap)                :: wrap
   type(ESMF_Config)          :: cf
   type(ESMF_Alarm)           :: Final_End_Time_Alarm

! Begin...

   call ESMF_GridCompGet ( GC, name=COMP_NAME, config=cf, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // "Initialize"

! Retrieve internal state

   call ESMF_GridCompGetInternalState ( GC, wrap, rc )
   ThisIntSt => wrap%intST

!  Initialize GSI
!  --------------
   if ( GEOS_AM_I_ROOT() ) print *, 'Initialize GSI'      
   call ESMF_GridCompInitialize ( ThisIntSt%gsiGC, &
                                  ThisIntSt%GSIimp, &
                                  ThisIntSt%GSIexp, &
                                  clock, &
                                  rc=status )
   VERIFY_(STATUS) 

   RETURN_(ESMF_SUCCESS)

   end subroutine Initialize

!-------------------------------------------------------------------------
! !IROUTINE: Run -- Run method

! !INTERFACE:

   subroutine Run ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

   type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
   type(ESMF_State),    intent(inout) :: IMPORT ! Import state
   type(ESMF_State),    intent(inout) :: EXPORT ! Export state
   type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
   integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: This method sets up the import and export fields
!               for regridding and runs the GSI g.c.  

! !REVISION HISTORY:
!
!  23Oct2006  Cruz  Initial code.
!
!EOP
!-------------------------------------------------------------------------

! Locals

   type (GEOS_GenericState),  pointer  :: STATE
   character(len=ESMF_MAXSTR) :: IAm 
   integer                    :: STATUS
   character(len=ESMF_MAXSTR) :: COMP_NAME
   type (ESMF_Alarm)          :: ALARM
   type(ESMF_Config)          :: cf
   type(ESMF_Time)            :: currT   
   integer                    :: numBKGs,BKGinterval 
   integer                    :: yy, mm, dd, hh, min, sec
   integer                    :: nymd, nhms
   integer                    :: n
   character(len=ESMF_MAXSTR) :: bkgpath
   character(len=ESMF_MAXSTR) :: DateStamp
   character(len=ESMF_MAXSTR) :: expid, aname
   character(len=ESMF_MAXSTR) :: upa_name, upa_tmpl, upa_file
   character(len=ESMF_MAXSTR) :: dynfile
   type(TGEOSaana), pointer   :: ThisIntSt
   type(Twrap)                :: wrap

! Begin...

   Iam = "Run"
   call ESMF_GridCompGet ( GC, name=COMP_NAME, config=cf, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // trim(Iam) 
   
   call ESMF_GridCompGetInternalState ( GC, wrap, rc )
   ThisIntSt => wrap%intST

! include this ifdef temporarily to signify that the GEOS g.c. is a file model
! The entire block should be removed once we can run a single executable.

#ifdef USING_GEOS_FILE_MODEL

   CALL ESMF_ConfigGetAttribute(cf, numBKGs,   label = 'nfldsig:',   &
        rc=status); VERIFY_(status)
   call ESMF_ConfigGetAttribute( cf, expid, label ='expid:', rc = status )
   VERIFY_(status)
   call ESMF_ConfigGetAttribute( cf, upa_name, label ='upa_name:', rc = status )
   VERIFY_(status)
   call ESMF_ConfigGetAttribute( cf, upa_tmpl, label ='upa_tmpl:', rc = status )
   VERIFY_(status)
   call ESMF_ConfigGetAttribute( cf, bkgpath, label ='bkgpath:', rc = status )
   VERIFY_(status)

#endif

   do n=1, numBKGs

      call ESMF_ClockGet ( clock, currTime=CurrT,  rc=STATUS ) ; VERIFY_(STATUS)
      call ESMF_TimeGet(currT, yy=YY, mm=MM, dd=DD, h=HH, m=MIN, s=SEC, rc=rc) 
      if(GEOS_AM_I_ROOT()) print *, "The FCSTTIME is ", YY, "/", MM, "/", DD, &
        " ", HH, ":", MIN, ":", SEC

! include this ifdef temporarily to signify that the GEOS g.c. is a file model
! The entire block should be removed once we can run a single executable.

#ifdef USING_GEOS_FILE_MODEL

      ! create filenames based on date stamp
      call get_DateStamp ( clock, DateStamp, expid,  &
           rc=status )
      VERIFY_(STATUS) 
      read(DateStamp( 1: 8),'(i8.8)') nymd
      read(DateStamp(10:15),'(i6.6)') nhms
      if ( GEOS_AM_I_ROOT() ) print *, 'NYMD, NHMS = ',nymd,nhms
      upa_file = trim(bkgpath)//trim(expid)//"."// trim(upa_name) // "." //trim(upa_tmpl)
      call StrTemplate ( dynfile, upa_file, 'GRADS', xid=expid, nymd=nymd, nhms=nhms, &
           stat=status )
      VERIFY_(STATUS) 
      if ( GEOS_AM_I_ROOT() ) print *, 'dynfile: ',trim(dynfile)

      !  Read GEOS EXPORT  from files
      !  ----------------------------
      call ESMF_ioRead  ( dynfile, currT, ThisIntSt%GEOSdynBUN1, rc=status, &
           verbose=.true., force_regrid=.true., &
           only_vars='ps')
      VERIFY_(status)
      
      !  Create bundle from file to hold GEOS IMPORT state fields
      !  --------------------------------------------------------
      call  ESMF_ioRead ( dynfile, currT, ThisIntSt%GEOSdynBUN2, rc=status, &
           verbose=.true., noread=.true., &
           only_vars='ps')
      VERIFY_(status)

      !   Extract fields from bundles and use them to populate states
      !   -----------------------------------------------------------
      call ESMFL_Bundle2State (ThisIntSt%GEOSdynBUN1, EXPORT)
      call ESMFL_Bundle2State (ThisIntSt%GEOSdynBUN2, IMPORT)  

#endif

      !   Regrid EXPORT       to GSIstate 
      !          halowidth=0  -> halowidth=1
      !   ----------------------------  
      if(GEOS_AM_I_ROOT()) then
        call ESMF_StatePrint(EXPORT)
      end if
      call ESMFL_Regrid (EXPORT, ThisIntSt%GSIimp, rc=status)
      VERIFY_(status)  
      if(GEOS_AM_I_ROOT()) then
        call ESMF_StatePrint(ThisIntSt%GSIimp)
      end if

      !   RUN GSI
      !   -------
      call ESMF_GridCompRun ( ThisIntSt%gsiGC, &
                              ThisIntSt%GSIimp, &
                              ThisIntSt%GSIexp, &
                              clock, rc=status )
      VERIFY_(STATUS)

   end do ! numBKGs

!   Regrid GSIstate     to IMPORT 
!          halowidth=1  -> halowidth=0
!   ----------------------------
   call ESMFL_Regrid (ThisIntSt%GSIexp, IMPORT, rc=status)
   VERIFY_(status)  

   RETURN_(ESMF_SUCCESS)

   end subroutine Run

!-------------------------------------------------------------------------
! !IROUTINE: Finalize -- Finalize method

! !INTERFACE:

   subroutine Finalize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

   type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
   type(ESMF_State),    intent(inout) :: IMPORT ! Import state
   type(ESMF_State),    intent(inout) :: EXPORT ! Export state
   type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
   integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: 
 

!! !REVISION HISTORY:
!
!  23Oct2006  Cruz  Initial code.
!
!EOP
!-------------------------------------------------------------------------
! Locals

   character(len=ESMF_MAXSTR)          :: IAm 
   integer                             :: STATUS
   character(len=ESMF_MAXSTR)          :: COMP_NAME
   type (GEOS_GenericState),  pointer  :: STATE
   type(TGEOSaana), pointer            :: ThisIntSt
   type(Twrap)                         :: wrap

! Begin...

   call ESMF_GridCompGet ( GC, name=COMP_NAME,  RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // "Finalize"
 
   call ESMF_GridCompGetInternalState ( GC, wrap, rc )
   ThisIntSt => wrap%intST

!  Finalize GSI
!  ------------
   if ( GEOS_AM_I_ROOT() ) print *, 'Finalize GSI'      
   call ESMF_GridCompFinalize ( ThisIntSt%gsiGC, &
        ThisIntSt%GSIimp, &
        ThisIntSt%GSIexp, &
        clock, &
        rc=status )
   VERIFY_(STATUS) 

   call GEOSaana_Destroy_ ( ThisIntSt, rc=status)
   VERIFY_(STATUS)

   deallocate(ThisIntSt, stat=status)
   VERIFY_(STATUS)

   RETURN_(ESMF_SUCCESS)

   end subroutine Finalize

!----------------------------------------------------------------------
!
! !IROUTINE: GEOSaana_Create -- 
 
! !INTERFACE:

   subroutine GEOSaana_Create_ ( gc, this, rc ) 

! !ARGUMENTS:

   type(ESMF_GridComp), intent(in)       :: gc
   type(TGEOSaana), pointer        :: this 
   integer, intent(out), optional  :: rc

! !DESCRIPTION:
!         
!----------------------------------------------------------------------

   type(ESMF_Config)  :: cf
   type(ESMF_VM)      :: vm
   type(ESMF_Grid)    :: geosGrid
   integer            :: status
   character(len=*), parameter :: IAm = 'GEOSaana_Create'

! start

   call ESMF_GridCompGet ( gc,  config=cf, RC=STATUS )
   VERIFY_(STATUS)

! Create  (sub-components use same VM)
! -----------------------------------   
   call ESMF_VMGetGlobal ( VM,  RC=STATUS )
   VERIFY_(STATUS)

! Create GSI coupling states
! --------------------------
   this%GSIexp  = ESMF_StateCreate (name="GSI_Export",  &   
        stateintent = ESMF_STATEINTENT_EXPORT, & 
        RC=STATUS ) 
   VERIFY_(STATUS)
   this%GSIimp  = ESMF_StateCreate (name="GSI_Import",  &
        stateIntent = ESMF_STATEINTENT_IMPORT, &
        RC=STATUS )
   VERIFY_(STATUS)

!  Create GSI gridded component
!  ----------------------------
   this%gsiGC = ESMF_GridCompCreate(vm, &
        name='GSI Gridded Component', &
        configfile="GSI_GridComp.rc", &
!        gridcomptype=ESMF_ATM, &
        rc=status)
   VERIFY_(status)
  
! include this ifdef temporarily to signify that the GEOS g.c. is a file model
! The entire block should be removed once we can run a single executable.

#ifdef USING_GEOS_FILE_MODEL

!   Create GEOS grid 
!   ----------------
    GeosGrid = GEOSGridCreate_ ( vm, cf, rc=status )
    VERIFY_(status) 

!  Create empty bundles to hold GEOS fields
!  ----------------------------------------

   This%GEOSdynBUN1 = ESMF_BundleCreate ( name='geos_dyn_bundle1', &
        grid=geosGrid, rc=status )
   VERIFY_(status)
   This%GEOSdynBUN2 = ESMF_BundleCreate ( name='geos_dyn_bundle2', &
        grid=geosGrid, rc=status )
   VERIFY_(status)

#endif

   RETURN_(ESMF_SUCCESS)

   end subroutine GEOSaana_Create_

!----------------------------------------------------------------------
!
! !IROUTINE: GEOSaana_Destroy_ --
 
! !INTERFACE:

   subroutine GEOSaana_Destroy_ ( this, rc )

! !ARGUMENTS:

   type(TGEOSaana), intent(inout) :: this
   integer, intent(out), optional :: rc

! !DESCRIPTION:
!
!----------------------------------------------------------------------

   integer :: status
   character(len=*), parameter :: IAm = 'GEOSaana_Destroy'

   call ESMF_GridCompDestroy  ( this%GSIgc, rc=status)
   VERIFY_(STATUS)
   call ESMF_StateDestroy  ( this%GSIimp, rc=status )
   VERIFY_(STATUS)
   call ESMF_StateDestroy  ( this%GSIexp, rc=status)
   VERIFY_(STATUS)
!   call ESMF_ClockDestroy  ( this%clock, rc=status)
!   VERIFY_(STATUS)

! include this ifdef temporarily to signify that the GEOS g.c. is a file model
! The entire block should be removed once we can run a single executable.

#ifdef USING_GEOS_FILE_MODEL

   call ESMF_BundleDestroy ( This%GEOSdynBUN1, rc=STATUS )
   VERIFY_(STATUS)
   call ESMF_BundleDestroy ( This%GEOSdynBUN2, rc=STATUS )
   VERIFY_(STATUS)

#endif

   RETURN_(ESMF_SUCCESS)

   end subroutine GEOSaana_Destroy_

! include this ifdef temporarily to signify that the GEOS g.c. is a file model
! The entire block should be removed once we can run a single executable.

#ifdef USING_GEOS_FILE_MODEL

!-------------------------------------------------------------------
   function GEOSGridCreate_ ( vm, cf, rc) result(grid)
!-------------------------------------------------------------------

    type (ESMF_VM),    intent(IN   ) :: VM
    type(ESMF_Config), intent(INOUT) :: cf
    integer, optional, intent(OUT)   :: rc
    type (ESMF_Grid)                 :: grid

! Local vars
    integer                                 :: status
    character(len=ESMF_MAXSTR), parameter   :: IAm='GEOSGridCreate'

    type (ESMF_DELayout)                 :: layout
    integer                         :: IM,JM,LM
    integer                         :: L
    integer                         :: NX, NY
    integer, allocatable            :: IMXY(:), JMXY(:)
    character(len=ESMF_MAXSTR)      :: gridname
    real(ESMF_KIND_R8)              :: minCoord(3)
    real(ESMF_KIND_R8)              :: deltaX, deltaY, deltaZ
    real                            :: LON0, LAT0

    real :: pi, d2r

! grid create

   call ESMF_ConfigGetAttribute( cf, im, label ='GEOS IM:', rc = status )
   VERIFY_(status)
   call ESMF_ConfigGetAttribute( cf, jm, label ='GEOS JM:', rc = status )
   VERIFY_(status)
   call ESMF_ConfigGetAttribute( cf, LM,       label ='GEOS KM:', rc = status )
   VERIFY_(status)
   call ESMF_ConfigGetAttribute( cf, NX,       label ='NX:', rc = status )
   VERIFY_(status)
   call ESMF_ConfigGetAttribute( cf, NY,       label ='NY:', rc = status )
   VERIFY_(status)

     pi  = 4.0 * atan ( 1.0 )
    d2r  = pi / 180.
    LON0 = -180  * d2r
    LAT0 = -90.0 * d2r
    deltaX = 2.0*pi/IM
    deltaY = pi/(JM-1)
    deltaZ = 1.0

! Get the IMXY and JMXY vectors
! -----------------------------
    allocate( imxy(0:nx-1), jmxy(0:ny-1), stat=status)
    VERIFY_(STATUS)
    call GEOS_GET_LOCAL_DIMS ( IM, imxy, nx ) 
    call GEOS_GET_LOCAL_DIMS ( JM, jmxy, ny )
    if ( GEOS_Am_I_Root() ) then
       print *, 'GEOS im,jm,km  = ',IM, JM,LM
       print *, 'nx : imxy = ', nx, ' : ', imxy
       print *, 'ny : jmxy = ', ny, ' : ', jmxy
    endif
    
! Define South-West Corner of First Grid-Box
! ------------------------------------------
    minCoord(1) = LON0 - deltaX/2
    minCoord(2) = LAT0 - deltaY/2
    minCoord(3) = deltaZ/2.
    
    layout = ESMF_DELayoutCreate(vm, deCountList=(/NX, NY/), rc=status)
    VERIFY_(STATUS)
    
    grid = ESMF_GridCreateHorzLatLonUni(         &
         counts = (/IM, JM/),                    &
         minGlobalCoordPerDim=minCoord(1:2),     &
         deltaPerDim=(/deltaX, deltaY /),        &
         horzStagger=ESMF_Grid_Horz_Stagger_A,   &
         periodic=(/ESMF_TRUE, ESMF_FALSE/),     &
         name='GEOS Grid', rc=status)
    VERIFY_(STATUS)

    call ESMF_GridAddVertHeight(grid,            &
         delta=(/(deltaZ, L=1,LM) /),            &
         vertStagger=ESMF_GRID_VERT_STAGGER_TOP, &
         rc=status)
    VERIFY_(STATUS)
   
    call ESMF_GridDistribute(grid,               &
         deLayout=layout,                        &
         countsPerDEDim1=imxy,                   &
         countsPerDEDim2=jmxy,                   &
         rc=status)
    VERIFY_(STATUS)
   
    deallocate(imxy, jmxy, stat=status)
    VERIFY_(STATUS)
    
    RETURN_(STATUS)
    
   end function GEOSGridCreate_

! taken from GEOShistory:

!-------------------------------------------------------------------
   subroutine get_DateStamp (clock, DateStamp, expid, offset, rc)
!-------------------------------------------------------------------
    type (ESMF_Clock)                 :: clock
    character(len=ESMF_MAXSTR)        :: DateStamp
    character(len=ESMF_MAXSTR)        :: expid
    type(ESMF_TimeInterval), optional :: offset
    integer, optional                 :: rc

! locals
    type(ESMF_Time)                   :: currentTime
    type(ESMF_TimeInterval)           :: TimeStep
    character(len=ESMF_MAXSTR)        :: TimeString
    integer                           :: secs
    character(len=ESMF_MAXSTR)        :: TimeStamp
    character                         :: String(ESMF_MAXSTR)

    character*4 year
    character*2 month
    character*2 day
    character*2 hour
    character*2 minute
    character*2 second

    equivalence ( string(01),TimeString )
    equivalence ( string(01),year       )
    equivalence ( string(06),month      )
    equivalence ( string(09),day        )
    equivalence ( string(12),hour       )
    equivalence ( string(15),minute     )
    equivalence ( string(18),second     )

    call ESMF_ClockGet (clock, currTime=currentTime, rc=rc)
    if (present(offset)) then
       currentTime = currentTime - offset
    end if
    call ESMF_TimeGet  (currentTime, timeString=TimeString, rc=rc)
    call ESMF_ClockGet (clock, timeStep=TimeStep, rc=rc)
    call ESMF_TimeIntervalGet (TimeStep, S=secs, rc=rc)

    DateStamp = year // month // day // '_' // hour // minute // second // 'z'
    TimeStamp = '   Date: ' // year // '/' // month // '/' // day
    TimeStamp = trim(TimeStamp) // '  Time: ' // timestring(12:19)

    if (.not. present(OFFSET)) then
       call WRITE_PARALLEL ( 'Expid: ' // trim(expid) // trim(TimeStamp) )
    endif

   end subroutine get_DateStamp

#endif


end module GEOSaana_GridCompMod
