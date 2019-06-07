! $Id$

#include "MAPL_Generic.h"

!=============================================================================

module GEOS_PChemPertGridCompMod

!BOP

! !MODULE: GEOS_PChemGridCompMod

! !DESCRIPTION: 
!   {\tt GEOS\_PChem} is a proxy component for Aerochem that implements the 
!    specification or simple parameterization of the Aerochem species. It works
!    on three types of species: chemical species (oxygen, nitrous oxide, CFC-11,
!    CFC-12, CFC-22, methane, water vapor), diagnostic species (age-of-air),
!    and aerosols (arbitrary).
!
!    Each of the chemical species can be treated in one
!    of two ways: parameterized prediction from tabled zonally-symmetric production and 
!    loss (P-L) data, or specifycation from zonally-symmetric values (see Resources section
!    for how to control thios behavior).  A single flat file containing 
!    both the P-L and climatology data {\it must} be provided (see Resources section).
!    Aerosols are set to 3-dimensional climatological values.
!    The ``age-of-air'' is predicted by setting the surface values of this tracer
!    to the to zero and advancing other levels by dt.
!    All of these quantities except water vapor are INTERNAL state variables of
!    {\tt GEOS\_PChem}. Water vapor is assumed to be a Friendly Import and {\tt GEOS\_PChem}
!    leaves it unmodified below the tropopause, or the 200 hPa level if the tropopause
!    is below this level. 
!
!
!    For chemical species, the production rate is tabled directly.
!    For the loss, a rate coefficient is tabled.  Using Odd-oxygen $O_x$ as an example,
!    the species are updated as follows:
!    $$
!    \frac{\partial O_x}{\partial t} = 
!        \dot{Q}_{o} - \kappa_{o} O_x
!    $$
!    where $O_x$ is the specific mass of odd oxygen, $ \dot{Q}_{o}$ is the
!    odd oxygen production rate, $\kappa_{o}$ is the tabled rate 
!    coefficient for odd oxygen loss. This is finite differenced in time as:
!    $$
!      O_x^{n+1} =  \frac{O_x^{n} + \Delta t \dot{Q}_{o} }{1 + \Delta t \kappa_{o}}
!    $$
!
!    The component reads the 12 monthly tables of the zonally averaged climatology
!    of the production rates and loss frequencies 
!    and interpolates them to the meridional locations of the natural grid in Initialize.
!    These are saved in the private internal state, which is static. The five species
!    are updated and kept in INTERNAL, an ESMF state attached to the
!    GEOS GENERIC object in the component. If no restart is specified for the
!    INTERNAL, the species are initialized to zero.
!
!    Ozone is diagnosed from $O_x$ by assuming that it accounts for all
!    $O_x$ at pressures greater than 100 Pa (1 hPa) during the day and at all 
!    pressures at night. For those daylit cells where pressures are less than 1 hPa, 
!    we assume that the ozone fraction in $O_x$ decreases exponentially with decreasing
!    pressure.
!    
!    Aerosols are read from 3-dimensional data files that have to be on model levels
!    but may be on any regular lat-lon grid. These are hdf files and horizontal
!    interpolation is done through CFIO. The aerosols are read into a bundle that
!    is exported, and the number and names of the aerosols in the bundle are set
!    from the CFIO file. 
!\newline
!
!  !RESOURCES:
!     >> RUN\_DT:  none    real    seconds
!        Heartbeat. {\tt GEOS\_PChem} is called all the time.
!     >> pchem\_clim:                'pchem\_clim.dat'  string    none 
!           Zonally-symmetric chemistry data flat file.  
!     >> AEROCLIM:                  no\_aerosols       string    none 
!           Aerosol monthly climatology hdf file.  
!     >> AEROCLIMDEL:               no\_aerosols       string    kg kg-1 
!           Aerosol month to month differences hdf file.  
!     >> AEROCLIMYEAR:              2002              integer   year 
!           Year used in timestamp of the two hdf aerosol files.  
!     >> {\it{name}}\_FIXED\_VALUE:    use\_file          real      pppv 
!           Constant value at which to fix chemical species {\it{name}}.  If not specified, 
!         the file data will be used. If specified, {\it{name}}\_RELAXTIME: is ignored. 
!          {\it{name}} can be any of OX, CH4, N2O, CFC11, CFC12, CFC22.  
!     >> {\it{name}}\_RELAXTIME:       0.0                real       seconds 
!         Timescale of relaxation to climatology on file for chemical species {\it{name}}. 
!         For values $<= 0$, the P-L parameterization will be used. To hold at the file's 
!         zonally-symmetric climatology, use a small positive number. 
!         {\it{name}} can be any of OX, CH4, N2O, CFC11, CFC12, CFC22, H2O.  
!     >> {\it{name}}\_PCRIT:          1.e+16             real      Pa 
!           Pressure of level above which the relaxation to climatology is done. This is 
!         ignored if {\it{name}}\_RELAXTIME: is ignored or if {\it{name}}\_RELAXTIME: is $<= 0$. 
!         {\it{name}} can be any of OX, CH4, N2O, CFC11, CFC12, CFC22.  
!     >> {\it{name}}\_DELP:           1.e-16             real      Pa 
!           Pressure interval over which the relaxation to climatology is ramped-in. This is 
!         ignored if {\it{name}}\_RELAXTIME: is ignored or if {\it{name}}\_RELAXTIME: is $<= 0$ 
!         {\it{name}} can be any of OX, CH4, N2O, CFC11, CFC12, CFC22.  
!     >> {\it{name}}\_FRIENDLIES:     self              string    none 
!           String of colon separated component names to which this species is Friendly. 
!         {\it{name}} can be any of OX, CH4, N2O, CFC11, CFC12, CFC22  
!     >> AOA\_FRIENDLIES:           'DYNAMICS:TURBULENCE'  string     none 
!           String of colon separated component names to which Age-of-Air is Friendly.  
!
! !USES:

  use ESMF
  use MAPL_Mod
  use Chem_Mod
  use ESMF_CFIOFileMOD
  use MAPL_CFIOMOD
  
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!EOP

!=============================================================================

  type T_Pchem_STATE
     private
     integer                             :: NLATS, NLEVS
     real, pointer, dimension(:)         :: LATS => null()
     real, pointer, dimension(:)         :: LEVS => null()
     real, pointer, dimension(:,:,:,:,:) :: MNPL => null()
     integer                             :: OX       = 1
     integer                             :: N2O      = 2
     integer                             :: CFC11    = 3
     integer                             :: CFC12    = 4
     integer                             :: CH4      = 5
     integer                             :: CFC22    = 6
     integer                             :: H2O      = 7
     integer                             :: NSPECIES = 7
  end type T_Pchem_STATE

  type Pchem_WRAP
     type (T_Pchem_STATE), pointer :: PTR
  end type Pchem_WRAP

contains

!=============================================================================

!BOP

! !IROUTINE: SetServices

! !DESCRIPTION: Sets Initialize and Run services. 
! \newline
!

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

! Local derived type aliases

    type (ESMF_Config          )            :: CF

    type (T_PCHEM_STATE), pointer           :: PCHEM_state 
    type (T_PCHEM_STATE)                    :: DUMMY
    type (PCHEM_wrap)                       :: wrap
    character(len=ESMF_MAXSTR)              :: OXFRIENDLY
    character(len=ESMF_MAXSTR)              :: N2OFRIENDLY
    character(len=ESMF_MAXSTR)              :: CFC11FRIENDLY
    character(len=ESMF_MAXSTR)              :: CFC12FRIENDLY
    character(len=ESMF_MAXSTR)              :: CFC22FRIENDLY
    character(len=ESMF_MAXSTR)              :: CH4FRIENDLY
    character(len=ESMF_MAXSTR)              :: AOAFRIENDLY
    character(len=ESMF_MAXSTR)              :: AEROFRIENDLY

    character(len=ESMF_MAXSTR)              :: FRIENDLIES
    character(len=ESMF_MAXSTR)              :: providerName

    INTEGER :: n

    type(Chem_Registry)                     :: chemReg   

    character(len=*), parameter             :: thisrc='Chem_Registry_apert.rc'

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

!   Start by loading the Chem Registry
!   ----------------------------------
    chemReg = Chem_RegistryCreate ( STATUS, rcfile=thisrc )
    VERIFY_(STATUS)

!   If not doing PChem, use GEOS Generic stubs from this point on
!   -------------------------------------------------------------
    if ( .NOT. chemReg%doing_PC ) then
       call MAPL_GenericSetServices ( GC, RC=STATUS )
       VERIFY_(STATUS)
       call Chem_RegistryDestroy ( chemReg, RC=STATUS )
       VERIFY_(STATUS)
       if (MAPL_AM_I_ROOT()) & 
           print *, trim(Iam)//': not ACTIVE, defaulting to GG stubs...'
       RETURN_(ESMF_SUCCESS)
    end if       

! Set the Initialize and Run entry point
! --------------------------------------


! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate( PCHEM_state, stat=STATUS )
    VERIFY_(STATUS)

    WRAP%PTR => PCHEM_STATE
    PCHEM_STATE = DUMMY
 
! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC, 'Pchem_state', WRAP, STATUS )
    VERIFY_(STATUS)

! Get the configuration
! ---------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

! Set the state variable specs.
! -----------------------------

     FRIENDLIES = trim(COMP_NAME)

     call ESMF_ConfigGetAttribute(CF, OXFRIENDLY, Label='OX_FRIENDLIES:'      ,&
                                  default=FRIENDLIES                          ,&
                                                                     RC=STATUS )
     VERIFY_(STATUS)

     call ESMF_ConfigGetAttribute(CF, N2OFRIENDLY, Label='N2O_FRIENDLIES:'    ,&
                                  default=FRIENDLIES                          ,&
                                                                     RC=STATUS )
     VERIFY_(STATUS)

     call ESMF_ConfigGetAttribute(CF, CFC11FRIENDLY, Label='CFC11_FRIENDLIES:',&
                                  default=FRIENDLIES                          ,&
                                                                     RC=STATUS )
     VERIFY_(STATUS)

     call ESMF_ConfigGetAttribute(CF, CFC12FRIENDLY, Label='CFC12_FRIENDLIES:',&
                                  default=FRIENDLIES                          ,&
                                                                     RC=STATUS )
     VERIFY_(STATUS)

     call ESMF_ConfigGetAttribute(CF, CFC22FRIENDLY, Label='CFC22_FRIENDLIES:',&
                                  default=FRIENDLIES                          ,&
                                                                     RC=STATUS )
     VERIFY_(STATUS)

     call ESMF_ConfigGetAttribute(CF, CH4FRIENDLY, Label='CH4_FRIENDLIES:'    ,&
                                  default=FRIENDLIES                          ,&
                                                                     RC=STATUS )
     VERIFY_(STATUS)

     call ESMF_ConfigGetAttribute(CF, AOAFRIENDLY, Label='AOA_FRIENDLIES:'    ,&
                                  default=FRIENDLIES                          ,&
                                                                     RC=STATUS )
     VERIFY_(STATUS)

     call ESMF_ConfigGetAttribute(CF, AEROFRIENDLY, Label='AERO_FRIENDLIES:'  ,&
                                  default=FRIENDLIES                          ,&
                                                                     RC=STATUS )
     VERIFY_(STATUS)


!BOP

! !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'PLE',                                       &
        LONG_NAME  = 'air_pressure',                              &
        UNITS      = 'Pa',                                        &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,                          &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME ='Q',                                          &
        LONG_NAME  ='specific_humidity',                          &
        UNITS      ='1',                                          &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,                        &  
                                                       RC=STATUS  )
     VERIFY_(STATUS)                                                                          

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'TROPP',                             &
        LONG_NAME          = 'tropopause_pressure',               &
        UNITS              = 'Pa',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


! !INTERNAL STATE:

! For odd-oxygen only:  Ox is the second member of the analysis bundle, and 
! if the ANALYSIS_OX_PROVIDER is PCHEM, Ox must be friendly to "ANALYSIS".
! If Ox is not in TRANA, AGCM will fail after it counts the number of members 
! in TRANA [ASSERT(NumFriendly == 2)].
! --------------------------------------------------------------------------

     CALL ESMF_ConfigGetAttribute(CF, providerName, Default="PCHEM", &
                                  Label="ANALYSIS_OX_PROVIDER:", RC=STATUS )
     VERIFY_(STATUS)
     
     IF( providerName == "PCHEM" .AND. (INDEX(OXFRIENDLY,"ANALYSIS") == 0) ) THEN
      IF(MAPL_AM_I_ROOT()) THEN
       PRINT *," "
       PRINT *,TRIM(Iam)//": OX_FRIENDLIES in AGCM.tmpl invalid."
       PRINT *,"      You must at least specificy:  OX_FRIENDLIES:ANALYSIS"
       PRINT *,"      You may also add other component if desired. Example: DYNAMICS:ANALYSIS."
       PRINT *," "
      END IF
      STATUS = 1
      VERIFY_(STATUS)
     END IF

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'OX',                                &
        LONG_NAME          = 'odd_oxygen_mixing_ratio',           &
        UNITS              = 'pppv',                              &
        FRIENDLYTO         = OXFRIENDLY,                          &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'N2O',                               &
        LONG_NAME          = 'nitrous_oxide_mixing_ratio',        &
        UNITS              = 'pppv',                              &
        FRIENDLYTO         = N2OFRIENDLY,                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CFC11',                             &
        LONG_NAME          = 'CFC11_(CCl3F)_mixing_ratio',        &
        UNITS              = 'pppv',                              &
        FRIENDLYTO         = CFC11FRIENDLY,                       &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CFC12',                             &
        LONG_NAME          = 'CFC12_(CCl2F2)_mixing_ratio',       &
        UNITS              = 'pppv',                              &
        FRIENDLYTO         = CFC12FRIENDLY,                       &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CFC22',                             &
        LONG_NAME          = 'CFC22_mixing_ratio',                &
        UNITS              = 'pppv',                              &
        FRIENDLYTO         = CFC22FRIENDLY,                       &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CH4',                               &
        LONG_NAME          = 'methane_mixing_ratio',              &
        UNITS              = 'pppv',                              &
        FRIENDLYTO         = CH4FRIENDLY,                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'AOA',                               &
        LONG_NAME          = 'age_of_air',                        &
        UNITS              = 'days',                              &
        FRIENDLYTO         = AOAFRIENDLY,                         &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'AERO',                              &
        LONG_NAME          = 'aerosol_mass_mixing_ratios',        &
        UNITS              = 'kg kg-1',                           &
        FRIENDLYTO         = AEROFRIENDLY,                        &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        DATATYPE           = MAPL_BundleItem,                     &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'AEROTEND',                          &
        LONG_NAME          = 'aerosol_mass_mixing_ratio_tendencies',&
        UNITS              = 'kg kg-1 s-1',                       &
        FRIENDLYTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        DATATYPE           = MAPL_BundleItem,                     &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

! !EXPORT STATE:

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'OX_TEND',                           &
        LONG_NAME          = 'tendency_of_odd_oxygen_mixing_ratio_due_to_chemistry', &
        UNITS              = 'kg kg-1 s-1',                       &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'H2O_TEND',                          &
        LONG_NAME          = 'tendency_of_water_vapor_mixing_ratio_due_to_chemistry', &
        UNITS              = 'kg kg-1 s-1',                       &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'OX_PROD',                           &
        LONG_NAME          = 'tendency_of_odd_oxygen_mixing_ratio_due_to_production', &
        UNITS              = 'pppv s-1',                          &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'OX_LOSS',                           &
        LONG_NAME          = 'tendency_of_odd_oxygen_mixing_ratio_due_to_loss', &
        UNITS              = 'pppv s-1',                          &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'N2O_PROD',                          &
        LONG_NAME          = 'tendency_of_nitrous_oxide_mixing_ratio_due_to_production', &
        UNITS              = 'pppv s-1',                          &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'N2O_LOSS',                          &
        LONG_NAME          = 'tendency_of_nitrous_oxide_mixing_ratio_due_to_loss', &
        UNITS              = 'pppv s-1',                          &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CFC11_PROD',                        &
        LONG_NAME          = 'tendency_of_CFC11_mixing_ratio_due_to_production', &
        UNITS              = 'pppv s-1',                          &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CFC11_LOSS',                        &
        LONG_NAME          = 'tendency_of_CFC11_mixing_ratio_due_to_loss', &
        UNITS              = 'pppv s-1',                          &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CFC12_PROD',                        &
        LONG_NAME          = 'tendency_of_CFC12_mixing_ratio_due_to_production', &
        UNITS              = 'pppv s-1',                          &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CFC12_LOSS',                        &
        LONG_NAME          = 'tendency_of_CFC12_mixing_ratio_due_to_loss', &
        UNITS              = 'pppv s-1',                          &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CFC22_PROD',                        &
        LONG_NAME          = 'tendency_of_CFC22_mixing_ratio_due_to_production', &
        UNITS              = 'pppv s-1',                          &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CFC22_LOSS',                        &
        LONG_NAME          = 'tendency_of_CFC22_mixing_ratio_due_to_loss', &
        UNITS              = 'pppv s-1',                          &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CH4_PROD',                          &
        LONG_NAME          = 'tendency_of_methane_mixing_ratio_due_to_production', &
        UNITS              = 'pppv s-1',                          &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'CH4_LOSS',                          &
        LONG_NAME          = 'tendency_of_methane_mixing_ratio_due_to_loss', &
        UNITS              = 'pppv s-1',                          &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'H2O_PROD',                          &
        LONG_NAME          = 'tendency_of_specific_humidity_due_to_production', &
        UNITS              = 's-1',                               &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'H2O_LOSS',                          &
        LONG_NAME          = 'tendency_of_specific_humidity_due_to_loss', &
        UNITS              = 's-1',                               &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'O3',                                &
        LONG_NAME          = 'ozone_mass_mixing_ratio',           &
        UNITS              = 'kg/kg',                             &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'O3PPMV',                            &
        LONG_NAME          = 'ozone_volume_mixing_ratio',         &
        UNITS              = 'ppmv',                              &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TO3',                               &
        LONG_NAME          = 'total_column_ozone',                &
        UNITS              = 'Dobsons',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TTO3',                              &
        LONG_NAME          = 'tropospheric_column_ozone',         &
        UNITS              = 'Dobsons',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DUST',                              &
        LONG_NAME          = 'mineral_dust_mixing_ratio',         &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'SALT',                              &
        LONG_NAME          = 'sea_salt_mixing_ratio',             &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'SO4',                               &
        LONG_NAME          = 'sulfate_aerosol_mixing_ratio',      &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'BC',                                &
        LONG_NAME          = 'black_carbon_aerosol_mixing_ratio', &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'OC',                                &
        LONG_NAME          = 'organic_carbon_aerosol_mixing_ratio',&
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


!EOP

! Generic Set Services
! --------------------

    call MAPL_GenericSetServices ( GC,RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

end module GEOS_PChemPertGridCompMod
