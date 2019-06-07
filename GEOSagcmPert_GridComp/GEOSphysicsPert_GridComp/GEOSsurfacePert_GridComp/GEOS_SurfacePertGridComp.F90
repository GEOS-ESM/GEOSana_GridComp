
!  $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_Surface   -- A composite component for the surface components.

! !INTERFACE:

module GEOS_SurfacePertGridCompMod

! !USES:

  use ESMF
  use MAPL_Mod

!  use GEOS_LakeGridCompMod,      only : LakeSetServices     => SetServices
!  use GEOS_LandiceGridCompMod,   only : LandiceSetServices  => SetServices
!  use GEOS_SaltwaterGridCompMod, only : OceanSetServices    => SetServices
!  use GEOS_LandGridCompMod,      only : LandSetServices     => SetServices

  implicit none
  private

  type( ESMF_VM ) :: VMG

! !PUBLIC ROUTINES:

  public SetServices

!=============================================================================

! !DESCRIPTION:
! 
!   {\tt GEOS\_Surface} is a light-weight gridded component that implements the
!      interface to the tiled surface components. The surface computational components
!      (LAND, LAKE, OCEAN, LANDICE) are its children. All of {\tt GEOS\_Surface}'s imports and exports
!      are in the atmospheric model's grid. In {\tt GEOS\_Surface} these are transformed to the
!      exchange grid, and the relevant portions of the exchange grid are passed to
!      each of the children. The children's results are them replaced in
!      the full exchange grid and transformed back to the atmospheric grid.
!
!      {\tt GEOS\_Surface} has two run stages, as do its children. These are meant
!      to interface with the two stages of {\tt GEOS\_Turbulence}. During the first run
!      stage, the children all produce surface exchange coefficients, and during the
!      second, they update the surface state and produce final values of the fluxes.
!
!      {\tt GEOS\_Surface} keeps a Private Internal State called 'SURF_state' in the 
!      component object. In this state it saves the tranforms between the atmospheric
!      grid and each of the children's exchange grids. This should be done more 
!      elegantly once ESMF has exchange grid support. It also has a \gg Internal State
!      that is used to communicate between the two run methods. These internal states
!      do not need to be saved in restarts.
!
!      The four children of {\tt GEOS\_Surface} are given the names:
!      'LAKE', which treats inland freshwater bodies; 'LANDICE', which treats permanent
!      glaciers; 'LAND', which treats all other land surface types, both bare and vegetated,
!      as well as vegetated wetlands not considered freshwater bodies; and  'SALTWATER', which
!      performs the surface calculations for all ocean areas. All four operate in lists
!      of tiles that are nonoverlapping subsets of the exhange grid, and their union---the
!      full exchange grid---tiles the entire sphere.
!
!      By default MAPL_Generic tries to resolve Imports and Exports among
!      the children; but the children of {\tt GEOS\_Surface} do not talk directly to each other,
!      and all communication between them would need to be performed by {\tt GEOS\_Surface} manipulating
!      their Import and Export states. Currently they do not communicate, but this will
!      chage when runoff routing is implemented.

!EOP


     
!  integer ::        LAKE
!  integer ::     LANDICE
!  integer ::       OCEAN 
!  integer ::        LAND 
!  integer, parameter :: NUM_CHILDREN = 4

!  character(len=ESMF_MAXSTR), pointer :: GCNames(:)
!  integer                    :: CHILD_MASK(NUM_CHILDREN)


   contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer,             intent(  OUT) :: RC  ! return code

! !DESCRIPTION: This version uses the GEOS\_GenericSetServices, which in addition
!                to setting default IRF methods, also allocates
!   our instance of a generic state and puts it in the 
!   gridded component (GC). Here we override the Initialize and Run methods.
!   The Run method is a two-stage method that implemets the interaction
!   between the 2-stage children representing the various surface types and the 2-stage
!   turbulence run methods.
!
!   
!   Note that, in addition to its explicit exports,
!   the entire internal state, which is used to communicate between the two run stages, 
!   is exported using the ``friendly-to-self'' mechanism.

!EOP

!=============================================================================
!
! ErrLog Variables


    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Locals

    integer                                 :: I
    type (MAPL_MetaComp    ), pointer       :: MAPL

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

!BOP

! Imports are read-only quantities computed by other gridded components.

! Note that the turbulence fluxes appearing in the import state are
! the values computed by the first run stage of turbulence using fixed
! surface conditions. The Export versions of these fluxes are the final
! values actually used in the surface budgets. The same applies to some
! of the radiative fluxes, for which the values exported here are those
! actually used in the budget. 


!  !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_pressure',                  &
        UNITS              = 'Pa',                                &
        SHORT_NAME         = 'PS',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_air_temperature',           &
        UNITS              = 'K',                                 &
        SHORT_NAME         = 'TA',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_air_specific_humidity',     &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'QA',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_wind_speed',                &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'SPEED',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'UA',                                        &
        LONG_NAME  = 'eastward_wind_bottom_level',                &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                            &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'VA',                                        &
        LONG_NAME  = 'northward_wind_bottom_level',               &
        UNITS      = 'm s-1',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,                          &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_layer_height',              &
        UNITS              = 'm',                                 &
        SHORT_NAME         = 'DZ',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface geopotential height',       &
        UNITS              = 'm2 sec-2',                          &
        SHORT_NAME         = 'PHIS',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'sensible_heat_flux',                &
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'SH',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'eastward_surface_stress_on_air',    &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUX',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'northward_surface_stress_on_air',   &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUY',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'evaporation',                       &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'EVAP',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'dewfall',                           &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'DEWL',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'frostfall',                         &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'FRSL',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_sensible_heat_wrt_dry_static_energy',&
        UNITS              = 'W m-2 K-1',                         &
        SHORT_NAME         = 'DSH',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_eastward_surface_stress_wrt_Us', &
        UNITS              = 'N s m-3',                           &
        SHORT_NAME         = 'DFU',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_northward_surface_stress_wrt_Us', &
        UNITS              = 'N s m-3',                           &
        SHORT_NAME         = 'DFV',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_evaporation_wrt_QS',  &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'DEVAP',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_dewfall_wrt_QS',  &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'DDEWL',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'derivative_of_frostfall_wrt_QS',  &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'DFRSL',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'liquid_water_convective_precipitation', &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'PCU',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'liquid_water_large_scale_precipitation', &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'PLS',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'snowfall',                          &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'SNO',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
        SHORT_NAME         = 'DRPARN',                            &
        LONG_NAME          = 'normalized_surface_downwelling_par_beam_flux', &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
        SHORT_NAME         = 'DFPARN',                            &
        LONG_NAME          = 'normalized_surface_downwelling_par_diffuse_flux', &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         SHORT_NAME         = 'DRNIRN'                      ,&
         LONG_NAME          = 'normalized_surface_downwelling_nir_beam_flux',&
         UNITS              = '1'                           ,&
         DIMS               = MAPL_DimsHorzOnly,             &
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         SHORT_NAME         = 'DFNIRN'                      ,&
         LONG_NAME          = 'normalized_surface_downwelling_nir_diffuse_flux',&
         UNITS              = '1'                           ,&
         DIMS               = MAPL_DimsHorzOnly,             &
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         SHORT_NAME         = 'DRUVRN'                      ,&
         LONG_NAME          = 'normalized_surface_downwelling_uvr_beam_flux',&
         UNITS              = '1'                           ,&
         DIMS               = MAPL_DimsHorzOnly,             &
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         SHORT_NAME         = 'DFUVRN'                      ,&
         LONG_NAME          = 'normalized_surface_downwelling_uvr_diffuse_flux',&
         UNITS              = '1'                           ,&
         DIMS               = MAPL_DimsHorzOnly,             &
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
        SHORT_NAME         = 'LWDNSRF',                           &
        LONG_NAME          = 'surface_downwelling_longwave_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
        SHORT_NAME         = 'ALW',                               &
        LONG_NAME          = 'linearization_of_surface_upwelling_longwave_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
        SHORT_NAME         = 'BLW',                               &
        LONG_NAME          = 'linearization_of_surface_upwelling_longwave_flux', &
        UNITS              = 'W m-2 K-1',                         &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

!  !EXPORT STATE:

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_visible_beam',   &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBVR',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_visible_diffuse',&
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBVF',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_nearinfrared_beam', &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBNR',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'surface_albedo_for_nearinfraed_diffuse', &
        UNITS              = '1',                                 &
        SHORT_NAME         = 'ALBNF',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'EMIS',                              &
        LONG_NAME          = 'surface_emissivity',                &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'Z0',                                &
        LONG_NAME          = 'surface_roughness',                 &
        UNITS              = 'm',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'MOU10M',                            &
        LONG_NAME          = 'zonal 10m wind from MO sfc',        &
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'MOV10M',                            &
        LONG_NAME          = 'meridional 10m wind from MO sfc',   &
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'MOT2M',                             &
        LONG_NAME          = 'temperature 2m wind from MO sfc',   &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'MOQ2M',                             &
        LONG_NAME          = 'humidity 2m wind from MO sfc',      &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'Z0H',                               &
        LONG_NAME          = 'surface_roughness_for_heat',        &
        UNITS              = 'm',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'RI',                                &
        LONG_NAME          = 'surface_bulk_richardson_number',    &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'RE',                                &
        LONG_NAME          = 'surface_reynolds_number',           &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'FRACI',                             &
        LONG_NAME          = 'ice_covered_fraction_of_tile',      &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'QDWL',                              &
        LONG_NAME          = 'surface_liquid_condensate',         &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'QFRL',                              &
        LONG_NAME          = 'surface_ice_condensate',            &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'SHAT',                              &
        LONG_NAME          = 'effective_surface_dry_static_energy',&
        UNITS              = 'm+2 s-2',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELUS',                             &
        LONG_NAME          = 'change_of_surface_eastward_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELVS',                             &
        LONG_NAME          = 'change_of_surface_northward_velocity',&
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELSS',                             &
        LONG_NAME          = 'change_of_surface_dry_static_energy',&
        UNITS              = 'm+2 s-2',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)
     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELTS',                             &
        LONG_NAME          = 'change_of_surface_skin_temperature',&
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DELQS',                             &
        LONG_NAME          = 'change_of_surface_specific_humidity',&
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DLQLL',                             &
        LONG_NAME          = 'change_of_surface_liquid_condensate',&
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'DLQIL',                             &
        LONG_NAME          = 'change_of_surface_frozen_condensate',&
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'FRLAND',                            &
        LONG_NAME          = 'fraction_of_land',                  &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'FRLANDICE',                         &
        LONG_NAME          = 'fraction_of_land_ice',              &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'FRLAKE',                            &
        LONG_NAME          = 'fraction_of_lake',                  &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'FROCEAN',                           &
        LONG_NAME          = 'fraction_of_ocean',                 &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'USTAR',                             &
        LONG_NAME          = 'surface_velocity_scale',            &
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TSTAR',                             &
        LONG_NAME          = 'surface_temperature_scale',         &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'QSTAR',                             &
        LONG_NAME          = 'surface_moisture_scale',            &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'BSTAR',                             &
        LONG_NAME          = 'surface_bouyancy_scale',            &
        UNITS              = 'm s-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'TSOIL1',                            &
        LONG_NAME          = 'soil_temperatures_layer_1'         ,&
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_land_snowcover',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ASNOW'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_temperature_of_snow',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPSNOW'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_temperature_of_saturated_zone',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPSAT'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_temperature_of_unsaturated_zone',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPUNST'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_temperature_of_wilted_zone',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPWLT'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_saturated_zone',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'FRSAT'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_unsaturated_zone',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'FRUST'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_wilting_zone',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'FRWLT'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'snow_water_equivalent_depth'       ,&
        UNITS              = 'mm'                                ,&
        SHORT_NAME         = 'SNOMAS'                            ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'surface_soil_wetness'              ,&
        UNITS              = '1'                                 ,&
        SHORT_NAME         = 'WET1'                              ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'root_zone_soil_wetness'            ,&
        UNITS              = '1'                                 ,&
        SHORT_NAME         = 'WET2'                              ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'leaf_area_index'                   ,&
        UNITS              = '1'                                 ,&
        SHORT_NAME         = 'LAI'                               ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'greeness_fraction'                 ,&
        UNITS              = '1'                                 ,&
        SHORT_NAME         = 'GRN'                               ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'canopy_height'                     ,&
        UNITS              = 'm'                                 ,&
        SHORT_NAME         = 'Z2CH'                              ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'root_length'                       ,&
        UNITS              = 'mm'                                ,&
        SHORT_NAME         = 'ROOTL'                             ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'sensible_heat_flux_from_turbulence',&
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'SH',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'eastward_surface_stress',           &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUX',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'northward_surface_stress',          &
        UNITS              = 'N m-2',                             &
        SHORT_NAME         = 'TAUY',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'evaporation_from_turbulence',       &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'EVAP',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '10-meter_eastward_wind',                                &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'U10M',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '10-meter_northward_wind',                               &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'V10M',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'equivalent_neutral_10-meter_eastward_wind',             &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'U10N',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'equivalent_neutral_10-meter_northward_wind',            &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'V10N',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '50-meter_eastward_wind',                                &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'U50M',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '50-meter_northward_wind',                               &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'V50M',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '10-meter_air_temperature',                              &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'T10M',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '10-meter_specific_humidity',                            &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'Q10M',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '2-meter_eastward_wind',                                 &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'U2M',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '2-meter_northward_wind',                                &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'V2M',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '2-meter_air_temperature',                               &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'T2M',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = '2-meter_specific_humidity',                             &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'Q2M',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_air_temperature',                               &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'TA',                                                    &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_air_specific_humidity',                         &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'QA',                                                    &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_eastward_wind',                                 &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'UA',                                                    &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_northward_wind',                                &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'VA',                                                    &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       SHORT_NAME = 'GUST',                                                  &
       LONG_NAME  = 'gustiness',                                             &
       UNITS      = 'm s-1',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       SHORT_NAME = 'VENT',                                                  &
       LONG_NAME  = 'surface_ventilation_velocity',                          &
       UNITS      = 'm s-1',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
       LONG_NAME          = 'land_water_ice_flag',               &
       UNITS              = '0-1-2',                             &
       SHORT_NAME         = 'LWI',                               &
       DIMS               = MAPL_DimsHorzOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
                                                      RC=STATUS  )
    VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                             ,&
        LONG_NAME          = 'snow_depth'                        ,&
        UNITS              = 'm'                                 ,&
        SHORT_NAME         = 'SNOWDP'                            ,&
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'eastward_stress_over_water',&
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUXW'                     ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'northward_stress_over_water',&
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUYW'                     ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'eastward_stress_over_ice',  &
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUXI'                     ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        LONG_NAME          = 'northward_stress_over_ice',  &
        UNITS              = 'N m-2'                     ,&
        SHORT_NAME         = 'TAUYI'                     ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'open_water_upward_sensible_heat_flux' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SHWTR'                     ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'sea_ice_upward_sensible_heat_flux' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SHICE'                     ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'open_water_latent_energy_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'HLATWTR'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'sea_ice_latent_energy_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'HLATICE'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)
     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'open_water_net_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'LWNDWTR'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'sea_ice_net_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'LWNDICE'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'open_water_net_downward_shortwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SWNDWTR'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                     ,&
        LONG_NAME          = 'sea_ice_net_downward_shortwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SWNDICE'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'ocean_snowfall'            ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'SNOWOCN'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
        LONG_NAME          = 'ocean_rainfall'            ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'RAINOCN'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)



  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'evaporation'               ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'EVAPOUT'                   ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'upward_sensible_heat_flux' ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'SHOUT'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'runoff_flux'               ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'RUNOFF'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'interception_loss_energy_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPINT'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'baresoil_evap_energy_flux' ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPSOI'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'transpiration_energy_flux' ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPVEG'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snowpack_evaporation_energy_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPICE'                    ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'baseflow_flux'             ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'BASEFLOW'                  ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_runoff_flux'       ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'RUNSURF'                   ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'EVLAND',                    &
    LONG_NAME          = 'Evaporation_land',          &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LHLAND',                    &
    LONG_NAME          = 'Latent_heat_flux_land',     &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SHLAND',                    &
    LONG_NAME          = 'Sensible_heat_flux_land',   &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SWLAND',                    &
    LONG_NAME          = 'Net_shortwave_land',        &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LWLAND',                    &
    LONG_NAME          = 'Net_longwave_land',         &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GHLAND',                    &
    LONG_NAME          = 'Ground_heating_land',       &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SMLAND',                    &
    LONG_NAME          = 'Snowmelt_flux_land',        &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TWLAND',                    &
    LONG_NAME          = 'Avail_water_storage_land',  &
    UNITS              = 'kg m-2',                    &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TELAND',                    &
    LONG_NAME          = 'Total_energy_storage_land', &
    UNITS              = 'J m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TSLAND',                    &
    LONG_NAME          = 'Total_snow_storage_land',   &
    UNITS              = 'kg m-2',                    &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'DWLAND',                    &
    LONG_NAME          = 'rate_of_change_of_total_land_water',&
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'DHLAND',                    &
    LONG_NAME          = 'rate_of_change_of_total_land_energy',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPLAND',                    &
    LONG_NAME          = 'rate_of_spurious_land_energy_source',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPWATR',                    &
    LONG_NAME          = 'rate_of_spurious_land_water_source',&
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPSNOW',                    &
    LONG_NAME          = 'rate_of_spurious_snow_energy',&
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snowmelt_flux'             ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'SMELT'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_outgoing_longwave_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'HLWUP'                     ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                    ,&
    LONG_NAME          = 'surface_net_downward_longwave_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'LWNDSRF'                   ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
    VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                    ,&
    LONG_NAME          = 'surface_net_downward_shortwave_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'SWNDSRF'                   ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
    VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'total_latent_energy_flux'  ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'LHFX'                      ,&
    DIMS               = MAPL_DimsHorzOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'ACCUM',                     &
    LONG_NAME          = 'net_ice_accumulation_rate', &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                         &
    SHORT_NAME         = 'ITY',                       &
    LONG_NAME          = 'vegetation_type',           &
    UNITS              = '1',                         &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                         &
    SHORT_NAME         = 'NITY',                      &
    LONG_NAME          = 'NCEP_vegetation_type',      &
    UNITS              = '1',                         &
    DIMS               = MAPL_DimsHorzOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

! !INTERNAL STATE:

!  These are here only because they are passed between run1 and run2.
!  They don't need to be saved in restarts. Note they are all exported
!  by being made friendly to self.
!  Some may be needed by turbulence, but not in a Friendly way; others
!  are only diagnostics.


     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'QS',                                &
        LONG_NAME          = 'surface_specific_humidity',         &
        UNITS              = '1',                                 &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'TS',                                &
        LONG_NAME          = 'surface_skin_temperature',          &
        UNITS              = 'K',                                 &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CT',                                &
        LONG_NAME          = 'surface_exchange_coefficient_for_heat', &
        UNITS              = 'kg m-2 s-1',                        &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CQ',                                &
        LONG_NAME          = 'surface_exchange_coefficient_for_moisture', &
        UNITS              = 'kg m-2 s-1',                        &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CM',                                &
        LONG_NAME          = 'surface_exchange_coefficient_for_momentum', &
        UNITS              = 'kg m-2 s-1',                        &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'CN',                                &
        LONG_NAME          = 'surface_neutral_drag_coefficient',  &
        UNITS              = '1',                                 &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'THAT',                              &
        LONG_NAME          = 'effective_surface_skin_temperature',&
        UNITS              = 'K',                                 &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'QHAT',                              &
        LONG_NAME          = 'effective_surface_specific_humidity',&
        UNITS              = '1',                                 &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'UHAT',                              &
        LONG_NAME          = 'effective_surface_eastward_velocity',&
        UNITS              = 'm s-1',                             &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        SHORT_NAME         = 'VHAT',                              &
        LONG_NAME          = 'effective_surface_northward_velocity',&
        UNITS              = 'm s-1',                             &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC,                           &
        LONG_NAME          = 'air_density_at_surface',            &
        UNITS              = 'kg m-3',                            &
        SHORT_NAME         = 'RHOS',                              &
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddInternalSpec(GC                           ,&
        LONG_NAME          = 'zero_plane_displacement_height'    ,&
        UNITS              = 'm'                                 ,&
        SHORT_NAME         = 'D0'                                ,&
        FriendlyTO         = trim(COMP_NAME),                     &
        DIMS               = MAPL_DimsHorzOnly                   ,&
        VLOCATION          = MAPL_VLocationNone                  ,&
                                                       RC=STATUS  ) 

     VERIFY_(STATUS)

!EOP


! Call SetServices for children
!------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS )
    VERIFY_(STATUS)
 
    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_SurfacePertGridCompMod


