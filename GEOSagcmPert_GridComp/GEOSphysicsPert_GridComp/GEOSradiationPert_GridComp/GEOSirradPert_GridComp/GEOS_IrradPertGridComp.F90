
!  $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_Irrad -- A Module to compute longwaves radiative transfer 
!          through a cloudy atmosphere

! !INTERFACE:

module GEOS_IrradPertGridCompMod

! !USES:

  use ESMF
  use MAPL_Mod
!  use GEOS_UtilsMod
  
  implicit none
  private

! !PUBLIC ROUTINES:

  public SetServices

!=============================================================================

! !DESCRIPTION:
! 
!   {\tt Irrad} is a light-weight gridded component to compute longwave 
! radiative fluxes. It operates on the ESMF grid that appears in the
! gridded component passed to its {\tt Initialize} method. Unlike
! heavier gridded components, it does not enforce its own grid.
! The only restrictions are that it be a 3-dimensional grid
! in which one dimension is aligned with the vertical coordinate and
! only the horizontal dimensions are decomposed.
!
!   The radiative transfer calculation is based on M-D Chou's IRRAD routine.
! A full documentation of the code may be found in
! "A Thermal Infrared Radiation Parameterization for Atmospheric Studies"
! M.-D. Chou et al., NASA/TM-2001-104606, Vol. 19, 55 pp, 2003.
! Based on the 1996-version of the Air Force Geophysical Laboratory HITRAN data
! base (Rothman et al., 1998), the parameterization includes the absorption due
! to major gaseous absorption (water vapor, CO2 , O3 ) and most of the minor 
! trace gases (N2O, CH4 , CFC's), as well as clouds and aerosols. The thermal
! infrared spectrum is divided into nine bands and a subband. To achieve a high
! degree of accuracy and speed, various approaches of computing the transmission
! function are applied to different spectral bands and gases. The gaseous 
! transmission function is computed either using the k-distribution method or 
! the table look-up method. To include the effect of scattering due to clouds 
! and aerosols, the optical thickness is scaled by the single-scattering albedo
! and asymmetry factor. The optical thickness, the single-scattering albedo, 
! and the asymmetry factor of clouds are parameterized as functions of the ice
! and water content and the particle size.

!   All outputs are optional and are filled only if they have been
! initialized by a coupler. 
!
!   The net (+ve downward) fluxes are returned at the layer
! interfaces, which are indexed from the top of the atmosphere (L=0)
! to the surface. It also computes the sensitivity of net downward flux to 
! surface temperature and emission by the surface.
! The full transfer calculation, including the linearization w.r.t. the surface temperature,
! is done intermitently, on the component's main time step and its results are 
! kept in the internal state. Exports are refreshed each heartbeat based on the
! latest surface temperature.
!
!   Radiation should be called either before or after thos components
!    (usually SURFACE and DYNAMICS) that use its fluxes and modify
!    its inputs. If it is called before, the intemittent refresh should
!    occur during the first step of the radiation cycle, while if it
!    is called after, it should occur during the last step. The behavior
!    of the component needs to be somewhat different in these two cases
!    and so a means is provided, through the logical attribute CALL_LAST in
!    configuration, of telling the component how it is being used. The 
!    default is CALL_LAST = "TRUE". 
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

! Local derived type aliases

    type (ESMF_Config          )            :: CF

    integer      :: MY_STEP
    integer      :: ACCUMINT
    real         :: DT

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

! Get the configuration
! ---------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

! Get the intervals; "heartbeat" must exist
! -----------------------------------------

    call ESMF_ConfigGetAttribute( CF, DT, Label="RUN_DT:"                          , RC=STATUS)
    VERIFY_(STATUS)

! Refresh interval defaults to heartbeat. This will also be read by
!  MAPL_Generic and set as the component's main time step.
! -----------------------------------------------------------------

    call ESMF_ConfigGetAttribute( CF, DT, Label=trim(COMP_NAME)//"_DT:", default=DT, RC=STATUS)
    VERIFY_(STATUS)

    MY_STEP = nint(DT)

! Averaging interval defaults to the refresh interval.
!-----------------------------------------------------

    call ESMF_ConfigGetAttribute(CF, DT, Label=trim(COMP_NAME)//'Avrg:', default=DT, RC=STATUS)
    VERIFY_(STATUS)

    ACCUMINT = nint(DT)

! Set the state variable specs.
! -----------------------------

!BOP

!  !IMPORT STATE:

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'PLE',                               &
        LONG_NAME          = 'air_pressure',                      &
        UNITS              = 'Pa',                                &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationEdge,                  &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'T',                                 &
        LONG_NAME          = 'air_temperature',                   &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'QV',                                &
        LONG_NAME          = 'specific_humidity',                 &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'QL',                                &
        LONG_NAME          = 'mass_fraction_of_cloud_liquid_water_in_air', &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'QI',                                &
        LONG_NAME          = 'mass_fraction_of_cloud_ice_in_air', &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'QR',                                &
        LONG_NAME          = 'mass_fraction_of_rain_water_in_air',&
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'RL',                                &
        LONG_NAME          = 'effective_radius_of_cloud_liquid_water_particles',      &
        UNITS              = 'm',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'RI',                                &
        LONG_NAME          = 'effective_radius_of_cloud_ice_particles',   &
        UNITS              = 'm',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'RR',                                &
        LONG_NAME          = 'effective_radius_of_rain_particles',&
        UNITS              = 'm',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'O3',                                &
        LONG_NAME          = 'ozone_mixing_ratio',                &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CH4',                               &
        LONG_NAME          = 'methane_concentration',             &
        UNITS              = 'pppv',                              &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'N2O',                               &
        LONG_NAME          = 'nitrous_oxide_concentration',       &
        UNITS              = 'pppv',                              &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CFC11',                             &
        LONG_NAME          = 'CFC11_concentration',               &
        UNITS              = 'pppv',                              &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CFC12',                             &
        LONG_NAME          = 'CFC12_concentration',               &
        UNITS              = 'pppv',                              &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'CFC22',                             &
        LONG_NAME          = 'CFC22_concentration',               &
        UNITS              = 'pppv',                              &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'FCLD',                              &
        LONG_NAME          = 'cloud_area_fraction_in_atmosphere_layer', &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'TS',                                &
        LONG_NAME          = 'surface_skin_temperature',          &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'EMIS',                              &
        LONG_NAME          = 'surface_emissivity',                &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        AVERAGING_INTERVAL = ACCUMINT,                            &
        REFRESH_INTERVAL   = MY_STEP,                             &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'PREF',                              &
        LONG_NAME          = 'reference_air_pressure',            &
        UNITS              = 'Pa',                                &
        DIMS               = MAPL_DimsVertOnly,                   &
        VLOCATION          = MAPL_VLocationEdge,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

! Instantaneous TS is used only for updating the IR fluxes due to TS change

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'TSINST',                            &
        LONG_NAME          = 'surface_skin_temperature',          &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'aerosols',                                              &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'AERO',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       DATATYPE   = MAPL_BundleItem,                                         &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

!  !EXPORT STATE:

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'FLX',                                       &
        LONG_NAME  = 'net_downward_longwave_flux_in_air',         &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'FLC',                                       &
        LONG_NAME  = 'net_downward_longwave_flux_in_air_assuming_clear_sky',  &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'FLA',                                       &
        LONG_NAME  = 'net_downward_longwave_flux_in_air_assuming_clear_sky_and_no_aerosol',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'SFCEM',                                     &
        LONG_NAME  = 'longwave_flux_emitted_from_surface',        &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'DSFDTS',                                    &
        LONG_NAME  = 'sensitivity_of_longwave_flux_emitted_from_surface_to_surface_temperature', &
        UNITS      = 'W m-2 K-1',                                 &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'TSREFF',                                    &
        LONG_NAME  = 'surface_temperature',                       &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'OLR',                                       &
        LONG_NAME  = 'upwelling_longwave_flux_at_toa',            &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'OLC',                                       &
        LONG_NAME  = 'upwelling_longwave_flux_at_toa_assuming_clear_sky',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'OLA',                                       &
        LONG_NAME  = 'upwelling_longwave_flux_at_toa_assuming_clear_sky_and_no_aerosol',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'FLNS',                                      &
        LONG_NAME  = 'surface_net_downward_longwave_flux',        &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'FLNSC',                                     &
        LONG_NAME  = 'surface_net_downward_longwave_flux_assuming_clear_sky',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'FLNSA',                                     &
        LONG_NAME  = 'surface_net_downward_longwave_flux_assuming_clear_sky_and_no_aerosol',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'LWS',                                       &
        LONG_NAME  = 'surface_downwelling_longwave_flux',         &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'LCS',                                       &
        LONG_NAME  = 'surface_downwelling_longwave_flux_assuming_clear_sky',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'LAS',                                       &
        LONG_NAME  = 'surface_downwelling_longwave_flux_assuming_clear_sky_and_no_aerosol',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'CLDTMP',                                    &
        LONG_NAME  = 'cloud_top_temperature',                     &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'CLDPRS',                                    &
        LONG_NAME  = 'cloud_top_pressure',                        &
        UNITS      = 'Pa',                                        &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
        SHORT_NAME = 'TAUIR',                                     &
        LONG_NAME  = 'longwave_cloud_optical_thickness_at_800_cm-1',&
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)


!  Irrad does not have a "real" internal state. To update the net_longwave_flux
!  due to the change of surface temperature every time step, we keep 
!  several variables in the internal state.


!  !INTERNAL STATE:

    call MAPL_AddInternalSpec(GC,                            &
        SHORT_NAME = 'FLX',                                       &
        LONG_NAME  = 'net_downward_longwave_flux_in_air',         &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                            &
        SHORT_NAME = 'FLC',                                       &
        LONG_NAME  = 'net_downward_longwave_flux_in_air_for_clear_sky(INTERNAL)',  &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                            &
        SHORT_NAME = 'FLA',                                       &
        LONG_NAME  = 'net_downward_longwave_flux_in_air_for_clear_sky_and_no_aerosol',  &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                            &
        SHORT_NAME = 'DFDTS',                                     &
        LONG_NAME  = 'sensitivity_of_net_downward_longwave_flux_in_air_to_surface_temperature', &
        UNITS      = 'W m-2 K-1',                                 &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                            &
        SHORT_NAME = 'DFDTSC',                                    &
        LONG_NAME  = 'sensitivity_of_net_downward_longwave_flux_in_air_to_surface_temperature_for_clear_sky',&
        UNITS      = 'W m-2 K-1',                                 &
        DIMS       = MAPL_DimsHorzVert,                           &
        VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                            &
        SHORT_NAME = 'SFCEM',                                     &
        LONG_NAME  = 'longwave_flux_emitted_from_surface',        &
        UNITS      = 'W m-2',                                     &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                            &
        SHORT_NAME = 'TS',                                        &
        LONG_NAME  = 'surface_temperature',                       &
        UNITS      = 'K',                                         &
        DIMS       = MAPL_DimsHorzOnly,                           &
        VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

!EOP

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,    name="-LW_DRIVER"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--IRRAD"      ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="--MISC"       ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-UPDATE_FLX"  ,RC=STATUS)
    VERIFY_(STATUS)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

end module GEOS_IrradPertGridCompMod

