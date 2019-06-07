!  $Id$

#include "MAPL_Generic.h"

module GEOS_SolarPertGridCompMod

!=============================================================================
!BOP

! !MODULE: GEOS_SolarGridCompMod -- Computes solar radiation fluxes in a cloudy atmosphere


! !DESCRIPTION:
! 
! {\tt GEOS\_SolarGridCompMod} is an ESMF/MAPL gridded component that performs
!  a broadband calculation of shortwave radiative fluxes for use as a
!  solar radiation parameterization in atmospheric models on a sphere. \newline
! 
! {\em Scientific Basis:} The radiative transfer calculation is based on M-D Chou shortwave
! parameterization. The basic reference for the scheme is:
! Chou and Suarez 1999: A Solar Radiation Parameterization for Atmospheric Studies,
! NASA-TM-1999-104606, Vol 15. An updated version of this report can be
! found in SolarDoc.pdf in this directory. \newline
!
! The parameterization treats direct and diffuse fluxes of solar 
! radiation in eight spectral bands:
! \begin{verbatim}
!        in the uv region :
!           index  1 for the 0.225-0.285 micron band
!           index  2 for the 0.175-0.225;0.285-0.300 micron band
!           index  3 for the 0.300-0.325 micron band
!           index  4 for the 0.325-0.4 micron band
!        in the par region :
!           index  5 for the 0.4-0.690 micron band
!        in the infrared region :
!           index  6 for the 0.690-1.220 micron band
!           index  7 for the 1.220-2.270 micron band
!           index  8 for the 2.270-3.850 micron band
! \end{verbatim}
! It includes gaseous absorption due to water vapor, ozone, CO$_2$, and
! molecular oxygen and the effects of molecular scattering,
! as well as multiple scattering due to clouds and aerosols. 
!
! It allows clouds to occur in any layer and 
! horizontal cloud cover fractions must be specified 
! for all layers; clear layers
! simply have a fraction of zero. Vertically, the layers are 
! assumed to be filled by cloud. To simplify the treatment of
! cloud effects, the model layers,
! are grouped into three super layers. Effective cloud properties are then
! parameterized by assuming that clouds are maximally overlapped within the super layers
! and randomly overlapped between the super layers.  The  
! optical properties of cloud particles depend on the liquid, ice, and rain mixing ratios,
! as well as on spatially dependent effective radii for the three species.
! These are all inputs to the component. \newline
!
!  The parameterization can include the effects of an arbitrary
!  number of aerosol species.
!  Aerosol optical thickness, single-scattering albedo, and asymmetry
!  factor must be determined as functions of height and spectral band
!  for each species. \newline
!
!
! {\em Code Implementation:} \newline
! 
!  {\tt GEOS\_SolarGridCompMod} is an encapsulation of Chou's plug-compatible
!  SORAD Fortran routine in a MAPL/ESMF gridded component (GC). 
!  It follows the standard rules for an ESMF/MAPL GCs.
!  It operates on the ESMF grid that appears in the
!  gridded component. This grid must
!  be present in the GC and properly initialized before Initialize
!  is called. The only restrictions on the grid are that it be 3-dimensional
!  with two horizontal and one vertical dimension and
!  with only the horizontal dimensions decomposed. The vertical dimension
!  is also assumed to the the thrid dimension of the Fortran arrays and
!  is indexed from the top down. No particular vertical coordinate is assumed,
!  rather the 3-dimensional field of air pressure at the layer interfaces is 
!  a required Import. \newline
!
!  This module contains only SetServices and Run methods. 
!  The Initialize and Finalize methods
!  being defaulted to the MAPL\_Generic versions. 
!  The SetServices method is the only public
!  entity. There are no public types or data. \newline
!
!  The contents of the Import, Export, and Internal States are explicitly
!  described in SetServices and in tables in this documentation.
!  All quantities in these states are in either ESMF Fields or Bundles,
!  and all share a common grid---the ESMF grid in the gridded component
!  at the time Initialize (MAPL\_GenericInitialize, in this case) was called.
!  All outputs appearing in the Export state are optional and are 
!  filled only if they have been allocated. All filled Exports are valid
!  for the time interval on the GC's clock when the run method is invoked.
!  Imports can be from either an instantaneous or a time-averaged state of the
!  atmosphere. All Imports are read-only; none are Friendly.
!  Most imports are simple ESMF Fields containing 2- or 
!  3-dimensional quantities, such as temperature and humidity, needed in
!  the flux calculation. Non-cloud aerosol amounts are the exception; they 
!  appear in an ESMF Bundle.  \newline 
!
!  The net (+ve downward) fluxes on the Export state are defined at the layer
!  interfaces, which are indexed from the top of the atmosphere (L=0)
!  to the surface. Incident fluxes
!  at the surface also appear in the Export state; these are separated
!  into direct (beam) and diffuse fluxes for three spectral bands
!  (uv, par, nir), as defined in the table above.  \newline
!
!  The full transfer calculation is done infrequently and 
!  its results kept in the Internal state. 
!  The frequency of full calculations is controlled
!  by an alarm whose interval can be set
!  from a value in the configuration and whose origin is taken as the 
!  beginning of the run.
!  For the full calculations, solar fluxes are computed based on
!  mean zenith angles averaged over sun positions for a 
!  given period (the long interval, which can be specified in 
!  the configuration) beyond the
!  current time on the input clock. On every call to the Run method,
!  whatever the state of the alarm that controls the full calculation,
!  the sun's position
!  is updated to the mean position for the clock's current interval
!  and fluxes are updated based on normalized fluxes computed during
!  the previous full transfer calculation, but using
!  the TOA insolation for the current time on the clock. Because of this
!  intermittent scheme, checkpoint-restart sequences are seamless
!  only when interrupted at the time of the full calculation.\newline
!
!  The calculation relies in MAPL's Astronomy layer, which in turn 
!  assumes that the ESMF grid can be queried for latitude and longitude
!  coordinates. \newline
!
! {\em Configuration:} \newline
!  
!  Like all MAPL GCs, {\tt GEOS\_SolarGridCompMod} assumes that the configuration
!  in the ESMF GC is open and treats it as an environment from which it can
!  {\em at any time} read control information. It uses MAPL rules for scanning this
!  configuration.
!\begin{verbatim}
!
!VARIABLE             DESCRIPTION           UNITS      DEFAULT   NOTES
!
!RUN_DT:              Short time interval   (seconds)  none   
!DT:                  Long time interval    (seconds)  RUN_DT
!AVGR:                Averaging interval    (seconds)  DT
!PRS_LOW_MID_CLOUDS:  Interface pressure    (Pa)       70000.
!                     between the low and
!                     middle cloud layers
!PRS_MID_HIGH_CLOUDS: Interface pressure    (Pa)       40000.
!                     between the high and
!                     middle cloud layers
!SOLAR_CONSTANT:                            (W m-2)    1365.
!CO2:                 CO2 concentration     (ppmv)     350.e-6   If negative,
!                                                                use historical CO2
!
!\end{verbatim}
!
!
! !BUGS:
!
!\end{verbatim} 
!\begin{itemize} 
!    \item Aerosol properties for each aerosol in the Bundle are obtained by
!    calling a global method (Get\_AeroOptProp) that must recognize
!    the aerosol by its Field name in the Bundle. This is a placeholder
!    for a scheme in which each Field carries with it a method for computing
!    its aerosol's optical properties.
!
!    \item The grid must have two horizontal dimensions and they must be the inner dimensions
!    of Fortran arrays.
!
!    \item The load-balancing relies on the grid describing a sphere. Everything
!    works for non-spherical grids but the load-balancing should be disabled
!    and this can be done only by going into the code.
!\end{itemize} 
! \begin{verbatim}
!
!     

! !USES:

  use ESMF
  use MAPL_Mod
!  use AeroOptPropTableMod

  implicit none
  private

! !PUBLIC ROUTINES:

  public SetServices


!EOP

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a MAPL\_MetaComp and putting it in the 
!   gridded component (GC). Here we only need to register the Run method with ESMF and
!   register the state variable specifications with MAPL. \newline
!   

!EOP

!=============================================================================

! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type (ESMF_Config          )            :: CF

! Locals

    integer      :: RUN_DT
    integer      :: MY_STEP
    integer      :: ACCUMINT
    real         :: DT
    integer      :: I

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // 'SetServices'

! Get the configuration
! ---------------------

    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

! Get the intervals; "heartbeat" must exist
! -----------------------------------------

    call ESMF_ConfigGetAttribute(CF, DT, Label="RUN_DT:"                          , RC=STATUS)
    VERIFY_(STATUS)

    RUN_DT = nint(DT)

! Refresh interval defaults to heartbeat.
! ---------------------------------------

    call ESMF_ConfigGetAttribute(CF, DT, Label=trim(COMP_NAME)//"_DT:", default=DT, RC=STATUS)
    VERIFY_(STATUS)

    MY_STEP = nint(DT)

! Averaging interval defaults to refresh interval.
!-------------------------------------------------

    call ESMF_ConfigGetAttribute(CF, DT, Label=trim(COMP_NAME)//'Avrg:',default=DT,RC=STATUS)
    VERIFY_(STATUS)

    ACCUMINT = nint(DT)


! Set the state variable specs.
! -----------------------------

!BOP

! !IMPORT STATE:

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'air_pressure',                                          &
       UNITS      = 'Pa',                                                    &
       SHORT_NAME = 'PLE',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                         &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'air_temperature',                                       &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'T',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                         &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'specific_humidity',                                     &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'QV',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                         &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'mass_fraction_of_cloud_liquid_water_in_air',            &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'QL',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                         &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'mass_fraction_of_cloud_ice_in_air',                     &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'QI',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                         &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'mass_fraction_of_rain_water_in_air',                    &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'QR',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                         &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'effective_radius_of_cloud_liquid_water_particles',      &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'RL',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                         &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'effective_radius_of_cloud_ice_particles',               &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'RI',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                         &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'effective_radius_of_rain_particles',                    &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'RR',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                         &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'odd_oxygen_mixing_ratio',                               &
       UNITS      = 'kg/kg',                                                 &
       SHORT_NAME = 'OX',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                         &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'ozone_mixing_ratio',                                    &
       UNITS      = 'mol/mol',                                               &
       SHORT_NAME = 'O3',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                         &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'cloud_area_fraction',                                   &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'FCLD',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                         &
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

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'surface_albedo_for_visible_beam',                       &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBVR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'surface_albedo_for_visible_diffuse',                    &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBVF',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'surface_albedo_for_near_infrared_beam',                 &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBNR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = 'surface_albedo_for_near_infrared_diffuse',              &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBNF',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                                         ,&
       LONG_NAME  = 'surface_net_downward_shortwave_flux'                   ,&
       UNITS      = 'W m-2'                                                 ,&
       SHORT_NAME = 'SWNDSRF'                                               ,&
       DIMS       = MAPL_DimsHorzOnly                                       ,&
       VLOCATION  = MAPL_VLocationNone                                      ,&
                                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       SHORT_NAME = 'PREF',                                                  &
       LONG_NAME  = 'reference_air_pressure',                                &
       UNITS      = 'Pa',                                                    &
       DIMS       = MAPL_DimsVertOnly,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)


!  Solar does not have a "real" state. We keep an internal variable
!  for each variable produced by solar during the compute steps. 
!  Versions of these, weighted by the appropriate TOA insolation,
!  are returned at each time step.

!  !INTERNAL STATE:

    call MAPL_AddInternalSpec(GC,                                       &
       LONG_NAME  ='normalized_net_downward_shortwave_flux_in_air',          &
       UNITS      ='1',                                                      &
       SHORT_NAME ='FSWN',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                       &
       LONG_NAME  ='normalized_net_downward_shortwave_flux_in_air_assuming_clear_sky',&
       UNITS      ='1',                                                      &
       SHORT_NAME ='FSCN',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                       &
       LONG_NAME  = 'normalized_surface_downwelling_ultraviolet_beam_flux',  &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DRUVRN',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                       &
       LONG_NAME  = 'normalized_surface_downwelling_ultraviolet_diffuse_flux', &
       UNITS      = '1',                                                 &
       SHORT_NAME = 'DFUVRN',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                       &
       LONG_NAME  = 'normalized_surface_downwelling_par_beam_flux',          &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DRPARN',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                       &
       LONG_NAME  = 'normalized_surface_downwelling_par_diffuse_flux',       &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DFPARN',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                       &
       LONG_NAME  = 'normalized_surface_downwelling_nearinfrared_beam_flux', &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DRNIRN',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                       &
       LONG_NAME  = 'normalized_surface_downwelling_nearinfrared_diffuse_flux', &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DFNIRN',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                         &
       LONG_NAME  ='normalized_net_downward_shortwave_flux_in_air_assuming_no_aerosol', &
       UNITS      ='1',                                                      &
       SHORT_NAME ='FSWNAN',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                         &
       LONG_NAME  ='normalized_net_downward_shortwave_flux_in_air_assuming_clear_sky_and_no_aerosol',&
       UNITS      ='1',                                                      &
       SHORT_NAME ='FSCNAN',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)


!  !EXPORT STATE:
  

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  ='net_downward_shortwave_flux_in_air',                     &
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSW',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  ='net_downward_shortwave_flux_in_air_assuming_clear_sky',  &
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSC',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  ='net_downward_shortwave_flux_in_air_assuming_no_aerosol', &
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSWNA',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  ='net_downward_shortwave_flux_in_air_assuming_clear_sky_and_no_aerosol',&
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSCNA',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                       &
       LONG_NAME  = 'normalized_surface_downwelling_ultraviolet_beam_flux',  &
       UNITS      = '1',                                                 &
       SHORT_NAME = 'DRUVRN',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
       LONG_NAME  = 'normalized_surface_downwelling_ultraviolet_diffuse_flux', &
       UNITS      = '1',                                                 &
       SHORT_NAME = 'DFUVRN',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
       LONG_NAME  = 'normalized_surface_downwelling_par_beam_flux',          &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DRPARN',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
       LONG_NAME  = 'normalized_surface_downwelling_par_diffuse_flux',       &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DFPARN',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
       LONG_NAME  = 'normalized_surface_downwelling_nearinfrared_beam_flux', &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DRNIRN',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
       LONG_NAME  = 'normalized_surface_downwelling_nearinfrared_diffuse_flux', &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DFNIRN',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_downwelling_ultraviolet_beam_normal_flux',      &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DRNUVR',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_downwelling_par_beam_normal_flux',              &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DRNPAR',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_downwelling_nearinfrared_beam_normal_flux',     &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DRNNIR',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_downwelling_ultraviolet_beam_flux',             &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DRUVR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_downwelling_ultraviolet_diffuse_flux',          &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DFUVR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_downwelling_par_beam_flux',                     &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DRPAR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_downwelling_par_diffuse_flux',                  &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DFPAR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_downwelling_nearinfrared_beam_flux',            &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DRNIR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_downwelling_nearinfrared_diffuse_flux',         &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DFNIR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'cloud_area_fraction',                                   &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'FCLD',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'cloud_area_fraction_for_low_clouds',                    &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'CLDLO',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'cloud_area_fraction_for_middle_clouds',                 &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'CLDMD',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'cloud_area_fraction_for_high_clouds',                   &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'CLDHI',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'total_cloud_area_fraction',                             &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'CLDTT',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'optical_thickness_of_low_clouds',                       &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'TAULO',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'optical_thickness_of_middle_clouds',                    &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'TAUMD',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'optical_thickness_of_high_clouds(EXPORT)',              &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'TAUHI',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'optical_thickness_of_all_clouds',                       &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'TAUTT',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'optical_thickness_for_ice_clouds',                      &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'TAUCLI',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'optical_thickness_for_liquid_clouds',                   &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'TAUCLW',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_net_downward_shortwave_flux_assuming_clear_sky',&
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSCS',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_net_downward_shortwave_flux',                   &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSRS',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_net_downward_shortwave_flux_assuming_clear_sky_and_no_aerosol',&
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSCSNA',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_net_downward_shortwave_flux_assuming_no_aerosol',&
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSRSNA',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_incoming_shortwave_flux',                       &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSF',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_incoming_shortwave_flux_assuming_clear_sky',    &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSFC',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_incoming_shortwave_flux_assuming_clean_sky',    &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSFNA',                                               &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_incoming_shortwave_flux_assuming_clear_clean_sky',&
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSFCNA',                                              &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_outgoinging_shortwave_flux',                       &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSUF',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_outgoing_shortwave_flux_assuming_clear_sky',    &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSUFC',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_outgoing_shortwave_flux_assuming_clean_sky',    &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSUFNA',                                               &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_outgoing_shortwave_flux_assuming_clear_clean_sky',&
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSUFCNA',                                              &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)
!END CAR

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'toa_outgoing_shortwave_flux',                           &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'OSR',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'toa_outgoing_shortwave_flux_assuming_clear_sky',        &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'OSRCLR',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'toa_outgoing_shortwave_flux_no_aerosol',                &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'OSRNA',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'toa_outgoing_shortwave_flux_no_aerosol__clear_sky',     &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'OSRCNA',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)
!END CAR

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'toa_net_downward_shortwave_flux',                       &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSR',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'toa_net_downward_shortwave_flux_assuming_clear_sky',    &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSC',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'toa_net_downward_shortwave_flux_assuming_no_aerosol',   &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSRNA',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'toa_net_downward_shortwave_flux_assuming_clear_sky_and_no_aerosol',&
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSCNA',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'toa_incoming_shortwave_flux',                           &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRTP',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_albedo',                                        &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBEDO',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_albedo_for_visible_beam',                       &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBVR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_albedo_for_visible_diffuse',                    &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBVF',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_albedo_for_near_infrared_beam',                 &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBNR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'surface_albedo_for_near_infrared_diffuse',              &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBNF',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'cosine_of_the_solar_zenith_angle',                      &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'COSZ',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'mean_cosine_of_the_solar_zenith_angle',                 &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'MCOSZ',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'aerosol_optical_thickness_in_0.225-0.285_band',               &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TAUA1',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'aerosol_optical_thickness_in_0.175-0.225_0.285-0.300_band',&
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TAUA2',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'aerosol_optical_thickness_in_0.300-0.325_band',         &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TAUA3',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'aerosol_optical_thickness_in_0.325-0.4_band',           &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TAUA4',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'aerosol_optical_thickness_in_ 0.4-0.690_band',          &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TAUA5',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'aerosol_optical_thickness_in_0.690-1.220_band',         &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TAUA6',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'aerosol_optical_thickness_in_1.220-2.270_band',         &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TAUA7',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'aerosol_optical_thickness_in_2.270-3.850_band',         &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TAUA8',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'dust_optical_thickness_in_ 0.4-0.690_band',             &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TAUDU',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'salt_optical_thickness_in_ 0.4-0.690_band',             &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TAUSS',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'sulfate_optical_thickness_in_ 0.4-0.690_band',          &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TAUSO',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'black_carbon_optical_thickness_in_ 0.4-0.690_band',     &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TAUBC',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'organic_carbon_optical_thickness_in_ 0.4-0.690_band',   &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TAUOC',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'total_dust_optical_thickness_in_ 0.4-0.690_band',       &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TTAUDU',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'total_salt_optical_thickness_in_ 0.4-0.690_band',       &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TTAUSS',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'total_sulfate_optical_thickness_in_ 0.4-0.690_band',    &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TTAUSO',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'total_black_carbon_optical_thickness_in_ 0.4-0.690_band',&
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TTAUBC',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'total_organic_carbon_optical_thickness_in_ 0.4-0.690_band',&
       UNITS      = '1',                                                     &
       SHORT_NAME = 'TTAUOC',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'total_dust_aerosol_loading',                            &
       UNITS      = 'kg m-2',                                                &
       SHORT_NAME = 'TDUST',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'total_sea_salt_aerosol_loading',                        &
       UNITS      = 'kg m-2',                                                &
       SHORT_NAME = 'TSALT',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'total_sulfate_aerosol_loading',                         &
       UNITS      = 'kg m-2',                                                &
       SHORT_NAME = 'TSO4',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'total_black_carbon_aerosol_loading',                    &
       UNITS      = 'kg m-2',                                                &
       SHORT_NAME = 'TBC',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                         &
       LONG_NAME  = 'total_organic_carbon_aerosol_loading',                  &
       UNITS      = 'kg m-2',                                                &
       SHORT_NAME = 'TOC',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
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

!EOP

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC, name="REFRESH"   ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="-SORAD"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="-BALANCE"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="-MISC"     ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="UPDATE"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="-AEROSOLS" ,RC=STATUS)
    VERIFY_(STATUS)
    
! Set Run method and use generic init and final methods
! -----------------------------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)  
  end subroutine SetServices


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_SolarPertGridCompMod

