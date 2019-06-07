esma_set_this ()

set (SRCS
  GEOS_RadiationPertGridComp.F90
  MieInterface.F90
  cloudradcouple.F90
  irrad.F90
  irrad_ad.F90
  irrad_tl.F90
  irradconstants.F90
  radconstants.F90
  sorad.F90
  sorad_ad.F90
  sorad_tl.F90
  soradconstants.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL_Base Chem_Base GEOS_Shared GEOS_PertShared GEOSdynamicsPert_GridComp
  )
