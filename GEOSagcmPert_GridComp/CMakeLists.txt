esma_set_this ()

esma_add_library (${this}
  SRCS GEOS_AgcmPertGridComp.F90 GEOS_AgcmPertADTester.F90
  DEPENDENCIES MAPL_Base
  SUBCOMPONENTS  GEOSdynamicsPert_GridComp GEOSphysicsPert_GridComp
  )
