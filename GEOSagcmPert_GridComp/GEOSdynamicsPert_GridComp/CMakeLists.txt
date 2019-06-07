esma_set_this ()

set (SRCS
  AppGridCreate.F90
  CubeHalo.F90
  CubeToCubeRegridder.F90
  CubeToLatLon.F90
  CubeToLatLonRegridder.F90
  CubedSphereGridFactory.F90
  DynCorePert_GridComp.F90
  GetWeights.F90
  GetWeightsC2C.F90
  LatLonToCubeRegridder.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL_Base GEOS_PertShared
  SUBCOMPONENTS fvdycorepert
  )

