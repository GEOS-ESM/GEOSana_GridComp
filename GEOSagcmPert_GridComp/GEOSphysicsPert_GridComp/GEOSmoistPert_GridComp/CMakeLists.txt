esma_set_this ()

set (SRCS 
  GEOS_MoistPertGridComp.F90
  cloud.F90
  cloud_ad.F90
  cloud_tl.F90
  convection.F90
  convection_ad.F90
  convection_tl.F90
  qsat_util.F90
  )
esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL_Base GEOS_Shared GEOS_PertShared fvdycorepert
  )
