esma_set_this ()

esma_add_library (${this}
  SRCS GEOS_TurbulencePertGridComp.F90 bldriver.F90 blsimp.F90
  DEPENDENCIES MAPL_Base GEOS_Shared fvdycorepert GEOS_PertShared
  )
