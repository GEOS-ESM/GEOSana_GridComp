esma_set_this()

set (alldirs
  GEOSaana_GridComp
  GEOSoana_GridComp
  )

esma_add_library(${this}
  SRCS GEOS_AnaGridComp.F90
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES MAPL
  INCLUDES ${INC_ESMF})


