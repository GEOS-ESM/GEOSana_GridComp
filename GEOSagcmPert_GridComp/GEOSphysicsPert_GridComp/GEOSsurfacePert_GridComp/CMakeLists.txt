esma_set_this ()

set (ALLDIRS
#  GEOSlandicePert_GridComp
#  GEOSlakePert_GridComp
#  GEOSlandPert_GridComp
#  GEOSsaltwaterPert_GridComp
  )

esma_add_library (${this}
  SRCS GEOS_SurfacePertGridComp.F90
  DEPENDENCIES MAPL_Base
#  SUBCOMPONENTS ${ALLDIRS}
  )
