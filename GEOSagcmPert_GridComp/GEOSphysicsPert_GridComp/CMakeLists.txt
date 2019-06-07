esma_set_this ()

set (ALLDIRS
  GEOSchemPert_GridComp
  GEOSmoistPert_GridComp
  GEOSsurfacePert_GridComp
  GEOSturbulencePert_GridComp
  GEOSgwdPert_GridComp
  GEOSradiationPert_GridComp
  )


esma_add_library (${this}
  SRCS GEOS_PhysicsPertGridComp.F90
  DEPENDENCIES MAPL_Base
  SUBCOMPONENTS ${ALLDIRS}
)
