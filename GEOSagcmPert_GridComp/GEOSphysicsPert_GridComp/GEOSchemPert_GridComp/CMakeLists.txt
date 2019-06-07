esma_set_this ()
set (ALLDIRS
  GEOSpchemPert_GridComp
  #  GOCART_PertGridComp
  #  StratChem_PertGridComp
  #  GMIchem_PertGridComp
  )

esma_add_library (${this}
  SRCS GEOS_ChemPertGridComp.F90
  DEPENDENCIES MAPL_Base
  SUBCOMPONENTS ${ALLDIRS}
  )

