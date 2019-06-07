esma_set_this ()

set (SRCS
  GEOS_GwdPertGridComp.F90
  gw_drag_b.F90
  gw_drag_d.F90
  gw_drag_init.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL_Base GEOS_Shared GEOS_PertShared fvdycorepert
  )

