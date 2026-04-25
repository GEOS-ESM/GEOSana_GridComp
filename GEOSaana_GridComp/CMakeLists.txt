esma_set_this ()

set (alldirs
  GSI_GridComp
  GEOSgsi_Coupler
  )

esma_add_library (${this}
  SRCS GEOSaana_GridCompMod.F90
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES MAPL NCEP_sp_r8i4 NCEP_w3_r8i4 NCEP_bacio_r4i4 NCEP_bufr_r8i4
  )
