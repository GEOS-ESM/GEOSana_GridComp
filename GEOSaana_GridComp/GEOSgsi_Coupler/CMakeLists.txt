esma_set_this ()

set (SRCSO
  geos_pertState.F90
  geos_pertStateIO.F90
  geos_StateIO.F90
  geos_pertmod.F90
  geos_nstmod.F90
  geos_gsi_util.F90
  g5_pertmod.F90
  cplr_pertmod.F90
  cplr_ensemble.F90
  cplr_nst.F90
  cplr_timermod.F90
#  stub_GEOSagcmPert_GridCompMod.F90
  )

set (SRCSX
  crtm_aerosol.F90 cplr_Set_CRTM_Aerosol.F90
  crtm_cloud.F90   cplr_Set_CRTM_Cloud.F90
  )

set (SRCS ${SRCSO} ${SRCSX})

# This is equivalent to FOPT= in GNU Make
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set (CMAKE_Fortran_FLAGS_RELEASE "${FOPT3} ${BIG_ENDIAN} ${BYTERECLEN} ${TRACEBACK} ${ALIGNCOM}")
endif ()

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GEOSagcmPert_GridComp GSI_GridComp Chem_Base nc_diag_write NCEP_crtm fvgcm MAPL fvgcm
  )

target_compile_definitions (${this} PRIVATE GEOS_PERT G5PERTMOD SAVE_DYN_PROG_STATE)
