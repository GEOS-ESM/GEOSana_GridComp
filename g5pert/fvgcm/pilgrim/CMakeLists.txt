target_sources (fvgcm PRIVATE
  mod_comm_adtl.F90  mod_comm.F90  mod_pio.F90
  memstuff.c  mp_assign_to_cpu.c  pio_util.c
)

target_include_directories (fvgcm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
