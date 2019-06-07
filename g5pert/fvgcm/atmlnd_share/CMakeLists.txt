target_sources (fvgcm PRIVATE
  datetime.F90  gauaw_mod.F90  ioFileMod.F90     
  mpishorthand.F  string_utils.F90 endrun.F90    
  infnan.F90     linebuf_stdout.c  precision.F90 wrap_mpi.F90
)

target_include_directories (fvgcm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
