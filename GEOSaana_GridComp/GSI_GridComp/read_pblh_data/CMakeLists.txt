esma_set_this()

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set (CMAKE_Fortran_FLAGS_RELEASE "${FOPT3} ${ALIGNCOM}")
endif ()

ecbuild_add_executable(TARGET combine_pblh_mplnet.x SOURCES combine_pblh_mplnet.f90 LIBS GSI_GridComp NetCDF::NetCDF_Fortran OpenMP::OpenMP_Fortran) 
set_target_properties(combine_pblh_mplnet.x PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})
target_compile_definitions (combine_pblh_mplnet.x PRIVATE _REAL8_)

ecbuild_add_executable(TARGET combine_pblh_wp.x SOURCES combine_pblh_wp.f90 LIBS GSI_GridComp NetCDF::NetCDF_Fortran OpenMP::OpenMP_Fortran) 
set_target_properties(combine_pblh_wp.x PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})
target_compile_definitions (combine_pblh_wp.x PRIVATE _REAL8_)

ecbuild_add_executable(TARGET read_binary.x SOURCES read_binary.f90 LIBS GSI_GridComp NetCDF::NetCDF_Fortran OpenMP::OpenMP_Fortran) 
set_target_properties(read_binary.x PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})
target_compile_definitions (read_binary.x PRIVATE _REAL8_)

ecbuild_add_executable(TARGET combine_pblh_cats_cnn.x SOURCES combine_pblh_cats_cnn.f90 LIBS GSI_GridComp NetCDF::NetCDF_Fortran OpenMP::OpenMP_Fortran) 
set_target_properties(combine_pblh_cats_cnn.x PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})
target_compile_definitions (combine_pblh_cats_cnn.x PRIVATE _REAL8_)

ecbuild_add_executable(TARGET combine_pblh_cats_dtds.x SOURCES combine_pblh_cats_dtds.f90 LIBS GSI_GridComp NetCDF::NetCDF_Fortran OpenMP::OpenMP_Fortran) 
set_target_properties(combine_pblh_cats_dtds.x PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})
target_compile_definitions (combine_pblh_cats_dtds.x PRIVATE _REAL8_)

ecbuild_add_executable(TARGET combine_pblh_icesat2_dtds.x SOURCES combine_pblh_icesat2_dtds.f90 LIBS GSI_GridComp NetCDF::NetCDF_Fortran OpenMP::OpenMP_Fortran)
set_target_properties(combine_pblh_icesat2_dtds.x PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})
target_compile_definitions (combine_pblh_icesat2_dtds.x PRIVATE _REAL8_)

ecbuild_add_executable(TARGET combine_pblh_icesat2_cnn.x SOURCES combine_pblh_icesat2_cnn.f90 LIBS GSI_GridComp NetCDF::NetCDF_Fortran OpenMP::OpenMP_Fortran)
set_target_properties(combine_pblh_icesat2_cnn.x PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})
target_compile_definitions (combine_pblh_icesat2_cnn.x PRIVATE _REAL8_)
