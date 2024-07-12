esma_set_this()


if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set (CMAKE_Fortran_FLAGS_RELEASE "${FOPT3} ${ALIGNCOM}")
endif ()

ecbuild_add_executable(TARGET raob_pblh.x SOURCES calculate_pblh_from_raob_diag.f90 LIBS GSI_GridComp GEOSgsi_Coupler NetCDF::NetCDF_Fortran MPI::MPI_Fortran)

set_target_properties(raob_pblh.x PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})

# CMake has an OpenMP issue with NAG Fortran: https://gitlab.kitware.com/cmake/cmake/-/issues/21280
if (NOT CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
   target_link_libraries(raob_pblh.x OpenMP::OpenMP_Fortran)
endif ()

target_compile_definitions(raob_pblh.x PRIVATE _REAL8_)

