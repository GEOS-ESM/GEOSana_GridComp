if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set (CMAKE_Fortran_FLAGS_RELEASE "${FOPT3} ${ALIGNCOM}")
endif ()

ecbuild_add_executable(TARGET make_convinfo.x SOURCES make_convinfo.F90 m_convtmpl.F90 m_actvotype.F90 m_actvchan.F90 satinfo_util.F90)
ecbuild_add_executable(TARGET make_ozinfo.x SOURCES make_ozinfo.F90 m_oztmpl.F90 m_actvchan.F90 satinfo_util.F90)
ecbuild_add_executable(TARGET make_satinfo.x SOURCES main.F90 m_sitmpl.F90 m_actvchan.F90 satinfo_util.F90)
