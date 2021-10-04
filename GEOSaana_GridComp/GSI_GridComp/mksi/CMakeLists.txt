esma_set_this()

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set (CMAKE_Fortran_FLAGS_RELEASE "${FOPT3} ${ALIGNCOM}")
endif ()

  esma_add_library(${this}
    SRCS  m_actvchan.F90 satinfo_util.F90
    )

ecbuild_add_executable(TARGET make_convinfo.x SOURCES make_convinfo.F90 m_convtmpl.F90 m_actvotype.F90 LIBS ${this}) 
set_target_properties(make_convinfo.x PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})
ecbuild_add_executable(TARGET make_ozinfo.x SOURCES make_ozinfo.F90 m_oztmpl.F90 LIBS ${this})
set_target_properties(make_ozinfo.x PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})
ecbuild_add_executable(TARGET make_satinfo.x SOURCES main.F90 m_sitmpl.F90 LIBS ${this})
set_target_properties(make_satinfo.x PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})
ecbuild_add_executable(TARGET make_pcscaling.x SOURCES pcscaling_main.F90 m_pcscaling.F90 satinfo_util.F90)
set_target_properties(make_pcscaling.x PROPERTIES Fortran_MODULE_DIRECTORY ${esma_include}/${this})

install (PROGRAMS gsiinfo.pl DESTINATION bin)

install(
   DIRECTORY sidb/
   DESTINATION etc/gmao_satinfo.db
   )

install(
   DIRECTORY ozinfo.db/
   DESTINATION etc/gmao_ozinfo.db
   )

install(
   DIRECTORY convinfo.db/
   DESTINATION etc/gmao_convinfo.db
   )

install(
   DIRECTORY pcscaling.db/
   DESTINATION etc/gmao_pcscaling.db
   )
