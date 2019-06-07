esma_set_this ()

set (SRCS
      cnop2spg.F90 
      m_admtlm.F90 
      m_cnop.F90
      m_eigen.f90
      m_fsens2pert.F90
      m_gAtlm.F90
      m_gAadm.F90
      m_gQadm.F90
      m_gQtlm.F90
      m_fdiff.F90
      m_initadj.F90
      m_jacobian.F90
      m_postStats.F90
      m_shtrans.f90
      m_svnorms.F90 
      m_getModConst.F
      m_pertenergy.F90
      m_pseudo.F90
      m_simsv.F90
      m_svprojs.F90 
      m_poles.F90
      m_setsvecs.F
      m_sVecDef.F
      m_svspectra.F90 
      m_pertutil.F90 
      op.F90
      pl_vert_sp.f 
      pl_spec_sp.f 
      svecdrv.F90
)

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES fvgcm GMAO_transf arpack
  )
if (EXTENDED_SOURCE)
  set_target_properties (${this} PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
endif()

target_compile_definitions (${this} PRIVATE REAL8 LSMH_off VERSION=\"fvgcm\")
target_compile_definitions (${this} PRIVATE SPMD  GFIO CHECKPOINTING LINUX)
string(REPLACE " " ";" tmp ${FREAL8})
foreach(flag ${tmp})
  target_compile_options (${this} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
endforeach()
