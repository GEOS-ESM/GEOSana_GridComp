target_sources (fvgcm PRIVATE
  avgc.F90     controlq.F90   dependentq.F90   func.F90    initfunc.F90     pkez.F90                  sw_core.F90  vpol5.F90
  benergy.F90  d2a2.F90       d_split.F90      fvcore.F90  mapz_module.F90  nlmodel.F90  polavg.F90       setfunc.F90   tp_core.F90  xxall.F90
  cd_core.F90  d2a3d.F90      fill_module.F90  geopk.F90   m_iostate.F      pft2d.F90    postfunc.F90     setrig.F90     trac2d.F90
  control.F90  dependent.F90  fillz.F90        highp2.F90  m_pftSV.F90      pft_cf.F90   prognostics.F90  stepon.F       upol5.F90
  m_trajmng.F m_trjphys.F90   a2d3d.F90        m_delp2ps.F90  prognostics_q.F90
)

target_include_directories (fvgcm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
