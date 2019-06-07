esma_set_this ()

set (SRCS
  tools/fv_treat_da_inc.F90
  model/lin_cloud_microphys.F90
  model/fv_cmp.F90
  geos_utils/fill_corner.F90
  geos_utils/flow_proj.F90
  geos_utils/ghost_cubsph.F90
  geos_utils/grid_utils_mod.F90
  geos_utils/cub2latlon.F90
  geos_utils/cub2cub.F90
  model/fv_arrays.F90
  tools/fv_mp_mod.F90
  tools/fv_eta.F90
  tools/fv_timing.F90
  model/fv_grid_utils.F90
  tools/sorted_index.F90
  tools/fv_grid_tools.F90
  model/a2b_edge.F90
  model/fv_fill.F90
  tools/init_hydro.F90
  model/fv_mapz.F90
  tools/fv_surf_map.F90
  model/fv_update_phys.F90
  model/fv_sg.F90
  tools/test_cases.F90
  tools/fv_diagnostics.F90
  tools/external_ic.F90
  tools/fv_restart.F90
  model/fv_control.F90
  model/tp_core.F90
  model/sw_core.F90
  model/fastexp.c
  model/nh_utils.F90
  model/nh_core.F90
  model/dyn_core.F90
  model/fv_tracer2d.F90
  model/fv_dynamics.F90
  tools/fv_io.F90
  tools/sim_nc_mod.F90
  model/boundary.F90
  tools/external_sst.F90
  tools/fv_nudge.F90
  model/fv_nesting.F90
  model_tlmadm/a2b_edge_adm.F90
  model_tlmadm/a2b_edge_tlm.F90
  model_tlmadm/boundary_adm.F90
  model_tlmadm/boundary_tlm.F90
  model_tlmadm/dyn_core_adm.F90
  model_tlmadm/dyn_core_tlm.F90
  model_tlmadm/fv_arrays_nlm.F90
  model_tlmadm/fv_control_nlm.F90
  model_tlmadm/fv_dynamics_adm.F90
  model_tlmadm/fv_dynamics_tlm.F90
  model_tlmadm/fv_grid_utils_adm.F90
  model_tlmadm/fv_grid_utils_tlm.F90
  model_tlmadm/fv_mapz_adm.F90
  model_tlmadm/fv_mapz_tlm.F90
  model_tlmadm/fv_mp_adm.F90
  model_tlmadm/fv_mp_tlm.F90
  model_tlmadm/fv_nesting_adm.F90
  model_tlmadm/fv_nesting_tlm.F90
  model_tlmadm/fv_pressures_tlm.F90
  model_tlmadm/fv_pressures_adm.F90
  model_tlmadm/fv_restart_adm.F90
  model_tlmadm/fv_restart_tlm.F90
  model_tlmadm/fv_tracer2d_adm.F90
  model_tlmadm/fv_tracer2d_tlm.F90
  model_tlmadm/nh_core_adm.F90
  model_tlmadm/nh_core_tlm.F90
  model_tlmadm/nh_utils_adm.F90
  model_tlmadm/nh_utils_tlm.F90
  model_tlmadm/fv_sg_tlm.F90
  model_tlmadm/fv_sg_adm.F90
  model_tlmadm/sw_core_adm.F90
  model_tlmadm/sw_core_tlm.F90
  model_tlmadm/tp_core_adm.F90
  model_tlmadm/tp_core_tlm.F90
)

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES MAPL_Base GFDL_fms_r8 Tapenade
  )

target_compile_definitions (${this} PRIVATE MAPL_MODE SPMD TIMING USE_COND MAIST_CAPPA OLDMPP)
esma_fortran_generator_list(${this} ${FREAL8})


  
