#esma_set_this ()

add_subdirectory(mksi)
add_subdirectory(etc)

set (GSIGC_SRCS
  GSI_GridCompMod.F90
  )


set (SRCSX
  m_gsiversion.F90
  m_FletcherChecksum.F90
  )

set (SRCS_OBSVR
  convb_ps.f90
  convb_q.f90
  convb_t.f90
  convb_uv.f90
  converr_ps.f90
  converr_pw.f90
  converr_q.f90
  converr_t.f90
  converr_uv.f90
  convinfo.f90
  correlated_obsmod.F90
  obsmod.F90
  ozinfo.f90
  radinfo.f90
 )
  
set (SRCS_UTIL
  aircraftinfo.f90
  aircraftobsqc.f90
  aniso_ens_util.f90
  antcorr_application.f90
  blendmod.f90
  cloud_efr_mod.f90
  compact_diffs.f90
  compute_fact10.f90
  constants.f90
  control_vectors.f90
  derivsmod.f90
  egrid2agrid_mod.f90
  ens_spread_mod.f90
  fgrid2agrid_mod.f90
  fill_mass_grid2.f90
  fill_nmm_grid2.f90
  fpvsx_ad.f90
  general_commvars_mod.f90
  general_read_gfsatm.f90
  general_read_nmmb.f90 
  general_specmod.f90
  general_spectral_transforms.f90
  general_sub2grid_mod.f90
  general_transform.f90
  general_tll2xy_mod.f90
  general_write_gfsatm.f90
  gengrid_vars.f90
  genqsat.f90
  get_derivatives.f90
  get_derivatives2.f90
  get_gefs_ensperts_dualres.f90
  get_nmmb_ensperts.f90
  get_pseudo_ensperts.F90
  get_wrf_mass_ensperts_netcdf.F90
  get_wrf_nmm_ensperts.F90
  get_semimp_mats.f90
  getprs.f90
  getuv.f90
  gsi_4dvar.f90
  gsi_bias.f90
  gsi_bundlemod.F90
  gsi_chemguess_mod.F90
  gsi_metguess_mod.F90
  gsi_rfv3io_mod.f90
  gsi_unformatted.F90
  grdcrd.f90
  gridmod.F90
  gscond_ad.f90
  gsi_4dcouplermod.f90
  gsi_enscouplermod.f90
  gsi_nstcouplermod.f90
  gsi_io.f90
  gsi_nemsio_mod.f90
  guess_grids.F90
  half_nmm_grid2.f90
  hilbert_curve.f90
  hybrid_ensemble_parameters.f90
  insitu_info.f90
  intrp_msk.f90
  jcmod.f90
  jfunc.f90
  kinds.F90
  looplimits.f90
  m_dgeevx.F90
  m_rerank.f90
  m_stats.f90
  m_tick.F90
  mpeu_mpif.F90
  mpeu_util.F90
  mod_fv3_lola.f90
  mod_nmmb_to_a.f90
  mod_strong.f90
  mod_wrfmass_to_a.f90
  mp_compact_diffs_mod1.f90
  mp_compact_diffs_support.f90
  mpimod.F90
  mpl_allreduce.F90
  mpl_bcast.f90
  native_endianness.f90
  netcdf_mod.f90
  ncepgfs_io.f90
  ncepgfs_ghg.f90
  ncepnems_io.f90
  nstio_module.f90
  Nst_Var_ESMFMod.f90
  nlmsas_ad.f90
  omegas_ad.f90
  phil.f90
  phil1.f90
  plib8.f90
  polcarf.f90
  prad_bias.f90
  precond.f90
  precpd_ad.f90
  psichi2uv_reg.f90
  psichi2uvt_reg.f90
  raflib.f90
  rapidrefresh_cldsurf_mod.f90
  rdgrbsst.f90
  read_diag.f90
  rfdpar.f90
  rsearch.F90
  rtlnmc_version3.f90
  satthin.F90
  set_crtm_aerosolmod.f90
  set_crtm_cloudmod.f90
  sfc_model.f90
  simpin1.f90
  simpin1_init.f90
  smoothwwrf.f90
  sparsearr.f90
  state_vectors.f90
  stop1.f90
  stub_ensmod.f90
  stub_nstmod.f90
  stub_set_crtm_aerosol.f90
  stub_pertmod.F90
  stub_timermod.f90
  tendsmod.f90
  timermod.f90
  tintrp2a.f90
  tintrp3.f90
  tv_to_tsen.f90
  tune_pbl_height.f90
  turbl.f90
  turbl_ad.f90
  turbl_tl.f90
  turblmod.f90
  unfill_mass_grid2.f90
  unfill_nmm_grid2.f90
  unhalf_nmm_grid2.f90
  wind_fft.f90
  blockIO.c
  ${SRCSX}
  )

set (SRCS_OTHER
  get_gefs_for_regional.f90
  gsd_update_mod.f90
  read_wrf_mass_files.F90
  read_wrf_mass_guess.F90
  read_wrf_nmm_files.F90
  read_wrf_nmm_guess.F90
  regional_io.f90
  wrf_binary_interface.F90
  wrf_mass_guess_mod.f90
  wrf_netcdf_interface.F90
  wrwrfmassa.F90
  wrwrfnmma.F90
  )

set (SRCS_OBSNODES
  m_aerolNode.F90
  m_aeroNode.F90
  m_cldchNode.F90
  m_colvkNode.F90
  m_dwNode.F90
  m_gpsNode.F90
  m_gustNode.F90
  m_howvNode.F90
  m_lagNode.F90
  m_lcbasNode.F90
  m_lwcpNode.F90
  m_mitmNode.F90
  m_mxtmNode.F90
  m_o3lNode.F90
  m_ozNode.F90
  m_pblhNode.F90
  m_pcpNode.F90
  m_pm10Node.F90
  m_pm2_5Node.F90
  m_pmslNode.F90
  m_psNode.F90
  m_pwNode.F90
  m_qNode.F90
  m_radNode.F90
  m_rwNode.F90
  m_spdNode.F90
  m_sstNode.F90
  m_swcpNode.F90
  m_tcamtNode.F90
  m_td2mNode.F90
  m_tcpNode.F90
  m_tNode.F90
  m_uwnd10mNode.F90
  m_visNode.F90
  m_vwnd10mNode.F90
  m_wNode.F90
  m_wspd10mNode.F90
  )


set (SRCS_SOLVER
  adjtest.f90
  adjtest_obs.f90
  adjust_cloudobs_mod.f90
  aeroinfo.f90
  allsky_mod.f90
  anberror.f90
  anbkerror.f90
  anisofilter.f90
  anisofilter_glb.f90
  antest_maps0.f90
  antest_maps0_glb.f90
  atms_spatial_average_mod.f90
  balmod.f90
  berror.f90
  bias_predictors.f90
  bkerror.f90
  bkgcov.f90
  bkgvar.f90
  bkgvar_rewgt.f90
  bicg.f90
  bicglanczos.F90
  blacklist.f90
  buddycheck_mod.f90
  calc_fov_conical.f90
  calc_fov_crosstrk.f90
  calctends.f90
  calctends_ad.f90
  calctends_tl.f90
  calctends_no_ad.f90
  calctends_no_tl.f90
  chemmod.f90
  clw_mod.f90
  cmaq_routines.f90
  co_mop_ak.f90
  coinfo.f90
  combine_radobs.f90
  compute_derived.f90
  compute_qvar3d.f90
  control2model.f90
  control2model_ad.f90
  control2state.f90
  control2state_ad.f90
  converr.f90
  convthin.f90
  convthin_time.f90
  crtm_interface.f90
  cvsection.f90
  cwhydromod.f90
  dtast.f90
  deter_sfc_mod.f90
  en_perts_io.f90
  enorm_state.f90
  ensctl2model.f90
  ensctl2model_ad.f90
  ensctl2state.f90
  ensctl2state_ad.f90
  evaljgrad.f90
  evaljo.f90
  evalqlim.f90
  genstats_gps.f90
  gesinfo.f90
  getcount_bufr.f90
  getsiga.f90
  getvvel.f90
  gfs_stratosphere.f90
  glbsoi.f90
  grtest.f90
  gsd_terrain_match_surfTobs.f90
  gsdcloudanalysis4NMMB.F90
  gsdcloudanalysis.F90
  gsimod.F90
  gsisub.F90
  hybrid_ensemble_isotropic.F90
  inc2guess.f90
  init_jcdfi.f90
  int3dvar.f90
  intall.f90
  intaod.f90
  intcldch.f90
  intco.f90
  intdw.f90
  intgps.f90
  intgust.f90
  inthowv.f90
  intjcmod.f90
  intjo.f90
  intlag.f90
  intlcbas.f90
  intlwcp.f90
  intmitm.f90
  intmxtm.f90
  intoz.f90
  intpblh.f90
  intpcp.f90
  intpm10.f90
  intpm2_5.f90
  intpmsl.f90
  intps.f90
  intpw.f90
  intq.f90
  intrad.f90
  intrp2a.f90
  intrp3oz.f90
  intrw.f90
  intspd.f90
  intwspd10m.f90
  intsst.f90
  intswcp.f90
  intt.f90
  inttcamt.f90
  inttd2m.f90
  inttcp.f90
  intuwnd10m.f90
  intvis.f90
  intvwnd10m.f90
  intw.f90
  jgrad.f90
  lag_fields.f90
  lag_interp.f90
  lag_traj.f90
  lagmod.f90
  lanczos.F90
  loglcbas_to_lcbas.f90
  logcldch_to_cldch.f90
  logvis_to_vis.f90
  m_berror_stats.f90
  m_berror_stats_reg.f90
  m_cvgridLookup.F90
  m_extOzone.F90
  m_gsiBiases.f90
  m_distance.f90
  m_dtime.F90
  m_find.f90
  m_gpsrhs.F90
  m_latlonRange.F90
  m_rhs.F90
  m_sortind.f90
  m_uniq.f90
  m_obsdiagNode.F90
  m_obsdiags.F90
  m_obsHeadBundle.F90
  m_obsNode.F90
#  m_obsNode_scalar2D.F90
  m_obsLList.F90
  m_obsNodeTypeManager.F90
  mod_vtrans.f90
  model_ad.F90
  model_tl.F90
#  model2control.f90
  normal_rh_to_q.f90
  obserr_allsky_mw.f90
  obs_ferrscale.F90
  obs_para.f90
  obs_sensitivity.f90
  observer.F90
  oneobmod.F90
  patch2grid_mod.f90
  pcgsoi.f90
  pcgsqrt.f90
  pcp_k.f90
  pcpinfo.f90
  penal.f90
  prewgt.f90
  prewgt_reg.f90
  projmethod_support.f90
  prt_guess.f90
  q_diag.f90
  qcmod.f90
  qnewton3.f90
  radiance_mod.f90
  read_aerosol.f90
  read_ahi.f90
  read_airs.f90
  read_amsre.f90
  read_amsr2.f90
  read_anowbufr.f90
  read_atms.f90
  read_avhrr.f90
  read_avhrr_navy.f90
  read_bufrtovs.f90
  read_co.f90
  read_files.f90
  read_fl_hdob.f90
  read_gfs_ozone_for_regional.f90
  read_goesimg.f90
  read_goesndr.f90
  read_goesimgr_skycover.f90
  read_gmi.f90
  read_gps.f90
  read_guess.F90
  read_iasi.f90
  read_cris.f90
  read_l2bufr_mod.f90
  read_lag.f90
  read_lidar.f90
  read_Lightning.f90
  read_nsstbufr.f90
  read_mitm_mxtm.f90
  read_modsbufr.f90
  read_nasa_larc.f90
  read_NASA_LaRC_cloud.f90
  read_obs.F90
  read_ozone.f90
  read_pblh.f90
  read_pcp.f90
  read_prepbufr.f90
  read_radar.f90
  read_radarref_mosaic.f90
  read_rapidscat.f90
  read_saphir.f90
  read_satmar.f90
  read_satwnd.f90
  read_seviri.f90
  read_sfcwnd.f90
  read_ssmi.f90
  read_ssmis.f90
  read_tcps.f90
  read_wcpbufr.f90
  reorg_metar_cloud.f90
  setupaod.f90
  setupbend.f90
  setupcldch.f90
  setupco.f90
  setupdw.f90
  setupgust.f90
  setuphowv.f90
  setuplag.f90
  setuplcbas.f90
  setuplwcp.f90
  setupmitm.f90
  setupmxtm.f90
  setuppmsl.f90
  setupoz.f90
  setuppblh.f90
  setuppcp.f90
  setuppm10.f90
  setuppm2_5.f90
  setupps.f90
  setuppw.f90
  setupq.f90
  setuprad.f90
  setupref.f90
  setuprhsall.f90
  setuprw.f90
  setupspd.f90
  setupwspd10m.f90
  setupsst.f90
  setupswcp.f90
  setupt.f90
  setuptcamt.f90
  setuptd2m.f90
  setuptcp.f90
  setupuwnd10m.f90
  setupvis.f90
  setupvwnd10m.f90
  setupw.f90
  sfcobsqc.f90
  smooth_polcarf.f90
  smoothrf.f90
  smoothzrf.f90
  sqrtmin.f90
  ssmis_spatial_average_mod.f90
  sst_retrieval.f90
  statsco.f90
  statsconv.f90
  statsoz.f90
  statspcp.f90
  statsrad.f90
  stp3dvar.f90
  stpaod.f90
  stpcalc.f90
  stpcldch.f90
  stpco.f90
  stpdw.f90
  stpgps.f90
  stpgust.f90
  stphowv.f90
  stpjo.f90
  stpjcmod.f90
  stpoz.f90
  stplcbas.f90
  stplwcp.f90
  stpmitm.f90
  stpmxtm.f90
  stppmsl.f90
  stppblh.f90
  stppcp.f90
  stppm10.f90
  stppm2_5.f90
  stpps.f90
  stppw.f90
  stpq.f90
  stprad.f90
  stprw.f90
  stpspd.f90
  stpsst.f90
  stpswcp.f90
  stpt.f90
  stptcamt.f90
  stptd2m.f90
  stptcp.f90
  stpuwnd10m.f90
  stpvis.f90
  stpvwnd10m.f90
  stpw.f90
  stpwspd10m.f90
  strong_bal_correction.f90
  strong_baldiag_inc.f90
  strong_fast_global_mod.f90
  sub2fslab_mod.f90
  support_2dvar.f90
  test_obsens.f90
  tcv_mod.f90
  tpause.f90
  tpause_t.F90
  update_guess.f90
  update_geswtend.f90
  write_all.F90
  write_bkgvars_grid.f90
  xhat_vordivmod.f90
  zrnmi_mod.f90
  ${SRCS_OBSNODES}
  )

set (SRCS
  ${SRCS_UTIL} ${SRCS_OBSVR} ${SRCS_OTHER} ${SRCS_SOLVER} ${GSIGC_SRCS} prepbykx.f
  )

set (OBJS_OPENBUFR
   getcount_bufr.f90
   m_extOzone.F90
   read_aerosol.f90
   read_ahi.f90
   read_airs.f90
   read_amsre.f90
   read_amsr2.f90
   read_anowbufr.f90
   read_atms.f90
   read_avhrr.f90
   read_avhrr_navy.f90
   read_bufrtovs.f90
   read_cris.f90
   read_fl_hdob.f90
   read_goesimg.f90
   read_goesndr.f90
   read_goesimgr_skycover.f90
   read_gmi.f90
   read_gps.f90
   read_iasi.f90
   read_l2bufr_mod.f90
   read_lidar.f90
   read_Lightning.f90
   read_mitm_mxtm.f90
   read_modsbufr.f90
   read_nsstbufr.f90
   read_nasa_larc.f90
   read_NASA_LaRC_cloud.f90
   read_ozone.f90
   read_pblh.f90
   read_pcp.f90
   read_prepbufr.f90
   read_radar.f90
   read_radarref_mosaic.f90
   read_rapidscat.f90
   read_saphir.f90
   read_satmar.f90
   read_satwnd.f90
   read_seviri.f90
   read_sfcwnd.f90
   read_ssmi.f90
   read_ssmis.f90
   read_wcpbufr.f90
   oneobmod.F90
   )

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   if (CMAKE_BUILD_TYPE MATCHES Release)
      set (CMAKE_Fortran_FLAGS_RELEASE "")
      foreach (src ${SRCS})
         if (${src} MATCHES blockIO.c)
            set_source_files_properties (${src} PROPERTIES COMPILE_FLAGS "")
         elseif (${src} MATCHES prepbykx.f)
            set_source_files_properties (${src} PROPERTIES COMPILE_FLAGS "${GEOS_Fortran_Release_Flags} ${common_Fortran_flags} ${BYTERECLEN} ${EXTENDED_SOURCE} ${FP_MODEL_STRICT} ${ALIGNCOM}")
         elseif (${src} MATCHES m_gsiversion.F90)
            set_source_files_properties (${src} PROPERTIES COMPILE_FLAGS "${GEOS_Fortran_Release_Flags} ${common_Fortran_flags} ${BIG_ENDIAN} ${BYTERECLEN}")
         elseif (${src} MATCHES crtm_interface.f90)
            set_source_files_properties (${src} PROPERTIES COMPILE_FLAGS "${BIG_ENDIAN} ${BYTERECLEN} ${TRACEBACK} ${FOPT3}")
         elseif (${src} MATCHES hybrid_ensemble_isotropic.F90)
            set_source_files_properties (${src} PROPERTIES COMPILE_FLAGS "${BIG_ENDIAN} ${BYTERECLEN} ${TRACEBACK} ${FOPT1}")
         elseif (${src} MATCHES read_diag.f90)
            set_source_files_properties (${src} PROPERTIES COMPILE_FLAGS "${BIG_ENDIAN} ${BYTERECLEN} ${TRACEBACK} ${FOPT0}")
         elseif (${src} IN_LIST OBJS_OPENBUFR)
            set_source_files_properties (${src} PROPERTIES COMPILE_FLAGS "${GEOS_Fortran_Release_Flags} ${common_Fortran_flags} ${BYTERECLEN} ${FP_MODEL_STRICT} ${ALIGNCOM}")
         else ()
            set_source_files_properties (${src} PROPERTIES COMPILE_FLAGS "${GEOS_Fortran_Release_Flags} ${common_Fortran_flags} ${BIG_ENDIAN} ${BYTERECLEN} ${FP_MODEL_STRICT} ${ALIGNCOM}")
         endif ()
      endforeach ()
   endif ()
endif ()

set (gsiscripts
     analyzer gsidiags match_obcls_obsys.pl)

install(
   PROGRAMS ${gsiscripts} 
   DESTINATION bin
   )
  
esma_add_library (GSI_GridComp
  SRCS ${SRCS_UTIL} ${SRCS_OBSVR} ${SRCS_OTHER} ${SRCS_SOLVER} ${GSIGC_SRCS} 
  DEPENDENCIES MAPL_Base NCEP_crtm NCEP_nemsio GMAO_stoch nc_diag_read nc_diag_write NCEP_bufr_r8i4)

foreach (item HAVE_ESMF _SKIP_READ_OBS_CHECK_ _TIMER_ON_ _REAL8_)
   target_compile_definitions (GSI_GridComp PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${item}>)
endforeach()

foreach (item LINUX funder FortranByte=char FortranInt=int FortranLlong=long\ long)
   target_compile_definitions (GSI_GridComp PRIVATE $<$<COMPILE_LANGUAGE:C>:${item}>)
endforeach()

# MAT: I am not sure we can do the split libraries like GNU Make GEOS used to.
#      As an example, obsmod.F90 (Observer) depends on constants.f90 (Util). But
#      aircraftinfo.f90 (Util) depends on obsmod.F90 (Observer). This means 
#      GSI_Util depends on GSI_Observer and GSI_Observer depends on GSI_Util.

#esma_add_library (GSI_Util
#  SRCS ${SRCS_UTIL}
#  DEPENDENCIES GSI_Observer MAPL_Base NCEP_crtm NCEP_nemsio nc_diag_read GMAO_stoch)  
#
#esma_add_library (GSI_Observer
#  SRCS ${SRCS_OBSVR}
#  DEPENDENCIES MAPL_Base GSI_Util nc_diag_write)
#
#esma_add_library (GSI_Other
#  SRCS ${SRCS_OTHER}
#  DEPENDENCIES MAPL_Base GSI_Util nc_diag_write)
#
#esma_add_library (GSI_Solver
#  SRCS ${SRCS_SOLVER}
#  DEPENDENCIES MAPL_Base)
#
#foreach (t GSI_Util GSI_Observer GSI_Other GSI_Solver)
#  target_compile_definitions (${t} PUBLIC _REAL8_ _SKIP_READ_OBS_CHECK_ LINUX funder FortranByte=char FortranInt=int FortranLlong=long\ long)
#endforeach()

ecbuild_add_executable(TARGET prepbykx.x SOURCES prepbykx.f LIBS ${this})

set (etc_txt_files
   etc/gmao_global_aeroinfo.txt
   etc/gmao_global_anavinfo.txt
   etc/gmao_global_anavinfo_l132.txt
   etc/gmao_global_anavinfo_rcov.txt
   etc/gmao_global_blacklist.txt
   etc/gmao_global_cloudy_radiance_info.txt
   etc/gmao_global_coinfo.txt
   etc/gmao_global_convinfo.txt
   etc/gmao_global_insituinfo.txt
   etc/gmao_global_ozinfo.txt
   etc/gmao_global_pcpinfo.txt
   etc/gmao_global_satinfo.txt
   etc/gmao_global_scaninfo.txt
   etc/gmao_global_tlapmean.txt
   etc/gsi_atms_beamwidth.txt
   )

foreach (txt_file ${etc_txt_files})
   get_filename_component(barename ${txt_file} NAME_WE)
   set (tmp_name "${barename}.rc")
   install(
      FILES ${txt_file}
      DESTINATION etc
      RENAME ${tmp_name}
      )
endforeach ()

# MAT Have to handle some of the files above differently due to
#     the extra dots and I can't figure out how to do it programatically

install(
   FILES etc/gmao_global_hybens_info.x288y181l72.txt
   DESTINATION etc
   RENAME gmao_global_hybens_info.x288y181l72.rc
   )
install(
   FILES etc/gmao_global_hybens_info.x576y361l72.txt
   DESTINATION etc
   RENAME gmao_global_hybens_info.x576y361l72.rc
   )

install(
   FILES etc/gmao_global_hybens_info.x288y181l132.txt
   DESTINATION etc
   RENAME gmao_global_hybens_info.x288y181l132.rc
   )
install(
   FILES etc/gmao_global_hybens_info.x576y361l132.txt
   DESTINATION etc
   RENAME gmao_global_hybens_info.x576y361l132.rc
   )

set (etc_other_files
   etc/gmao_airs_bufr.tbl
   etc/gsi_fdda_1.rc.tmpl
   etc/gsi_fdda_2.rc.tmpl
   etc/gsi_fgat_1.rc.tmpl
   etc/gsi_fgat_2.rc.tmpl
   etc/gsi.rc.tmpl
   etc/gsi_sens.rc.tmpl
   etc/gsidiags.rc
   etc/obs.rc.tmpl
   )
install (
   FILES ${etc_other_files}
   DESTINATION etc
   )

file(GLOB sidb_files mksi/sidb/*)
install(
   FILES ${sidb_files}
   DESTINATION etc/gmao_satinfo.db
   )

file(GLOB ozinfo_files mksi/ozinfo.db/*)
install(
   FILES ${ozinfo_files}
   DESTINATION etc/gmao_ozinfo.db
   )

file(GLOB convinfo_files mksi/convinfo.db/*)
install(
   FILES ${convinfo_files}
   DESTINATION etc/gmao_convinfo.db
   )

set (MERRA2_files
   etc/MERRA2/gmao_global_satinfo.rc
   etc/MERRA2/gmao_global_scaninfo.rc
   etc/MERRA2/gsi.rc.tmpl
   etc/MERRA2/gsi_sens.rc.tmpl
   etc/MERRA2/sac.nl.tmpl
   )

install(
   FILES ${MERRA2_files}
   DESTINATION etc/gsi/MERRA2
   )

file(GLOB etc_py etc/*.py)
install (
   PROGRAMS ${etc_py}
   DESTINATION bin
   )


