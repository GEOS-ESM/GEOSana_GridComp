target_sources (fvgcm PRIVATE
  age_of_air.F90  diaginit.F90  findid.F90     init2dz.F90   merge.F90       ppme.F90        timingModule.F90  wrt3d.F90   zout.F90
  aoa_wrt.F90     diagout.F90   gmean4.F90     initmem.F90   minmax.F90      refout.F90      trimleft.F90      wrt3dr.F90  zsmean.F90
  atod.F90        dry_adj.F90   gmean.F90      initp3d.F90   outfld.F90      set_eta.F90     upper.F90         xpavg.F90
  avgp2.F90       drymadj.F90   grads_ctl.F90  len_trim.F90  par_vecsum.F90  slp_das.F90     vmax4.F90         zflip.F90
  blowup.F90      epvd.F90      hswf.F90       linint.F90    pmaxmin.F90     time_shift.F90  vmax.F90          zmean.F90
  timingModule.F90 adtimingmodule.F90 g_timingmodule.F90
)

target_include_directories (fvgcm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
