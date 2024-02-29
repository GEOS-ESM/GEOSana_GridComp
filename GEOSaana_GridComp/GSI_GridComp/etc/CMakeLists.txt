file(GLOB perlscripts *.pl)
file(GLOB pythonscripts *.py)

install(
   PROGRAMS ${perscripts} ${pythonscripts}
   DESTINATION bin)

set (etc_txt_files
   gmao_global_aeroinfo.txt
   gmao_global_anavinfo.txt
   gmao_global_anavinfo_l91.txt
   gmao_global_anavinfo_l132.txt
   gmao_global_anavinfo_l181.txt
   gmao_global_anavinfo_rcov.txt
   gmao_global_blacklist.txt
   gmao_global_cloudy_radiance_info.txt
   gmao_global_tgasinfo.txt
   gmao_global_convinfo.txt
   gmao_global_insituinfo.txt
   gmao_global_ozinfo.txt
   gmao_global_pcpinfo.txt
   gmao_global_rad_time_thin.txt
   gmao_global_satinfo.txt
   gmao_global_scaninfo.txt
   gmao_global_tlapmean.txt
   gsi_atms_beamwidth.txt
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
   FILES gmao_global_hybens_info.x288y181l72.txt
   DESTINATION etc
   RENAME gmao_global_hybens_info.x288y181l72.rc
   )
install(
   FILES gmao_global_hybens_info.x576y361l72.txt
   DESTINATION etc
   RENAME gmao_global_hybens_info.x576y361l72.rc
   )

install(
   FILES gmao_global_hybens_info.x288y181l91.txt
   DESTINATION etc
   RENAME gmao_global_hybens_info.x288y181l91.rc
   )
install(
   FILES gmao_global_hybens_info.x576y361l91.txt
   DESTINATION etc
   RENAME gmao_global_hybens_info.x576y361l91.rc
   )

install(
   FILES gmao_global_hybens_info.x288y181l132.txt
   DESTINATION etc
   RENAME gmao_global_hybens_info.x288y181l132.rc
   )
install(
   FILES gmao_global_hybens_info.x576y361l132.txt
   DESTINATION etc
   RENAME gmao_global_hybens_info.x576y361l132.rc
   )

install(
   FILES gmao_global_hybens_info.x288y181l181.txt
   DESTINATION etc
   RENAME gmao_global_hybens_info.x288y181l181.rc
   )
install(
   FILES gmao_global_hybens_info.x576y361l181.txt
   DESTINATION etc
   RENAME gmao_global_hybens_info.x576y361l181.rc
   )

set (etc_other_files
   gmao_airs_bufr.tbl
   gsi_fdda_1.rc.tmpl
   gsi_fdda_2.rc.tmpl
   gsi_fgat_1.rc.tmpl
   gsi_fgat_2.rc.tmpl
   gsi.rc.tmpl
   gsi_sens.rc.tmpl
   gsidiags.rc
   obs.rc.tmpl
   )
install (
   FILES ${etc_other_files}
   DESTINATION etc
   )

# We can't install this with DIRECTORY as in CVS
# only specific files were copied to install/etc
set (MERRA2_files
   MERRA2/gmao_global_satinfo.rc
   MERRA2/gmao_global_scaninfo.rc
   MERRA2/gsi.rc.tmpl
   MERRA2/gsi_sens.rc.tmpl
   MERRA2/sac.nl.tmpl
   )
install(
   FILES ${MERRA2_files}
   DESTINATION etc/gsi/MERRA2
   )

# Handle O3-only analysis settings
set (o3anasa_files
   o3anasa/gmao_global_anavinfo_rcov.rc
   o3anasa/gmao_global_anavinfo.rc
   o3anasa/gsi.rc.tmpl
   o3anasa/GSI_GridComp.rc.tmpl
   )
install(
   FILES ${o3anasa_files}
   DESTINATION etc/gsi/o3anasa
   )

# Handle TGAS-only analysis settings
set (coanasa_files
   coanasa/gmao_global_anavinfo_rcov.rc
   coanasa/gmao_global_anavinfo.rc
   coanasa/gsi.rc.tmpl
   coanasa/GSI_GridComp.rc.tmpl
   )
install(
   FILES ${coanasa_files}
   DESTINATION etc/gsi/coanasa
   )

# Handle GEOSIT settings in a similar way as MERRA2
set (GEOSIT_files
   GEOSIT/gmao_global_cloudy_radiance_info.rc
   GEOSIT/gmao_global_satinfo.rc
   GEOSIT/gsi.rc.tmpl
   GEOSIT/gsi_sens.rc.tmpl
   )
install(
   FILES ${GEOSIT_files}
   DESTINATION etc/gsi/GEOSIT
   )

# Handle GEOS-R21C settings in a similar way as MERRA2
set (R21C_files
   R21C/gmao_global_cloudy_radiance_info.rc
   R21C/gmao_global_convinfo.rc
   R21C/gmao_global_satinfo.rc
   R21C/gmao_global_tlapmean.rc
   R21C/gmao_global_ozinfo.rc
   R21C/gsi.rc.tmpl
   R21C/gsi_sens.rc.tmpl
   R21C/prepobs_errtable.global
   )
install(
   FILES ${R21C_files}
   DESTINATION etc/gsi/R21C
   )
file(GLOB etc_py *.py)
install (
   PROGRAMS ${etc_py}
   DESTINATION bin
   )


