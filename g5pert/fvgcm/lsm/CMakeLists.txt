target_sources (fvgcm PRIVATE
  allocbas.F  basindrv.F   cplavg.F     histend.F  histwrt.F  lsmdrv.F   lsmtvi.F   opnfil.F   restnam.F  soialb.F   suralb.F    vegconi.F
  allocht.F   basinini.F   diagtskin.F  histhan.F  infil.F    lsm.F      lsmzen.F   phenol.F   restwrt.F  soiconi.F  surbgc.F    vegtypi.F
  alloctc.F   basinmat.F   ecodyn.F     histini.F  lakconi.F  lsmhdr.F   map_i.F    phyconi.F  rst_lsm.F  soih2o.F   surphy.F    verh2o.F
  alloctv.F   bsslzr.F     fsurdat.F    histlst.F  laktem.F   lsmini.F   max_ovr.F  putfil.F   setarr.F   soilcol.F  surrad.F
  ao.F        canh2o.F     gaulat.F     histmlf.F  landtyp.F  lsmmain.F  memfree.F  radconi.F  slf2xy.F   soiltex.F  surtem.F
  ao_i.F      cell_area.F  getavu.F     histnam.F  landwat.F  lsmmap.F   memget.F   rank.F     snoalb.F   soitem.F   timedate.F
  area_ave.F  cell_box.F   getfil.F     histslf.F  lencs.F    lsmres.F   mlf2xy.F   rankn.F    snoconi.F  stomata.F  tridia.F
  atmini.F    cell_edge.F  histcrt.F    histtc.F   lsmctli.F  lsmtci.F   msspni.F   relavu.F   snoh2o.F   stoprun.F  twostr.F
)

target_include_directories (fvgcm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

