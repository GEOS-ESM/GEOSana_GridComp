target_sources (fvgcm PRIVATE
  addfld.F    cldclw.F    diagref.F  intmax.F     mgwintr.F          prealc.F  radinp.F    sstan.F         trcplk.F      whenfle.F
  aermix.F    cldefr.F    esat.F     iostop.F     mgworo.F           q1q2.F    radoz2.F    sstini.F        trcpth.F      whenflt.F
  albocean.F  cldems.F    esinti.F   isrchfgt.F   mgwprof.F          qmassa.F  radozn.F    sstint.F        trefoce.F     whenne.F
  aqsatd.F    cldfrc.F    fldlst.F   isrchfle.F   mtdlss.F           qneg3.F   radtpl.F    subrturb.F      trefsice.F    writeric.F
  aqsat.F     cldint.F    flxoce.F   kdpfnd.F     mvdiff.F           qneg4.F   rdharr.F    timing_stubs.F  tropopause.F  zenith.F
  attachnc.F  cldprp.F    flxsice.F  lcbas.F      navu.F             qsat.F    rdhdr.F     torgrid.F       tsinti.F      zint.F
  binf2c.F    cldsav.F    fmrgrid.F  lcdbas.F     nrl_gwd.F          qvdiff.F  readric.F   tphysaa.F       turb.F
  bldcom.F    closure.F   freemem.F  lenchr.F     ocez0.F            radabs.F  relunit.F   tphysac.F       vdiff.F
  bldfld.F    cmfmca.F    gestbl.F   lininterp.F  orbit.F            radclr.F  resetr.F    tphysbc.F       vdinti.F
  blk_dat.F   cond.F      getmem.F   ljust.F      oznini.F           radclw.F  sas_drv.F   trajmp.F        vdintr.F
  bnddyi.F    conv_ccm.F  getunit.F  m_diagsfc.F  oznint.F           radcsw.F  sas.F       trbflx.F        virtem.F
  buoyan.F    convtran.F  gffgch.F   mfinti.F     pbinti.F           radctl.F  smallest.F  trcab.F         vqsatd.F
  cainti.F    copyac.F    igtseq.F   mgwbgnd.F    pbldif.F           radded.F  srfoce.F    trcabn.F        vqsat.F
  ccmoce.F    dadadj.F    inifld.F   mgwdrag.F    pbldif_simpl.F     radems.F  srfsice.F   trcems.F        wheneq.F
  ccmslp.F    diagbma.F   initext.F  mgwinti.F    pbldif_simpl_tq.F  radini.F  srftsb.F    trcmix.F        whenfgt.F
  m_pmgrid.F  m_comsrf.F
)

target_include_directories (fvgcm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
