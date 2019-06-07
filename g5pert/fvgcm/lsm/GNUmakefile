#
# Makefile for ESMA components.
#

ifeq ("$(USER)","suarez")
.SILENT:
endif

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/..
endif


# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include ../GMAO_fvbase.mk
  #--include ./Makefile.dependency

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

esma_help :
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE) "


THIS := fvgcm
LIB   = ../lib$(THIS).a

esma_install install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.h        $(ESMAINC)/$(THIS)

esma_clean clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd]

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]

esma_doc doc:
	@echo "Target $@ not implemented yet in `pwd`"


#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS :=\
allocbas.F  basindrv.F   cplavg.F     histend.F  histwrt.F  lsmdrv.F   lsmtvi.F   opnfil.F   restnam.F  soialb.F   suralb.F    vegconi.F\
allocht.F   basinini.F   diagtskin.F  histhan.F  infil.F    lsm.F      lsmzen.F   phenol.F   restwrt.F  soiconi.F  surbgc.F    vegtypi.F\
alloctc.F   basinmat.F   ecodyn.F     histini.F  lakconi.F  lsmhdr.F   map_i.F    phyconi.F  rst_lsm.F  soih2o.F   surphy.F    verh2o.F\
alloctv.F   bsslzr.F     fsurdat.F    histlst.F  laktem.F   lsmini.F   max_ovr.F  putfil.F   setarr.F   soilcol.F  surrad.F\
ao.F        canh2o.F     gaulat.F     histmlf.F  landtyp.F  lsmmain.F  memfree.F  radconi.F  slf2xy.F   soiltex.F  surtem.F\
ao_i.F      cell_area.F  getavu.F     histnam.F  landwat.F  lsmmap.F   memget.F   rank.F     snoalb.F   soitem.F   timedate.F\
area_ave.F  cell_box.F   getfil.F     histslf.F  lencs.F    lsmres.F   mlf2xy.F   rankn.F    snoconi.F  stomata.F  tridia.F\
atmini.F    cell_edge.F  histcrt.F    histtc.F   lsmctli.F  lsmtci.F   msspni.F   relavu.F   snoh2o.F   stoprun.F  twostr.F

OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS)))

# This include must be called here
# --------------------------------
include ../GMAO_fvgcm.mk

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

#----------------------------------------------------
# Hack to prevent remaking dep files during cleaning
#----------------------------------------------------
ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
   -include $(DEPS)
endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
#.
