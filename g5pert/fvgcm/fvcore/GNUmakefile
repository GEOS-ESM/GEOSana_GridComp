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
  include ../GMAO_fvbase.mk               # Local definitions
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
	$(CP) -p *.mod      $(ESMAINC)/$(THIS)

esma_clean clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd]

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]

esma_doc doc:
	@echo "Target $@ not implemented yet in `pwd`"


#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS := \
avgc.F90     controlq.F90   dependentq.F90   func.F90    initfunc.F90     pkez.F90                  sw_core.F90  vpol5.F90\
benergy.F90  d2a2.F90       d_split.F90      fvcore.F90  mapz_module.F90  nlmodel.F90  polavg.F90       setfunc.F90   tp_core.F90  xxall.F90\
cd_core.F90  d2a3d.F90      fill_module.F90  geopk.F90   m_iostate.F      pft2d.F90    postfunc.F90     setrig.F90     trac2d.F90\
control.F90  dependent.F90  fillz.F90        highp2.F90  m_pftSV.F90      pft_cf.F90   prognostics.F90  stepon.F       upol5.F90 \
m_trajmng.F m_trjphys.F90   a2d3d.F90        m_delp2ps.F90  prognostics_q.F90

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
