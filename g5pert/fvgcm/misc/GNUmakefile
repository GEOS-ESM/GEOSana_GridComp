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
	$(CP) -p *.h            $(ESMAINC)/$(THIS)
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)

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
age_of_air.F90  diaginit.F90  findid.F90     init2dz.F90   merge.F90       ppme.F90        timingModule.F90  wrt3d.F90   zout.F90\
aoa_wrt.F90     diagout.F90   gmean4.F90     initmem.F90   minmax.F90      refout.F90      trimleft.F90      wrt3dr.F90  zsmean.F90\
atod.F90        dry_adj.F90   gmean.F90      initp3d.F90   outfld.F90      set_eta.F90     upper.F90         xpavg.F90\
avgp2.F90       drymadj.F90   grads_ctl.F90  len_trim.F90  par_vecsum.F90  slp_das.F90     vmax4.F90         zflip.F90\
blowup.F90      epvd.F90      hswf.F90       linint.F90    pmaxmin.F90     time_shift.F90  vmax.F90          zmean.F90\
timingModule.F90 adtimingmodule.F90 g_timingmodule.F90

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
