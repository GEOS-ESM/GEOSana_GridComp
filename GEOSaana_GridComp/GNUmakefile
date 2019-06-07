#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 09Jun2003  da Silva  First crack.
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

esma_help help:
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
	@echo "         SITE = $(SITE)"

THIS := $(shell basename `pwd`)
LIB   = lib$(THIS).a

ALLDIRS = GSI_GridComp GEOSgsi_Coupler

SUBDIRS = $(wildcard $(ALLDIRS))

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          install clean doc

export ESMADIR BASEDIR ARCH SITE

$(TARGETS):
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do                    \
	    ( cd $$d                            ;\
	      echo ""; echo Making $$t in `pwd`          ;\
	      $(MAKE) -e $$t ) \
	  done
	$(MAKE) local_$@

local_esma_install local_install : $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a            $(ESMALIB)
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)

local_esma_clean local_clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd]

local_esma_distclean local_distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]

#                  --------------------
#                  User Defined Targets
#                  --------------------

INC_GEOS = $(INC_MAPL_BASE) $(INC_GEOS_SHARED)

SRCS := GEOSaana_GridCompMod.F90
OBJS := $(addsuffix .o, $(basename $(SRCS))) 
DEPS := $(addsuffix .d, $(basename $(SRCS))) 

THIS_SP    = NCEP_sp_r8i4
THIS_W3    = NCEP_w3_r8i4
THIS_BACIO = NCEP_bacio_r4i4
THIS_BUFR  = NCEP_bufr_r8i4
THIS_CFIO  = GMAO_cfio_r4

FREAL = $(FREAL4)

INC_DIRS = . $(ESMAINC)/MAPL_Base $(ESMAINC)/GSI_GridComp $(ESMAINC)/GMAO_mpeu $(ESMAINC)/MAPL_cfio

MOD_DIRS = . $(INC_SDF) $(INC_ESMF) $(INC_CFIO) \
	     $(ESMAINC)/GEOS_Shared \
	     $(ESMAINC)/MAPL_Base \
	     $(ESMAINC)/GSI_GridComp \
	     $(INC_GFIO) $(INC_CFIO) $(INC_MPEU) \
	     $(INC_HERMES) $(INC_CRTM) $(INC_IRSSE) \
	     $(INC_SFCIO) $(INC_TRANSF) $(INC_IRUTIL) $(INC_RADTRANS) \
	     $(foreach dir,$(SUBDIRS),$(ESMAINC)/$(dir))

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FFLAGS = $(M)$(ESMAINC)/$(THIS_CFIO) $(M)$(ESMAINC)/GSI_GridComp

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS) 
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#              ------------------------------------------
#              Package Dependencies for Parallel Install
#              ------------------------------------------
  GEOSgsi_Coupler_install : GSI_GridComp_install

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
#.
