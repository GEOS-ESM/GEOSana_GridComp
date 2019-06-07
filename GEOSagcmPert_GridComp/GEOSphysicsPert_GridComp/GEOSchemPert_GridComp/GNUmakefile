#
# Makefile for ESMA components.
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/..
endif


# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # GMAO Generic stuff

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


THIS := $(shell basename `pwd`)
LIB   = $(LIB_THIS)
#                  --------------------------------
#                   Recurse Make in Sub-directories
#                  --------------------------------

ALLDIRS = GEOSpchemPert_GridComp GOCART_PertGridComp \
          StratChem_PertGridComp GMIchem_PertGridComp

SUBDIRS = $(wildcard $(ALLDIRS))

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          install clean distclean doc 

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

export ESMADIR BASEDIR ARCH SITE FREAL

$(TARGETS): 
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do			 \
	    ( cd $$d				;\
	      echo ""; echo Making $$t in `pwd`          ;\
	      $(MAKE) -e $$t ) \
	  done
	$(MAKE) local_$@

local_esma_install local_install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.mod          $(ESMAINC)/$(THIS)

local_esma_clean local_clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd]

local_esma_distclean local_distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd] *___.[fF]*

local_esma_doc local_doc:
	@echo "Target $@ not implemented yet in `pwd`"

#                  --------------------
#                  User Defined Targets
#                  --------------------

# If childrem not present, create stubs
# -------------------------------------
ifeq ( $(wildcard GMIchem_GridComp), $(null) )
        STUBS += GMIchem_GridCompMod___.F90
endif
ifeq ( $(wildcard GOCART_GridComp), $(null) )
        STUBS += GOCART_GridCompMod___.F90
endif
ifeq ( $(wildcard StratChem_GridComp), $(null) )
        STUBS += StratChem_GridCompMod___.F90
endif

#SRCS := GEOS_ChemGridComp.F90 GEOS_ChemEnvGridComp.F90 $(STUBS)
SRCS := GEOS_ChemPertGridComp.F90
OBJS := $(addsuffix .o, $(basename $(SRCS))) 
DEPS := $(addsuffix .d, $(basename $(SRCS))) 

NEED_DIRS = . $(INC_GMAO_SHARED) $(INC_ESMF) i#$(INC_GEOS_CHEM)
INC_DIRS = . $(NEED_DIRS) $(foreach dir,$(SUBDIRS),$(ESMAINC)/$(dir))
MOD_DIRS =   $(INC_DIRS)


USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FFLAGS = $(BIG_ENDIAN) 

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

# Stubs
# -----
%___.F90: $(STUB) 
	$(STUB) $* > $@ 

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.











  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

