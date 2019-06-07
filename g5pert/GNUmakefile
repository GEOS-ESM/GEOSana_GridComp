#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 16Sep2004  da Silva  First ESMA version of fvGCM build.
# 30Mar2005  Todling   Adapted to work with module FVGSI.
#
#................................................................

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../..
endif


# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies

# Wire-in BASEDIR since g5pert works w/ different version of basedir
# than actual GEOS-5 (ESMF lib version issues)
# ------------------------------------------------------------------
ifeq ($(MACH), ia64)
  BASEDIR = /share/ESMA/baselibs/v1_8r1p
endif

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


#                  --------------------------------
#                   Recurse Make in Sub-directories
#                  --------------------------------

SUBDIRS = fvgcm svec ptdrvs 

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          help install clean doc

export ESMADIR BASEDIR ARCH SITE FVGCM_RES FVGCM_TRACER

$(TARGETS): 
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do			 \
	    ( cd $$d				;\
	      echo ""; echo Making $$t in `pwd`          ;\
	      $(MAKE) -e $$t ) \
	  done

allclean:
	$(MAKE) clean
	make -f Makefile clean

#                  --------------------
#                  User Defined Targets
#                  --------------------

local_install:
	@echo Nothing to locally install here

dist:
	$(MAKE) distclean
	( cd ..  ;\
          $(TAR)  cvf evac-`date '+%d%b%Y'`.tar src ;\
          $(GZIP) evac-`date '+%d%b%Y'`.tar )

distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]
	$(MAKE) esma_distclean

binclean:
	-$(RM) -r $(ESMADIR)/Config 
	-$(RM) -r $(ESMADIR)/$(ARCH)

realclean:
	$(MAKE) distclean
	$(MAKE) binclean

#              ------------------------------------------
#              Package Dependencies for Parallel Install
#              ------------------------------------------

svec_install   : fvgcm_install
ptdrvs_install : fvgcm_install svec_install

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
#.
