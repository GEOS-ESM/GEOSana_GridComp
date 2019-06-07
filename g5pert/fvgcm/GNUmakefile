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
LIB   = lib$(THIS).a
FVGCMRES = b72
FVGCMTRC = 1

#                  --------------------------------
#                   Recurse Make in Sub-directories
#                  --------------------------------

ALLDIRS = \
csm_share atmlnd_share pilgrim ecmfft phys lsm misc drvs fvcore fvtlm fvadm fvptrace

SUBDIRS = $(wildcard $(ALLDIRS))

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          install clean doc 

export ESMADIR BASEDIR ARCH SITE

$(TARGETS): fvgcm.h
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do			 \
	    ( cd $$d				;\
	      echo ""; echo Making $$t in `pwd`          ;\
	      $(MAKE) -e $$t ) \
	  done
	$(MAKE) local_$@

local_esma_install local_install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a        $(ESMALIB)
	$(CP) -p *.h        $(ESMAINC)/$(THIS)
	$(CP) -p GMAO_fvbase.mk  $(ESMAINC)/$(THIS)
	$(CP) -p GMAO_fvgcm.mk   $(ESMAINC)/$(THIS)

local_esma_clean local_clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd] fvgcm.h

local_esma_distclean local_distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd] fvgcm.h

local_esma_doc local_doc:
	@echo "Target $@ not implemented yet in `pwd`"

fvgcm.h:
	./fvsetdim -t $(FVGCMTRC) $(FVGCMRES)

#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS := 
OBJS := $(addsuffix .o, $(basename $(SRCS))) 
DEPS := $(addsuffix .d, $(basename $(SRCS))) 
FDP   = 
AR_FLAGS += a

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FFLAGS = $(BIG_ENDIAN) 

vpath % $(MOD_DIRS)

$(LIB) lib : fvgcm.h $(OBJS) 
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)


#
# These includes much remain here given the funny
# dependence of stepon_ad/tl on m_trjmng and the 
# position of this latter in the fvcore
#
#--include fvcore/Makefile.dependency
#--include  fvtlm/Makefile.dependency
#--include  fvadm/Makefile.dependency
#--include fvptrace/Makefile.dependency

#              ------------------------------------------
#              Package Dependencies for Parallel Install
#              ------------------------------------------
#
# NOTE: (J.Stassi 3/09/2019) Build SUBDIRS directories one at a time, because
#       they all add object files to libfvgcm.a, and if two directories try
#       to add files at the same time, then files from one may not get there.

  atmlnd_share_install : csm_share_install
       pilgrim_install : atmlnd_share_install
        ecmfft_install : pilgrim_install
          phys_install : ecmfft_install fvgcm.h
           lsm_install : phys_install fvgcm.h
          misc_install : lsm_install
          drvs_install : misc_install
        fvcore_install : drvs_install
         fvtlm_install : fvcore_install
         fvadm_install : fvtlm_install
      fvptrace_install : fvadm_install

#          clm2_install : fvgcm.h
#          util_install : misc_install fvgcm.h

#----------------------------------------------------
# Hack to prevent remaking dep files during cleaning
#----------------------------------------------------
ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
   -include $(DEPS)
endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
#.
