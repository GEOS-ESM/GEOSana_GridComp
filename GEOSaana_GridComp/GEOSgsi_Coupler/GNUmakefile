#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 3mar2004  Zaslavsky  Initial imlementation.
# 20Oct2004  da Silva  Standardization
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------


THIS = $(shell basename `pwd`)
LIBO = lib$(THIS).a
LIBX = libGEOSgsi_Util.a
LIB  = $(LIBO) $(LIBX)
BIN  =

test:
	: test  = "$@"
	: SHELL = $(SHELL)
	: PWD   = $(shell pwd)
	: THIS  = $(THIS)
	: LIB   = $(LIB)
	: MODD	= $(MOD_DIRS)

# To de/activate GEOS_PERT-related routines
# -----------------------------------------
GEOS_PERT = -DGEOS_PERT				# not clear what it means in general
apert_presence  := $(shell ls $(ESMALIB)/libGEOSdynamicsPert_GridComp.a)
g5pert_presence := $(shell ls $(ESMALIB)/libfvgcm.a $(ESMALIB)/libsvec.a)

libpresence := $(apert_presence)		# default is apert, if present
ifeq ( $(wildcard $(libpresence)),$(null) )
   libpresence := $(g5pert_presence)		# or g5pert, if present, see below
   ifeq ( $(wildcard $(libpresence)),$(null) )
     GEOS_PERT += -DSTUB_AGCMPERT		# or use stub_apert, at least
   else
     GEOS_PERT += -DG5_PERTMOD			# set to use g5_pertmod, also
     GEOS_PERT += -DSAVE_DYN_PROG_STATE		# compile it with "#ifdef SAVE_DYN_PROG_STATE"
   endif
endif

esma_install install: $(LIB) $(BIN)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p *.a         $(ESMALIB)
	$(CP) -p *.mod       $(ESMAINC)/$(THIS)

esma_clean clean:
	$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]

esma_distclean distclean:
	$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]

esma_doc doc:
	@echo "Target $@ not implemented yet in `pwd`"


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
	@echo "        FREAL = $(FREAL)"
	@echo "    MODS_DIRS = $(MOD_DIRS)"
	@echo "   MMODS_DIRS = $(MMOD_DIRS)"
	@echo "        VPATH = $(VPATH)"

#                  --------------------
#                  User Defined Targets
#                  --------------------


SRCSO =	geos_pertState.F90 \
	geos_pertStateIO.F90 \
	geos_StateIO.F90 \
        geos_pertmod.F90 \
        geos_nstmod.F90 \
        geos_gsi_util.F90 \
	g5_pertmod.F90   \
	cplr_pertmod.F90 \
	cplr_ensemble.F90 \
	cplr_nst.F90 \
        cplr_timermod.F90 \
        stub_GEOSagcmPert_GridCompMod.F90

#       geos_pgcmtest.F90

SRCSX =	crtm_aerosol.F90 cplr_Set_CRTM_Aerosol.F90 \
        crtm_cloud.F90   cplr_Set_CRTM_Cloud.F90

SRCS = $(SRCSX) $(SRCSO)

##        gsi_4dcouplermod.F90 \
##        timermod.F90

OBJSX := $(addsuffix .o, $(basename $(SRCSX)))
OBJS  := $(addsuffix .o, $(basename $(SRCS)))
DEPS  := $(addsuffix .d, $(basename $(SRCS)))

_D = $(GEOS_PERT)
#_D += -DDEBUG_TRACE -DDEBUG_MPOUT	# switch on DEBUG_TRACE for tracing and
					# DEBUG_MPOUT for parallel stdout streams

FOPT       = $(FOPT3)
FREAL      = $(FREAL4) 
FPE        =

INC_APERT = $(ESMAINC)/GEOSagcmPert_GridComp $(ESMAINC)/GEOSdynamicsPert_GridComp $(ESMAINC)/fvdycorepert \
            $(ESMAINC)/GEOS_PertShared $(ESMAINC)/Tapenade

MOD_DIRS = . $(INC_MPEU) $(INC_HERMES) $(INC_APERT) $(INC_TRANSF) $(INC_GEOS) $(INC_FVGCM) \
             $(INC_GSI) $(INC_CRTM) $(INC_CHEM_BASE) \
             $(INC_MAPL_BASE) $(INC_CFIO) $(INC_ESMF) $(INC_SDF)
MMOD_DIRS = $(M).  $(M)$(INC_MPEU) $(M)$(INC_HERMES) $(M)$(INC_APERT) $(M)$(INC_TRANSF) $(M)$(INC_GEOS) \
                   $(M)$(INC_FVGCM) $(M)$(INC_GSI) $(M)$(INC_CRTM) $(M)$(INC_CHEM_BASE) \
                   $(M)$(INC_MAPL_BASE) $(M)$(INC_CFIO) $(M)$(INC_ESMF) $(M)$(INC_SDF)
USER_FDEFS = $(_D)
USER_FFLAGS = $(BIG_ENDIAN) $(BYTERECLEN) -traceback
USER_CFLAGS = -I . -Dfunder -DFortranByte=char -DFortranInt=int -DFortranLlong='long long' -O3
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FINCS = $(foreach dir,$(MOD_DIRS), $(I)$(dir))

vpath % $(MOD_DIRS) /usr/include

lib : $(LIBX) $(LIBO)

$(LIBX) : $(LIB_GSI_UTIL) $(OBJSX)
	$(CP) $(LIB_GSI_UTIL) $(LIBX)
	$(AR)         -d $(LIBX) stub_set_crtm_aerosol.o stub_set_crtm_cloud.o
	$(AR) $(ARFLAGS) $(LIBX) $(OBJSX)

$(LIBO) : $(LIB_GSI_UTIL) $(OBJS)
	$(CP) $(LIB_GSI_UTIL) $(LIBO)
	$(AR)         -d $(LIBO) stub_pertmod.o stub_ensmod.o stub_nstmod.o stub_timermod.o stub_set_crtm_aerosol.o stub_set_crtm_cloud.o
	$(AR) $(ARFLAGS) $(LIBO) $(OBJS)

%.x : $(LIBO) %.o 
	$(LD) $(LDFLAGS) -o $@ $*.o $(LIBO) $(LIB_SYS)


# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

m_mpif.mod: $(LIB_MPEU)

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
#.
