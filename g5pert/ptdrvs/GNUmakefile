#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 16Feb2004  Zaslavsky  Initial imlementation.
# 20Oct2004  da Silva   Standardization
# 20Dec2006  Todling    Added fvsubset to install; compiling dyn2rs5 no big_endian
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

  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
     include $(ESMAINC)/fvgcm/GMAO_fvbase.mk # Local definitions
     include $(ESMAINC)/fvgcm/GMAO_fvgcm.mk  # fvgcm-specific
  endif
  #--include ./Makefile.dependency

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------


THIS := $(shell basename `pwd`)
LIB  = 

BINS = fvmodel.x fvpert.x jactest.x fvsens.x fvsvec.x qpert.x qsens.x

UTS = 

esma_install install: $(LIB) $(BINS)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p *.[Hh] $(ESMAINC)/$(THIS)
	$(CP) -p *.[Mm][Oo][Dd] $(ESMAINC)/$(THIS)
	$(CP) -p fvsvec.rc.sample  $(ESMAETC)/fvsvec.rc
	$(CP) -p jactest.rc.sample $(ESMAETC)/jactest.rc
	$(CP) -p pseudo.rc.sample  $(ESMAETC)/pseudo.rc
	$(CP) -p tracer.rc.sample  $(ESMAETC)/tracer.rc
	$(CP) -p fvgcm.ccmflags.namelist $(ESMAETC)
	$(CP) -p *.tmpl  $(ESMAETC)
	$(CP) -p $(BINS) $(ESMABIN)

esma_clean clean:
	$(RM) *~ *.[aox] *.[Mm][Oo][Dd]
	$(RM) -r *.hdf *.grd *.[axol]  \
              *~ *.log *.ps *.dvi *.aux *.toc *.lot *.lof rii_files/

esma_distclean distclean:
	$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]
	$(RM) -r *.hdf *.grd *.[axol]  \
              *~ *.log *.ps *.dvi *.aux *.toc *.lot *.lof rii_files/
esma_help:
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
	@echo "       MPFLAG = $(MPFLAG)"

help:
	@ echo "Available targets:"
	@ echo "  make lib             (to create libods.a)"
	@ echo "  make lib MODOBJS=    (to omit f90 modules)"
	@ echo "  make test            (to test ODS library)"
	@ echo "  make doc             (to create documentation)"
	@ echo "  make all"  
	@ echo "  make export"  
	@ echo "  make clean"  
	@ echo "  make distclean"


#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS = m_fvpsas.F90

OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS))) $(addsuffix .d, $(basename $(BINS)))

MOD_DIRS   += $(INC_SVEC)
USER_FFLAGS = $(XUSER_FFLAGS) $(MOD_DIRS)

LIB_UTIL  = $(ESMALIB)/libutil.a
LIB_PUTIL = $(ESMALIB)/libputil.a

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)


# test and auxiliary programs
# ---------------------------

LIB_SYS += -L/opt/pd/gcc/gcc-3.3.1/lib -lstdc++

%.x : %.o $(OBJS) $(LIB_SVEC) $(LIB_FVGCM) 
	$(LD) $(LDFLAGS)  -o $@ $*.o $(OBJS) \
	      $(LIB_SVEC) $(LIB_ARPACK) $(LIB_PARPACK) $(LIB_PLANSO) $(LIB_LANSO) \
	      $(LIB_SVEC) $(LIB_FVGCM) $(LIB_ARPACK) $(LIB_PARPACK) $(LIB_HERMES) \
	      $(LIB_GFIO) $(LIB_UTIL) $(LIB_PUTIL) $(LIB_MPEU) $(LIB_SDF) $(LIB_MPI) $(LIB_SCI) $(LIB_SYS)

#----------------------------------------------------
# Hack to prevent remaking dep files during cleaning
#----------------------------------------------------
ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
   -include $(DEPS)
endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
#.
