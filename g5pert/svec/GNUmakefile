#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 16Feb2004  Zaslavsky  Initial imlementation.
# 20Oct2004  da Silva   Standardization
# 20Dec2006  Todling    Added fvsubset to install; compiling dyn2rs5 no big_endian
#

ifeq ("$(USER)","suarez")
.SILENT:
endif

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
LIB   = lib$(THIS).a

BINS = initadj.x fsens2pert.x simsv.x postStats.x svspectra.x \
       pertenergy.x pseudo.x

UTS = 

esma_install install: $(LIB) $(BINS)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)
	$(CP) -p *.a    $(ESMALIB)
	$(CP) -p *.[Hh] $(ESMAINC)/$(THIS)
	$(CP) -p *.[Mm][Oo][Dd] $(ESMAINC)/$(THIS)
	$(CP) simsv.rc.sample     $(ESMAETC)/simsv.rc
	$(CP) initadj.rc.sample   $(ESMAETC)/initadj.rc
	$(CP) postStats.rc.sample $(ESMAETC)/postStats.rc
	$(CP) pertenergy.rc.tmpl  $(ESMAETC)/pertenergy.rc.tmpl
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

SRCS := \
      cnop2spg.F90 \
      m_admtlm.F90 \
      m_cnop.F90\
      m_eigen.f90\
      m_fsens2pert.F90\
      m_gAtlm.F90\
      m_gAadm.F90\
      m_gQadm.F90\
      m_gQtlm.F90\
      m_fdiff.F90\
      m_initadj.F90\
      m_jacobian.F90\
      m_postStats.F90\
      m_shtrans.f90\
      m_svnorms.F90 \
      m_getModConst.F\
      m_pertenergy.F90\
      m_pseudo.F90\
      m_simsv.F90\
      m_svprojs.F90 \
      m_poles.F90\
      m_setsvecs.F\
      m_sVecDef.F\
      m_svspectra.F90 \
      m_pertutil.F90 \
      op.F90\
      pl_vert_sp.f \
      pl_spec_sp.f \
      svecdrv.F90

OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS))) $(addsuffix .d, $(basename $(BINS)))

USER_FFLAGS = $(XUSER_FFLAGS) $(MOD_DIRS) $(FOPTS)
DEFS += -D_DLAPACK_  -D_DARPACK_

MOD_DIRS += $(INC_ARPACK)
LIB_UTIL  = $(ESMALIB)/libutil.a
LIB_PUTIL = $(ESMALIB)/libputil.a

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)


# test and auxiliary programs
# ---------------------------

%.x : $(LIB) %.o 
	$(LD) $(LDFLAGS)  -o $@ $*.o \
	      $(LIB) $(LIB_ARPACK) $(LIB_PARPACK) $(LIB_PLANSO) $(LIB_LANSO) \
	      $(LIB) $(LIB_FVGCM) $(LIB_ARPACK) $(LIB_PARPACK) $(LIB_TRANSF) $(LIB_MPEU) $(LIB_HERMES) \
	      $(LIB_GFIO) $(LIB_UTIL) $(LIB_PUTIL) $(LIB_SDF) $(LIB_MPI) $(LIB_SCI) $(LIB_SYS)

#----------------------------------------------------
# Hack to prevent remaking dep files during cleaning
#----------------------------------------------------
ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
   -include $(DEPS)
endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
#.
