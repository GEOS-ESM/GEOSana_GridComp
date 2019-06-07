# By R. Todling, Jan 2007

# Define conversion of source to object
# and deactivate GNUmake auto-dependency built
# --------------------------------------------
#DEPS  =
#FDP   =
#ifeq ($(MACH), x86_64)
#     INC_MPI := $(filter /opt/scali%,$(FPATHS))
#     MYMPT   := $(subst include,lib,$(INC_MPI))
#else
##             /opt/sgi/mpt/1.9.1.0/include
#     MYMPT  = /opt/sgi/mpt/1.12.0.nas
#     LIB_MPI = -L$(MYMPT)/lib -lmpi
#     LIB_SYS = -lscs
#endif


# Define location of includes/modules
# MPT is wired in for now ... this will be generalization to platforms beyond PALM
# --------------------------------------------------------------------------------
MYMOD_DIRS = -I../ -I$(INC_HERMES) -I$(INC_TRANSF) -I$(INC_GFIO) -I$(INC_CFIO) -I$(INC_MPEU) \
	           -I$(INC_SDF) -I$(INC_MAPL_BASE) -I$(INC_GEOS_SHARED) \
                   -I$(ESMADIR)/src/g5pert/fvgcm/fvgcm \
                   -I$(ESMADIR)/src/g5pert/fvgcm/atmlnd_share \
                   -I$(ESMADIR)/src/g5pert/fvgcm/csm_share \
                   -I$(ESMADIR)/src/g5pert/fvgcm/pilgrim \
                   -I$(ESMADIR)/src/g5pert/fvgcm/misc \
                   -I$(ESMADIR)/src/g5pert/fvgcm/fvadm \
                   -I$(ESMADIR)/src/g5pert/fvgcm/fvcore \
                   -I$(ESMADIR)/src/g5pert/fvgcm/fvtlm \
                   -I$(ESMADIR)/src/g5pert/fvgcm/ecmfft \
                   -I$(ESMADIR)/src/g5pert/fvgcm/drvs \
                   -I$(ESMADIR)/src/g5pert/fvgcm/phys \
                   -I$(ESMADIR)/src/g5pert/fvgcm/fvtrace \
                   -I$(ESMADIR)/src/g5pert/fvgcm/lsm \
                   -I$(BASEDIR)/$(ARCH)/include/netcdf \
                   -I$(INC_MPI)


MYINC_DIRS = $(MYMOD_DIRS)
USER_FINCS := $(MYMOD_DIRS)

MOD_DIRS = $(INC_HERMES) $(INC_TRANSF) $(INC_GFIO) $(INC_CFIO) $(INC_MPEU) \
	   $(INC_SDF) $(INC_MAPL_BASE) $(INC_GEOS_SHARED) \
	   $(INC_FVGCM) \
           $(INC_MPI)


# Define compilation specific instructions
# ----------------------------------------
# OpenMP Only
# -----------
#DEFS = -DHIDE_SHR_MSG -DLINUX -DTIMING -DREAL8 -DLSMH_off -DVERSION="'fvgcm'" -DHIDE_MPI -DUSE_OPENMP -DSET_CPUS -DGFIO -DCHECKPOINTING -DLinux
#XUSER_FFLAGS = -extend_source -r8 -w -cm -O3 -cpp -openmp -convert big_endian $(DEFS) 
#LDFLAGS += -openmp

# MPI Only 
# --------
USER_FDEFS = -DTIMING  -DREAL8 -DLSMH_off -DVERSION="'fvgcm'"  -DSPMD -DGFIO -DCHECKPOINTING -DLINUX
USER_FDEFS =           -DREAL8 -DLSMH_off -DVERSION="'fvgcm'"  -DSPMD -DGFIO -DCHECKPOINTING -DLINUX
USER_FDEFS =           -DREAL8 -DLSMH_off -DSPMD -DGFIO -DCHECKPOINTING -DLINUX
XUSER_FFLAGS = -extend_source 180 -r8 -w -cm -O3 -cpp -mp -convert big_endian $(USER_FDEFS)

# Hybrid: MPI+OpenMP
# ------------------
#DEFS =           -DREAL8 -DLSMH_off -DVERSION="'fvgcm'"  -DSPMD -DGFIO -DCHECKPOINTING -DUSE_OPENMP -DLINUX -DLinux -DSET_CPUS
#XUSER_FFLAGS = -extend_source -r8 -w -cm -O3 -cpp -convert big_endian $(DEFS) -openmp
#LDFLAGS += -openmp

FPE  =
USER_FFLAGS = $(XUSER_FFLAGS) $(MYMOD_DIRS)
USER_CFLAGS = $(I). $(I)$(INC_SDF) -DLINUX -DFORTRANUNDERSCORE $(USER_FDEFS)
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 

LIB_MPEU = $(ESMALIB)/libGMAO_eu.a
LIB_MPEU = $(ESMALIB)/libGMAO_mpeu.a
