#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 03Oct2011  Todling    Adapted from GMAO_etc (Stassi's style)
#

# Get GNUmakefile name (requires gmake versions >= 3.80)
#-------------------------------------------------------
ifdef MAKEFILE_LIST
   GNUMAKEFILE := $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))
else
   GNUMAKEFILE := GNUmakefile
endif

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../../../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

########################################################################

ALLPERLSCRIPT = 

PERLSCRIPT = $(wildcard $(ALLPERLSCRIPT))

PERLPKG = $(wildcard *.pm)

#-------------------------------------------------------
# ALLSRC   = source files for binaries with no extension
# ALLSRCX  = source files for binaries with .x extension
# OTHERBIN = other files for ESMABIN directory
# ALLETC   = files for ESMAETC directory
#-------------------------------------------------------
ALLSRC = 

ALLSRCX = 

OTHERBIN = 

ALLETC = gmao_airs_bufr.tbl \
         gmao_global_aeroinfo.txt \
         gmao_global_anavinfo.txt \
         gmao_global_blacklist.txt \
         gmao_global_tgasinfo.txt \
         gmao_global_convinfo.txt \
         gmao_global_hybens_betainfo.l72.txt \
         gmao_global_hybens_locinfo.l72.txt \
         gmao_global_insituinfo.txt \
         gmao_global_ozinfo.txt \
         gmao_global_pcpinfo.txt \
         gmao_global_satinfo.txt \
         gmao_global_scaninfo.txt \
         gsi_atms_beamwidth.txt \
         gsi_fdda_1.rc.tmpl \
         gsi_fdda_2.rc.tmpl \
         gsi_fgat_1.rc.tmpl \
         gsi_fgat_2.rc.tmpl \
         gsi.rc.tmpl \
         gsi_sens.rc.tmpl \
	 gsidiags.rc \
         obs.rc.tmpl
 
########################################################################

bin        = $(BINS)
etc        = $(ETCFILES)
txt        = $(TXTFILES)
perlpkg    = $(PERLPKG)
perlscript = $(PERLSCRIPT)
PERLFILES  = $(PERLSCRIPT) $(PERLPKG)

SEDCMD = s^\@DASPERL^\$(PERL)^;s^\@DASSED^\$(SED)^;s^\@DASHELL^\$(SHELL)^

########################################################################

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

THIS := $(shell basename $(PWD))

SRC = # $(wildcard $(ALLSRC))
BIN = # $(basename $(SRC))

SRCX = # $(wildcard $(ALLSRCX))
BINX = # $(addsuffix .x, $(basename $(SRCX)))

BINO = $(wildcard $(OTHERBIN))

ALLDIRS = MERRA2

# all files for bin, etc, and txt
#---------------------------
BINS      = $(BIN) $(BINX) $(BINO)
ETCFILES  = $(wildcard $(ALLETC))
SUBDIRS   = $(wildcard $(ALLDIRS))

TARGETS = esma_install esma_clean esma_distclean esma_doc install clean distclean doc
export ESMADIR BASEDIR ARCH SITE

$(TARGETS):
	@ t=$@; argv="$(SUBDIRS)"  ;\
	  for d in $$argv; do  \
	    ( cd $$d  ;\
	      echo ""; echo Making $$t in `pwd`  ;\
	      $(MAKE) -e $$t )  \
	  done
	$(MAKE) local_$@

local_esma_install local_install: $(BINS)
	$(MKDIR) $(ESMAETC) $(ESMABIN) 

	@if [ "$(BINS)" != "  " ]; then (\
	   gmake -f \$(GNUMAKEFILE) bin \
	) fi

	@if [ "$(PERLFILES)" != " " ]; then (\
	   gmake -f \$(GNUMAKEFILE) perl \
	) fi

	@if [ "$(ETCFILES)" != "" ]; then (\
	   gmake -f \$(GNUMAKEFILE) etc \
	) fi

	@if [ "$(TXTFILES)" != "" ]; then (\
	   gmake -f \$(GNUMAKEFILE) etc \
	) fi


local_esma_clean local_clean:
	$(RM) -r $(BIN) *.[axol] *.[Mm][Oo][Dd] *~ Makefile*~ rii_files/ \
	         *.aux *.dvi *.ps  *.toc *.log
local_esma_distclean local_distclean:
	$(RM) -r $(BIN) *.[axodl] *.[Mm][Oo][Dd] *~ Makefile*~ rii_files/ \
	         *.aux *.dvi *.ps  *.toc *.log
local_esma_doc local_doc:
	@echo "Target $@ not implemented yet in `pwd`"

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
	@ echo "  make depend"  

########################################################################

bin:
	@echo ""
	@for file in $($@); do  \
	   ( echo "copy to \$$(ESMABIN): $$file"; \
	     $(CP) $$file $(ESMABIN))  \
	done
	@echo ""
	@echo "\$$(ESMABIN) = $(ESMABIN)"
	@echo ""

etc:
	@echo ""
	@for file in $(ALLETC); do \
	   ( case $$file in \
              *.sample)         F=`basename $$file .sample` ;;\
              *.txt)            F=`basename $$file .txt`.rc ;;\
              *)                F=$$f                    ;;\
              esac ;\
              echo "$(CP) -p $$file     $(ESMAETC)/$$F" ;\
              $(CP) -p $$file $(ESMAETC)/$$F )\
	done
	@echo ""
	@echo "\$$(ESMAETC) = $(ESMAETC)"
	@echo ""

perlpkg:
	@echo ""
	@for pkg in $($@); do  \
	   ( echo "perl package to \$$(ESMABIN): $$pkg"; \
	     $(CP) $$pkg $(ESMABIN)/$$pkg )  \
	done
	@echo ""
	@echo "\$$(ESMABIN) = $(ESMABIN)"
	@echo ""

perlscript:
	@echo ""
	@for scr in $($@); do  \
	   ( $(SED) -e "$(SEDCMD)"< $$scr > $(ESMABIN)/$$scr; \
	     echo "perl script to \$$(ESMABIN): $$scr";  \
	     chmod 755 $(ESMABIN)/$$scr ) \
	done
	@echo ""
	@echo "\$$(ESMABIN) = $(ESMABIN)"
	@echo ""

perl:
	@if [ "$(PERLSCRIPT)" != "" ]; then (\
	   gmake -f \$(GNUMAKEFILE) perlscript \
	) fi

	@if [ "$(PERLPKG)" != "" ]; then (\
	   gmake -f \$(GNUMAKEFILE) perlpkg \
	) fi

#                  --------------------
#                  User Defined Targets
#                  --------------------

FOPT = $(FOPT3)
USER_FFLAGS = $(BIG_ENDIAN) $(M). $(M)$(INC_MPEU) $(I)$(INC_SDF)


  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
