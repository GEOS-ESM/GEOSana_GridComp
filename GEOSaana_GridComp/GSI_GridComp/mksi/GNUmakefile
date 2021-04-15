#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 07Feb2008  Jing Guo   Based on GNUmakefile of GSI_GridComp
# 14Apr2012  Todling    Add gsiiinfo.pl to control info creation
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../../../../../..
endif
.PHONY: esma_install local_install      # Mark targets without action to avoid default %_install rule

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------


THIS = $(shell basename `pwd`)
BIN  =  gsiinfo.pl \
        make_convinfo.x	\
	make_ozinfo.x	\
	make_satinfo.x

DATA =	convinfo.db	\
	ozinfo.db	\
	sidb

DIRS = $(ESMAETC) $(ESMABIN)
DB_DIRS = $(ESMAETC)/gmao_satinfo.db	\
	  $(ESMAETC)/gmao_ozinfo.db	\
	  $(ESMAETC)/gmao_convinfo.db

	# My default make-target is not "install" for local convenience.
my-default: test

esma_install install: install_bin install_etc
local_install: install_bin install_etc

$(DB_DIRS):	$(ESMAETC)
$(ESMABIN) $(ESMAETC) $(DB_DIRS):
	$(MKDIR) $@

install_bin: $(ESMABIN) $(BIN)
	  $(CP) -p $(BIN) $(ESMABIN)/

install_etc: install_satinfo install_ozinfo install_convinfo
install_satinfo: $(ESMAETC)/gmao_satinfo.db sidb
	  $(CP) -prv sidb/* $(ESMAETC)/gmao_satinfo.db/
install_ozinfo: $(ESMAETC)/gmao_ozinfo.db ozinfo.db
	  $(CP) -prv ozinfo.db/* $(ESMAETC)/gmao_ozinfo.db/
install_convinfo: $(ESMAETC)/gmao_convinfo.db convinfo.db
	  $(CP) -prv convinfo.db/* $(ESMAETC)/gmao_convinfo.db/

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

#                  --------------------
#                  User Defined Targets
#                  --------------------


SRCS_SATINFO :=	main.F90	\
		m_sitmpl.F90
OBJS_SATINFO := $(addsuffix .o, $(basename $(SRCS_SATINFO)))

SRCS_OZINFO  :=	make_ozinfo.F90	\
		m_oztmpl.F90
OBJS_OZINFO  := $(addsuffix .o, $(basename $(SRCS_OZINFO)))

SRCS_CONVINFO := make_convinfo.F90 \
		 m_convtmpl.F90 m_actvotype.F90
OBJS_CONVINFO  := $(addsuffix .o, $(basename $(SRCS_CONVINFO)))

SRCS_OTHERS  :=	m_actvchan.F90	\
		satinfo_util.F90
OBJS_OTHERS  :=	$(addsuffix .o, $(basename $(SRCS_OTHERS)))

SRCS := $(SRCS_SATINFO) $(SRCS_OZINFO) $(SRCS_CONVINFO) $(SRCS_OTHERS)
OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS)))

make_satinfo.x: $(OBJS_SATINFO) $(OBJS_OTHERS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS_SATINFO) $(OBJS_OTHERS)

make_ozinfo.x: $(OBJS_OZINFO) $(OBJS_OTHERS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS_OZINFO) $(OBJS_OTHERS)

make_convinfo.x: $(OBJS_CONVINFO) $(OBJS_OTHERS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS_CONVINFO) $(OBJS_OTHERS)

list:
	@ echo $(SRCS)

test: test_make_satinfo test_make_ozinfo test_make_convinfo

test_make_satinfo: make_satinfo.x setup.nml
	: ## -- $@
	./make_satinfo.x <setup.nml

ozinfo.tmpl: ../etc/gmao_global_ozinfo.txt
	cp $< $@
test_make_ozinfo ozinfo.txt: ozinfo.tmpl make_ozinfo.x
	: ## -- $@: $?
	echo \&setup nymd=20021231 nhms=000000 info_tmpl=$< info_outf=ozinfo.txt notused_flag=-2 / | make_ozinfo.x
	: ## -- diff $< ozinfo.txt
	diff  $< ozinfo.txt; exit 0

convinfo.tmpl: ../etc/gmao_global_convinfo.txt
	cp $< $@
test_make_convinfo convinfo.txt: convinfo.tmpl make_convinfo.x
	: ## -- $@
	echo \&setup nymd=20021231 nhms=000000 info_tmpl=$< info_outf=convinfo.txt notused_flag=-2 / | make_convinfo.x
	: ## -- diff $< convinfo.txt
	diff  $< convinfo.txt; exit 0

conv_pfix := convinfo.db/
new-convinfo-tbls: convinfo.tmpl
	./convtbl_extract convinfo.tmpl $(conv_pfix)available.tbl $(conv_pfix)active.tbl

_D =

FOPT       = $(FOPT3)
FREAL      = $(FREAL4) 
FPE        =

USER_FDEFS	= $(_D)
USER_FFLAGS	= $(BIG_ENDIAN) $(BYTERECLEN)
USER_FFLAGS	=
USER_CFLAGS	=

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
#.
