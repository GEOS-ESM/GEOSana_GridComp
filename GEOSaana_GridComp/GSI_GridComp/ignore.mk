#!/usr/bin/env make
#
# NAME:
# 	ignore.mk - a makefile segment to ignore some prerequisites
#
# SYNOPSIS:	
# 	Within a GNUmakefile:
# 	> __ignored__	:= <list of prerequisites to be ignored.
# 	> -include ignore.mk
#
# REVISION HISTORY:
#	2013-05-14	Jing Guo	- initial version
#	2015-05-01	Jing Guo	- simpler, better, and enhanced
#	2015-06-01	Jing Guo	- minor changes in its descriptions.
#

ifdef __ignored__
# This mechanism let make() to believe $(__ignored__) prerequisites are
# present and up to date.
#

__IGNORED_DIR__	?= __ignored__
VPATH += $(__IGNORED_DIR__)
$(__ignored__): |$(__IGNORED_DIR__) ; touch $(__IGNORED_DIR__)/$@
$(__IGNORED_DIR__): ; mkdir -p $@

# Clean after myself, when a distclean target is invoked.
#

local_distckean distclean: __ignored__.distclean
__ignored__.distclean: ; rm -fr $(__IGNORED_DIR__)
endif
