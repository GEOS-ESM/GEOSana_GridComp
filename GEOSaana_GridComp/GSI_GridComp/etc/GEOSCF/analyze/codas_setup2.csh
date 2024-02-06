#!/bin/csh -f

# CODAS_SETUP2	Script for creating RC files from templates
#
# Notes:
# This file is typically run from the codas_run.j job.  However, it is exported
# by that job to a file called codas_setup2.ctmp (along with another file called
# codas_setup1.ctmp) so that in interactive sessions, you can modify the RC files
# and re-run codas_setup1.ctmp and codas_setup2.ctmp to get the appropriate RC
# files.
#
# Author(s):    Brad Weir (brad.weir@nasa.gov)
#
# Changelog:
# 2018/11/28	First crack
#==============================================================================#

/bin/mv GSI_GridComp.rc  GSI_GridComp.tmp
cat     GSI_GridComp.tmp | sed -e "s/^BEG_DATE: .*/BEG_DATE: @BEG_DATE/;         \
                                   s/^END_DATE: .*/END_DATE: @END_DATE/;         \
                                   s/^RECORD_REF_TIME: .*/RECORD_REF_TIME: @RECORD_REF_TIME/; \
                                   s/^RECORD_REF_DATE: .*/RECORD_REF_DATE: @RECORD_REF_DATE/" > GSI_GridComp.rc
