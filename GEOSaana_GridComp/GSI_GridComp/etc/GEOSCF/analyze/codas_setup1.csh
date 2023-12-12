#!/bin/csh -f

# CODAS_SETUP1	Script for creating RC files from templates
#
# Notes:
# This file is typically run from the codas_run.j job.  However, it is exported
# by that job to a file called codas_setup1.ctmp (along with another file called
# codas_setup2.ctmp) so that in interactive sessions, you can modify the RC files
# and re-run codas_setup1.ctmp and codas_setup2.ctmp to get the appropriate RC
# files.
#
# Author(s):    Brad Weir (brad.weir@nasa.gov)
#
# Changelog:
# 2018/10/04	Second-ish crack
# 2018/11/14	Added background type (itracer) column
# 2018/11/28	Added @EXPID field to be filled by codas_run job
#==============================================================================#

# Output debugging information?
set DEBUG  = 0

# Fill lists with info from codas.rc
# ---------------------------------
set tgases = ''
set provs  = ''
set types  = ''
set nbcovs = ''
set ncrtms = ''
set scales = ''
set asnums = ''
set blurbs = ''
set tganas = ''
set tgmons = ''
foreach line ("`cat codas.rc`")
   set firstword  = `echo "$line" | awk -F "," '{print $1}'`
   set firstchar  = `echo "$firstword" | cut -c1`
   set secondword = `echo "$line" | awk -F "," '{print $2}'`

   if ( "$firstword" == "::" ) goto done

   if ( "$firstchar" != "#" & "$firstchar" != "!" ) then
      if ( "$secondword" == "::" ) goto done

      set tgas  = `echo "$line" | awk -F "," '{print $1}' | sed -e "s/'//g"`
      set prov  = `echo "$line" | awk -F "," '{print $2}' | sed -e "s/'//g"`
      set type  = `echo "$line" | awk -F "," '{print $3}' | sed -e "s/'//g"`
      set nbcov = `echo "$line" | awk -F "," '{print $4}'`
      set ncrtm = `echo "$line" | awk -F "," '{print $5}'`
      set scale = `echo "$line" | awk -F "," '{print $6}'`
      set asnum = `echo "$line" | awk -F "," '{print $7}'`
      set blurb = `echo "$line" | awk -F "," '{print $8}' | sed -e "s/'//g"`

      set tgases = `echo $tgases $tgas`
      set provs  = `echo $provs  $prov`
      set types  = `echo $types  $type`
      set nbcovs = `echo $nbcovs $nbcov`
      set ncrtms = `echo $ncrtms $ncrtm`
      set scales = `echo $scales $scale`
      set asnums = `echo $asnums $asnum`
      set blurbs = `echo $blurbs $blurb`

#     Variables for output
#     --------------------
      if ( "$type" == "ANA" ) set tganas = `echo $tganas $tgas`
      if ( "$type" == "MON" ) set tgmons = `echo $tgmons $tgas`
   endif

   if ( "$firstword" == "BACKGROUND:" ) then
      set tgases = ''
      set provs  = ''
      set types  = ''
      set nbcovs = ''
      set ncrtms = ''
      set scales = ''
      set asnums = ''
      set blurbs = ''
      set tganas = ''
      set tgmons = ''
   endif
end

done:
   echo ""
   echo ""
   echo "CoDAS_SETUP: Analyzing  gases ... $tganas"
   echo "             Monitoring gases ... $tgmons"
   echo ""
   echo ""

   if ( "$DEBUG" == 1 ) then
      echo "provs  = $provs"
      echo "types  = $types"
      echo "nbcovs = $nbcovs"
      echo "ncrtms = $ncrtms"
      echo "scales = $scales"
      echo "asnums = $asnums"
      echo "blurbs = $blurbs"
   endif

# Fill in template files
# ----------------------
/bin/cp anavinfo.tmpl         anavinfo
/bin/cp GSI_GridComp.rc.tmpl  GSI_GridComp.rc
/bin/cp CoDAS_ExtData.rc.tmpl CoDAS_ExtData.rc

echo "# Increments to read" > ana2inc.rc

foreach num ( `seq 1 ${#tgases}` )
   set tgas  = ${tgases[$num]}
   set prov  = ${provs[$num]}
   set type  = ${types[$num]}
   set nbcov = ${nbcovs[$num]}
   set ncrtm = ${ncrtms[$num]}
   set scale = ${scales[$num]}
   set asnum = ${asnums[$num]}
   set blurb = ${blurbs[$num]}

#  Set field keywords
#  ------------------
   set gsikey  = ">>>TGASNAME<<<"
   set metkey  = ">>>MET<<<"
   set dtkey   = ">>>DERIV<<<"
   set tendkey = ">>>TEND<<<"
   set chemkey = ">>>CHEM<<<"
   set st8key  = ">>>STATE<<<"
   set conkey  = ">>>CONTROL<<<"
   set inckey  = ">>>INCREMENT<<<"

   set tlow = `echo $tgas | tr "[A-Z]" "[a-z]"`
   set tpad = "`printf  "%-8s" $tlow`"
   set tana = "CODAS_${tgas}_PROVIDER: ${prov}"

#  Hacks for ozone
   set tout = "$tgas"
   set tinc = "$tgas"
   if ( "$tout" == "OZ" )                   set tout = "ozone"
   if ( "$tinc" == "OZ" | "$tinc" == "O3" ) set tinc = "OX"

#  This should only add whitespace, not remove characters
#  ------------------------------------------------------
   set nbcov = "`printf  "%-4s" $nbcov`"
   set ncrtm = "`printf  "%-7s" $ncrtm`"
   set scale = "`printf  "%-8s" $scale`"
   set asnum = "`printf "%-11s" $asnum`"
   set blurb = "`printf "%-18s" $blurb`"

   if ( "$type" == "ANA" | "$type" == "MON" ) then
#     1. Create anavinfo file
#     -----------------------
      if ( "$tout" == "ozone" ) then
         set bund = "met_guess "

         cat     anavinfo | sed -e "/$dtkey/a   \ $tpad 72      met_guess"                     > anavinfo.tmp
         /bin/mv anavinfo.tmp anavinfo
         cat     anavinfo | sed -e "/$tendkey/a \ $tpad 72      met_guess"                     > anavinfo.tmp
         /bin/mv anavinfo.tmp anavinfo
         cat     anavinfo | sed -e "/$metkey/a  \ $tpad 72      $ncrtm    $blurb   $tout"      > anavinfo.tmp
         /bin/mv anavinfo.tmp anavinfo
      else
         set bund = "chem_guess"

         cat     anavinfo | sed -e "/$chemkey/a \ $tpad 72      1         $ncrtm $blurb $tout" > anavinfo.tmp
         /bin/mv anavinfo.tmp anavinfo
      endif
      cat        anavinfo | sed -e "/$st8key/a  \ $tpad 72      1     $bund   $tlow"           > anavinfo.tmp
      /bin/mv    anavinfo.tmp anavinfo
      cat        anavinfo | sed -e "/$conkey/a  \ $tpad 72      $nbcov  $asnum -1.0     state    $tlow     -1.00" > anavinfo.tmp
      /bin/mv    anavinfo.tmp anavinfo

#     2. Create GSI_GridComp.rc file
#     ------------------------------
      cat     GSI_GridComp.rc | sed -e "s/$gsikey/$tout,$gsikey/" > GSI_GridComp.tmp
      /bin/mv GSI_GridComp.tmp  GSI_GridComp.rc
   endif

   if ( "$type" == "ANA" ) then
#     3. Create CoDAS_ExtData.rc file
#     -------------------------------
      set tinc = "`printf "%-19s" "INC_$tinc"`"
      set head = "'mol/mol'           Y N  -                     0.0     "
      set tpup = "`printf "%-13s" $tout`"

      cat     CoDAS_ExtData.rc | sed -e "/$inckey/a  \ $tinc $head $scale $tpup inc.eta.nc4" > CoDAS_ExtData.tmp
      /bin/mv CoDAS_ExtData.tmp CoDAS_ExtData.rc

#     4. Create ana2inc.rc file
#     -------------------------
      echo ${tana} >> ana2inc.rc
   endif
end

# Remove key from GSI_GridComp.rc file
# ------------------------------------
cat     GSI_GridComp.rc | sed -e "s/,$gsikey//"                > GSI_GridComp.tmp
/bin/mv GSI_GridComp.tmp  GSI_GridComp.rc

# Update GSI_GridComp.rc.tmpl to have same EXPID
# -----------------------------------------------
cat     GSI_GridComp.rc | sed -e "s/^expid: .*/expid: @EXPID/" > GSI_GridComp.tmp
/bin/mv GSI_GridComp.tmp  GSI_GridComp.rc

# Create AGCM files for use with CoDAS
# ------------------------------------
cat     AGCM.rc ana2inc.rc > AGCM.rc.corr		# corrector step (increment applied)

# Rename all ozones to OX
# -----------------------
/bin/cp AGCM.rc.corr   AGCM.rc.tmp
cat     AGCM.rc.tmp  | sed -e "s/CODAS_OZ_PROVIDER:\(.*\)/CODAS_OX_PROVIDER:\1/" > AGCM.rc.corr
cat     AGCM.rc.corr | sed -e "s/CODAS_O3_PROVIDER:\(.*\)/CODAS_OX_PROVIDER:\1/" > AGCM.rc.tmp
/bin/mv AGCM.rc.tmp    AGCM.rc.corr

# Build OBSCLASS and OBSDIAGS from gsiparm.anl
# --------------------------------------------
setenv OBSCLASS ','
setenv OBSDIAGS ','
foreach line ("`cat gsiparm.anl`")
   set firstword  = `echo $line | awk '{print $1}'`
   set firstchar  = `echo $firstword | cut -c1`
   set secondword = `echo $line | awk '{print $2}'`

   if ( "$firstword" == "::" ) goto doneobs

   if ( "$firstchar" != "#" & "$firstchar" != "!" & "$OBSCLASS" != "," ) then
      set obsclas = `echo $line | sed -e "s/[[:space:]]\+/ /g" | cut -d' ' -f8`
      set obssis  = `echo $line | sed -e "s/[[:space:]]\+/ /g" | cut -d' ' -f4`
      set obsdiag = `echo diag_${obssis}`

#     Hack to not start with comma
      if ( "$OBSCLASS" == "" ) setenv OBSCLASS $obsclas
      if ( "$OBSDIAGS" == "" ) setenv OBSDIAGS $obsdiag

#     Add if not already in list
      if ( $obsclas !~ {$OBSCLASS} ) setenv OBSCLASS `echo $OBSCLASS,$obsclas`
      if ( $obsdiag !~ {$OBSDIAGS} ) setenv OBSDIAGS `echo $OBSDIAGS,$obsdiag`

      if ( "$secondword" == :: ) goto doneobs
   endif

#  Replace comma when it's time to begin
   if ( "$firstword" == OBS_INPUT:: ) then
      setenv OBSCLASS ''
      setenv OBSDIAGS ''
   endif
end

doneobs:
   setenv OBSDIAGS `echo $OBSDIAGS | sed -e "s/,/ /g"`
   $FVROOT/bin/append_gsigcrc.pl obsys.rc GSI_GridComp.rc
