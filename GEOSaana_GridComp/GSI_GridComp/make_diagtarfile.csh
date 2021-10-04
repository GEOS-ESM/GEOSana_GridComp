#!/bin/csh

# make_diagtarfile - create tar file of radiance diags for new radiance bc
#
# Usage:
#   make_diagtarfile.csh $EXPID nymd nhms tarfilename
# !REVISION HISTORY:
#
#  08Jul2014  Sienkiewicz   original version
#  02Jun2015  Sienkiewicz   using Joe Stassi method from 'gsidiags' to
#                               get instrument names
#  spring2017 Sienkiewicz   grep -v to remove pcp from master file list
#                               rather than using a custom rc file
#------------------------------------------------------------------

set EXPID = $1
set nymd = $2
set nhms = $3

set echorc_x    = $FVROOT/bin/echorc.x
set gsidiags_rc = $FVROOT/etc/gsidiags.rc

set hr = `echo $nhms | cut -c1-2`
set dtg = ${nymd}_${hr}z

set tarfile = $4
/bin/rm -f $tarfile

foreach sis ( `$echorc_x -rc $gsidiags_rc satlist | grep -v pcp` )
   set dfile = $EXPID.diag_${sis}_ges.${dtg}.bin
   if ( -s $dfile ) then
      if ( ! -s $tarfile ) then
        tar -cf $tarfile $dfile
      else
        tar -uf $tarfile $dfile
      endif
   endif
end
