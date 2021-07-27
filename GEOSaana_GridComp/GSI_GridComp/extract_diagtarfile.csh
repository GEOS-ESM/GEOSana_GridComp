#!/bin/csh

# extract_diagtarfile - create tar file of radiance diags for new radiance bc
#
# !REVISION HISTORY:
#
#  08Jul2014  Sienkiewicz   original version
#------------------------------------------------------------------

set tarfile = $1

if ( -s $tarfile ) then

  set filelist = `tar tf $tarfile`

  foreach file ( $filelist )
     set fname = `echo $file | cut -d. -f2 | sed 's/_ges//' `
     tar xf $tarfile $file
     mv $file $fname
  end
endif
