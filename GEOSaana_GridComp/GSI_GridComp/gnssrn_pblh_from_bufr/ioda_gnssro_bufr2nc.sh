#!/bin/sh

#-----------------------------------------------
yyyy=2012
mm=02
dd=01
hh=06
#-----------------------------------------------

# GNSSRO BUFR directory
export gpsro_bufr=/discover/nobackup/projects/gmao/input/dao_ops/ops/flk/ncep_g5obs/bufr/GPSRO
# JEDI build  (Dan's)
export build=/discover/nobackup/drholdaw/JediDev/develop/build-intel-release
source $build/modules
## excutable to Transfer BUFR to IODA NetCDF
export bufr2nc_x=${build}/bin/bufr2nc_fortran.x 

export locdir=`pwd`
[ ! -d testinput ]  && mkdir testinput
[ ! -d testrun ]  && mkdir testrun

cd ${locdir}/testinput
   yy=`echo $yyyy |cut -b3-4`
   filename=gdas1.${yy}${mm}${dd}.t${hh}z.gpsro.tm00
   ln -s ${gpsro_bufr}/Y$yyyy/M$mm/${filename}.bufr_d  ${filename}.bufr
cd ${locdir}
   $bufr2nc_x -i testinput -o testrun ${filename}.bufr 
