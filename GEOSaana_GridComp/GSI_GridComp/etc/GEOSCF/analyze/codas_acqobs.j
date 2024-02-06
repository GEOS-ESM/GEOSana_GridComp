#!/bin/csh -xvf

#######################################################################
#                     Batch Parameters for Archive Job
#######################################################################

#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --job-name=cdas_test_ACQO
#SBATCH --partition=datamove
#SBATCH --account=s1866
#SBATCH --output=@ACQOBS_O

#######################################################################
#                  System Environment Variables
#######################################################################

umask 022

limit stacksize unlimited

   setenv MV2_ON_DEMAND_THRESHOLD 8192
   setenv MV2_USE_SHMEM_ALLREDUCE 0
   setenv MV2_USE_SHMEM_COLL      0
   setenv MV2_USE_UD_HYBRID       0

#######################################################################
#           Architecture Specific Environment Variables
#######################################################################

setenv FVROOT @FVROOT
source $FVROOT/bin/g5_modules

#######################################################################
#                         Archive Commands
#######################################################################

set path = ( $FVROOT/bin $path )

