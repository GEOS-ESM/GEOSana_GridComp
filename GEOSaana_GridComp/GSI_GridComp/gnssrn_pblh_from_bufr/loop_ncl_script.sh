#!/bin/bash

module load ncl

export YYYY='2021'
export YY='21'
export MM='12'
export MM='12'

#---------------------------------------
# Outer loop for days in December 2021
#---------------------------------------
for DD in {1..31}; do

    if [ "$DD" -lt 10 ]; then
       export DD="0$DD"
    else
       export DD="$DD"
    fi

    #---------------------------------------
    # Inner loop for running every 6 hours
    #---------------------------------------
    for HH in "00" "06" "12" "18"; do
       export HH=$HH
       echo "HH=$HH"

       ncl ./read_nclfile_twopeaks.ncl >& Log_$YYYY$MM$DD$HH

    done
done
