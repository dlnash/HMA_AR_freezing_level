#!/bin/bash
######################################################################
# Filename:    download_ERA5.sh
# Author:      Deanna Nash dlnash@ucsb.edu
# Description: Script to download all the necessary ERA5 files for HMA_AR_freezing_level
#
######################################################################

### Activate bash

source /home/sbarc/students/nash/miniconda3/etc/profile.d/conda.sh 

### Activate conda env

conda activate cds

# names of configuration dictionaries to loop through
array=(
ivt_case_201002 # hourly ivt and precip for 2010 AR case
ivt_case_198901 # hourly ivt and precip for 1989 AR case
250z_case_201002 # hourly 250 hPa for 2010 AR case
250z_case_198901 # hourly 250 hPa for 1989 AR case
)

# now loop through each configuration dictionary to download the ERA5 data

for i in ${!array[*]}
do 
    inconfig="${array[$i]}"
    python getERA5_batch.py ${inconfig}
done