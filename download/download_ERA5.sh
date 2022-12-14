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
ivt # 6 hourly IVT climatology over HMA
prec # 6 hourly precip climatology over HMA
huvq # 6 hourly 250 hPa geopotential heights and wind
ivt_case_201002 # hourly ivt and precip for 2010 AR case
ivt_case_200201 # hourly ivt and precip for 2002 AR case
250z_case_201002 # hourly 250 hPa for 2010 AR case
250z_case_200201 # hourly 250 hPa for 2002 AR case
)

# now loop through each configuration dictionary to download the ERA5 data

for i in ${!array[*]}
do 
    inconfig="${array[$i]}"
    python getERA5_batch.py ${inconfig}
done