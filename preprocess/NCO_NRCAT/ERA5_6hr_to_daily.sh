#!/bin/bash
######################################################################
# Filename:    ERA5_6hr_to_daily.sh
# Author:      Deanna Nash dlnash@ucsb.edu
# Description: Script to calculate the daily mean of 6-hourly ERA5 data
#
######################################################################
### Activate bash

source /home/sbarc/students/nash/miniconda3/etc/profile.d/conda.sh 

### Activate conda env
conda activate nco-env

# read config file using python
yaml() {
    python -c "import yaml;print(yaml.safe_load(open('$1'))$2$3)"
}

array=(
ivt # 6 hourly IVT climatology over HMA
prec # 6 hourly precip climatology over HMA
huvq # 6 hourly 250 hPa geopotential heights and wind
)
# Loop through each of the variable sets
outer=1      # set outer loop counter
for i in ${!array[*]}
do
    echo "Pass $outer in outer loop."
    input2="['${array[$i]}']"
    data=$(yaml nco_config.yml ${input2} "['data']")
    start_yr=$(yaml nco_config.yml ${input2} "['start_yr']")
    end_yr=$(yaml nco_config.yml ${input2} "['end_yr']")
    domain=$(yaml nco_config.yml ${input2} "['domain']")
    resol=$(yaml nco_config.yml ${input2} "['resol']")
    datadir=$(yaml nco_config.yml ${input2} "['datadir']")
    outdir=$(yaml nco_config.yml ${input2} "['outdir']")

    # Loop to concatenate all files within month-year into 1 netcdf file
    # Begin inner loop (e.g. each year)
    inner=1    # reset inner loop counter
    cd $datadir
    for year in $(seq $start_yr $end_yr)
    do
        echo "Pass $inner in inner loop."
        
        infile="${data}_${region}_${resol}_6hr_${vars}_${year}.nc"
        echo "$infile"
        outfile="${outdir}out.${data}_${domain}_${resol}_daily_${vars}_${year}.nc"
        echo "$outfile"
        # calculate daily mean from 6hourly
        cdo daymean ${infile} ${outfile}

        let "inner+=1" # Increment inner loop counter
        echo "$year concatenation complete"
        echo           # Space between output blocks in pass of inner loop

    done
    
    # End of outer loop
    let "outer+=1" # Increment outer loop counter
    echo "${array[$i]} variable processing complete."
    echo         # Space between output blocks in pass of inner loop
    
done