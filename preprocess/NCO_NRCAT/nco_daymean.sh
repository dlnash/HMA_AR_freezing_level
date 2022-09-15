#!/bin/bash
######################################################################
# Filename:    nco_wrf_3hr_to_daily_mean.sh
# Author:      Deanna Nash dlnash@ucsb.edu
# Description: Script to get daily mean from 3-hourly WRF
#
# - Convert 3-hourly WRF data to daily
# - To run: change input parameters, then go to NCO_NRCAT directory and activate bash by typing "bash"
# - "conda activate nco-env", then run "bash nco_wrf_3hr_to_daily_mean.sh"
######################################################################

# Input parameters
varname="sr"
domain="d02"
start_yr=1979
end_yr=2015
datadir="/home/sbarc/students/nash/data/wrf_hasia/${domain}/${varname}/3hr/"
outdir="/home/sbarc/students/nash/data/wrf_hasia/${domain}/${varname}/daily/"

outer=1      # set outer loop counter

# Loop to concatenate all files within month-year into 1 netcdf file
# Begin outer loop (e.g. each year)
cd $datadir
for year in $(seq $start_yr $end_yr)
do
    echo "Pass $outer in outer loop."
    inner=1    # reset inner loop counter
    infile="${datadir}tmp_${year}.nc"
    echo "$infile"
    outfile="${outdir}out.wrf6km.${varname}.daily_${year}.nc"
    echo "$outfile"
    # calculate daily mean from 6hourly
    cdo daymean ${infile} ${outfile}

    let "outer+=1" # Increment outer loop counter
    echo "$year concatenation complete"
    echo           # Space between output blocks in pass of outer loop

done
# End of outer loop

exit 0