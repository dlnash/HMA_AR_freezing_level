#!/bin/bash
######################################################################
# Filename:    copy_figs.sh
# Author:      Deanna Nash dlnash@ucsb.edu
# Description: Script to copy final figures to one folder
#
######################################################################

# Input parameters
maindir="/home/sbarc/students/nash/repositories/HMA_AR_freezing_level/figs/" # main figure folder
finaldir="/home/sbarc/students/nash/repositories/HMA_AR_freezing_level/figs/final_figs_ch3/" # final figure folder

# fig names in main folder
array=(
wrf_geogrid_norris
composite_ar_types_bias
arfreq_AR_trend_mon
ivt_AR_trend_DJF_portrait
zero_isotherm_height_AR_trend_DJF_portrait
freezing_sr_diff_composite_portrait
ivt_precip_scatter_era5_max_3col
jan2002_era5_IVT_6hrly
feb2010_era5_IVT_6hrly
jan2002_WRF_prec_3hr
feb2010_WRF_prec_3hr
WRF_prec_summary
jan2002_summary_daily
feb2010_summary_daily
IVT_hovmoller
prec_hovmoller
time_height
)

# new names to be fig<name given in array2>
array2=(
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
)

for i in ${!array[*]}
do 
    infile="${maindir}${array[$i]}.png"
    outfile="${finaldir}fig${array2[$i]}.png"
#     echo "${infile} to ${outfile}"
    cp -v ${infile} ${outfile}
done

