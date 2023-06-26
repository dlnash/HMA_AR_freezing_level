#!/bin/bash
######################################################################
# Filename:    copy_figs.sh
# Author:      Deanna Nash dlnash@ucsb.edu
# Description: Script to copy final figures to one folder
#
######################################################################

# Input parameters
maindir="/home/sbarc/students/nash/repositories/HMA_AR_freezing_level/figs/" # main figure folder
finaldir="/home/sbarc/students/nash/repositories/HMA_AR_freezing_level/figs/final_figs/" # final figure folder

# fig names in main folder
array=(
wrf_geogrid_norris
ivt_AR_trend_DJF_portrait
zero_isotherm_height_AR_trend_DJF_portrait
freezing_sr_diff_composite_portrait
wrf_scatter_rosemax
WRF_prec_summary
synoptic_daily_jan1989
synoptic_daily_feb2010
cross-section_summary_IVT_WV-flux
cross-section_summary_prec_mixing_ratio
wvf_climatology


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
)



for i in ${!array[*]}
do 
    infile="${maindir}${array[$i]}.png"
    outfile="${finaldir}fig${array2[$i]}.png"
#     echo "${infile} to ${outfile}"
    cp -v ${infile} ${outfile}
done

# ### supplemental figs
# supp_array=(
# composite_ar_types_bias
# )

# ## new names to be given
# supp_array2=(
# 1
# )
# for i in ${!supp_array[*]}
# do 
#     infile="${maindir}${supp_array[$i]}.png"
#     outfile="${finaldir}figS${supp_array2[$i]}.png"
# #     echo "${infile} to ${outfile}"
#     cp -v ${infile} ${outfile}
# done

### zip files together
cd ${finaldir}
zip final_figs fig*.png