#!/bin/bash
######################################################################
# Filename:    copy_data_pack.sh
# Author:      Deanna Nash dlnash@ucsb.edu
# Description: Script to copy data from multiple directories to single data directory
#
######################################################################

# Input parameters
maindir="" # main figure folder
finaldir="" # final figure folder

# data names in original
array=(

################
### DOWNLOAD ###
################

## # OTHER downloaded data
# "/home/sbarc/students/nash/data/ar_catalog/globalARcatalog_ERA-Interim_1979-2019_v3.0.nc"
# "/home/sbarc/students/nash/data/CH2_generated_data/Global_Landslide_Catalog_Export.csv"

## ERA5 downloaded data
# "/home/sbarc/students/nash/data/ERA5/era5_hma_025deg_1hr_ivt_198901_1989.nc" 
# "/home/sbarc/students/nash/data/ERA5/era5_hma_025deg_1hr_250z_198901_1989.nc"
# "/home/sbarc/students/nash/data/ERA5/era5_hma_025deg_1hr_ivt_201002_2010.nc"
# "/home/sbarc/students/nash/data/ERA5/era5_hma_025deg_1hr_250z_201002_2010.nc"


#############
### DRYAD ###
#############

# ## WRF NORRIS (60 GB)
# "/home/hasia/2009/wrfout_d02_2010-02-04_03:00:00" 
# "/home/hasia/2009/wrfout_d01_2010-02-04_03:00:00" 
# "/home/hasia/1988/wrfout_d02_1988-12-31_03:00:00" 
# "/home/hasia/1988/wrfout_d02_1989-01-05_03:00:00"
# "/home/hasia/1988/wrfout_d01_1988-12-31_03:00:00" 
# "/home/hasia/1988/wrfout_d01_1989-01-05_03:00:00"

# ## WVF (846 KB)
# "/home/sbarc/students/nash/data/wrf_hasia/d02/wvflux_34.87N_72.66E.nc"
# "/home/sbarc/students/nash/data/wrf_hasia/d02/wvflux_34.09N_74.02E.nc"

## Zero Degree Isotherm (19 GB)
# "/home/sbarc/students/nash/data/wrf_hasia/d01/zerodegisotherm/3hr/tmp_*.nc"

# ## IVT (150 GB)
# "/home/sbarc/students/nash/data/wrf_hasia/d01/ivt/3hr/out.wrf.d01.ivt.3hr_*.nc"



####################
### PREPROCESSED ###
####################

## WRF preprocessed data
# "/home/sbarc/students/nash/data/wrf_hasia/d01/ivt/daily/out.wrf6km.ivt.daily_*.nc" # 19 GB
# "/home/sbarc/students/nash/data/wrf_hasia/d01/zerodegisotherm/daily/out.wrf6km.zerodegisotherm.daily_*.nc" # 3.6 GB
# "/home/sbarc/students/nash/data/wrf_hasia/d01/geopotential/daily/out.wrf6km.geopotential.daily_*.nc" # 9 GB
# "/home/sbarc/students/nash/data/wrf_hasia/d02/prec/daily/out.wrf6km.prec.daily_*.nc" # 6.5 GB

# "/home/sbarc/students/nash/data/wrf_hasia/d01/ivt/filtered_daily_mean_clim_ivt.nc"
# "/home/sbarc/students/nash/data/wrf_hasia/d01/ivt/daily_std_clim_ivt.nc"
# "/home/sbarc/students/nash/data/wrf_hasia/d01/ivt/daily_mean_clim_ivt.nc"
# "/home/sbarc/students/nash/data/wrf_hasia/d01/zerodegisotherm/filtered_daily_mean_clim_zerodegisotherm.nc"
# "/home/sbarc/students/nash/data/wrf_hasia/d01/zerodegisotherm/daily_std_clim_zerodegisotherm.nc"
# "/home/sbarc/students/nash/data/wrf_hasia/d01/zerodegisotherm/daily_mean_clim_zerodegisotherm.nc"
# "/home/sbarc/students/nash/data/wrf_hasia/d01/zerodegisotherm_ar_std_new.nc"
# "/home/sbarc/students/nash/data/wrf_hasia/d01/zerodegisotherm_ar_clim_new.nc"
# "/home/sbarc/students/nash/data/wrf_hasia/d01/geopotential/filtered_daily_mean_clim_geopotential.nc"
# "/home/sbarc/students/nash/data/wrf_hasia/d01/geopotential/daily_std_clim_geopotential.nc"
# "/home/sbarc/students/nash/data/wrf_hasia/d01/geopotential/daily_mean_clim_geopotential.nc"

)

# new names in given location
array2=(

################
### DOWNLOAD ###
################

# # OTHER
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/downloads/globalARcatalog_ERA-Interim_1979-2019_v3.0.nc"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/downloads/Global_Landslide_Catalog_Export.csv"

## ERA5 downloaded data
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/downloads/ERA5/era5_hma_025deg_1hr_ivt_198901.nc" 
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/downloads/ERA5/era5_hma_025deg_1hr_250z_198901.nc"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/downloads/ERA5/era5_hma_025deg_1hr_ivt_201002.nc"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/downloads/ERA5/era5_hma_025deg_1hr_250z_201002.nc"

#############
### DRYAD ###
#############

# # WRF NORRIS (60 GB)
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/dryad/WRF_norris/wrfout_d02_2010-02-04_03:00:00"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/dryad/WRF_norris/wrfout_d01_2010-02-04_03:00:00" 
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/dryad/WRF_norris/wrfout_d02_1988-12-31_03:00:00" 
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/dryad/WRF_norris/wrfout_d02_1989-01-05_03:00:00"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/dryad/WRF_norris/wrfout_d01_1988-12-31_03:00:00" 
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/dryad/WRF_norris/wrfout_d01_1989-01-05_03:00:00"

# ## WVF (846 KB)
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/dryad/wvf/wvflux_34.87N_72.66E.nc"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/dryad/wvf/wvflux_34.09N_74.02E.nc"

## Zero Degree Isotherm (19 GB)
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/dryad/zerodegisotherm/out.wrf.d01.zerodegisotherm.3hr_*.nc"

# ## IVT (150 GB)
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/dryad/ivt/out.wrf.d01.ivt.3hr_*.nc"

####################
### PREPROCESSED ###
####################

# "/home/sbarc/students/nash/data/HMA_freezing_level_data/preprocessed/ivt/out.wrf.d01.ivt.daily_*.nc" # 19 GB
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/preprocessed/zerodegisotherm/out.wrf.d01.zerodegisotherm.daily_*.nc" # 3.6 GB

# "/home/sbarc/students/nash/data/HMA_freezing_level_data/preprocessed/ivt/filtered_daily_mean_clim_ivt.nc"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/preprocessed/ivt/daily_std_clim_ivt.nc"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/preprocessed/ivt/daily_mean_clim_ivt.nc"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/preprocessed/zerodegisotherm/filtered_daily_mean_clim_zerodegisotherm.nc"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/preprocessed/zerodegisotherm/daily_std_clim_zerodegisotherm.nc"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/preprocessed/zerodegisotherm/daily_mean_clim_zerodegisotherm.nc"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/preprocessed/zerodegisotherm/ar_std_clim_zerodegisotherm.nc"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/preprocessed/zerodegisotherm/ar_mean_clim_zerodegisotherm.nc"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/preprocessed/geopotential/filtered_daily_mean_clim_geopotential.nc"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/preprocessed/geopotential/daily_std_clim_geopotential.nc"
# "/home/sbarc/students/nash/data/HMA_freezing_level_data/preprocessed/geopotential/daily_mean_clim_geopotential.nc"

)

# for i in ${!array[*]}
# do 
#     infile=${array[$i]}
#     outfile=${array2[$i]}
#     echo "${infile} to ${outfile}"
#     mv -- ${infile} ${outfile}
# done

################################################################################################################
### This section of code is for moving iteratively (files that have same name structure and year as wildcard ###
################################################################################################################



start_yr=1979
end_yr=2015
# in_path="/home/sbarc/students/nash/data/wrf_hasia/d01/geopotential/daily/"
in_path="/home/sbarc/students/nash/data/HMA_freezing_level_data/preprocessed/ivt/"
out_path="/home/sbarc/students/nash/data/HMA_freezing_level_data/preprocessed/ivt/"
for year in $(seq $start_yr $end_yr) 
do
    infile="${in_path}out.wrf.d01.ivt.3hr_${year}.nc"
    outfile="${out_path}out.wrf.d01.ivt.daily_${year}.nc"
    
    mv -- ${infile} ${outfile}
done