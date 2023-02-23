"""
Filename:    preprocess_wrf_climatology.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: preprocess chosen variables from 36 years of wrf output and save as yearly nc files
"""

## Imports
import os, sys
import glob
import numpy as np
import xarray as xr
import wrf
import netCDF4 as nc
from scipy.integrate import trapz
# Necessary to add cwd to path when script run
# by SLURM (since it executes a copy)
sys.path.append(os.getcwd())

# Path to modules
sys.path.append('../../modules')
# Import my modules
from wrf_funcs_preprocess import preprocess_prec, preprocess_ivt, preprocess_SR, preprocess_pressure_lev_var

## Set up paths
path_to_wrf = '/home/hasia/'
path_to_data = '/home/nash/DATA/data/'                                      # project data -- read only

# command for selecting variables **include XTIME,XLONG,XLAT in var list
# ncks -v varname1,varname2 input_file output_file

### Input ###

# ### GEOPOTENTIAL HEIGHT ###
# # variable to interpolate
# var1_name = 'z'
# # vertical variable to interpolate to (e.g. pressure, height)
# var2_name = 'pressure'
# # level to interpolate to
# levs = [1000, 850, 500, 250]
# output_varname = 'geopotential'
# domain = 'd01'

### ZERO DEGREE ISOTHERM ###
# variable to interpolate
var1_name = 'z'
# vertical variable to interpolate to (e.g. pressure, height)
var2_name = 'tc'
# level to interpolate to
levs = [0]
output_varname = 'zerodegisotherm'
domain = 'd01'

# ### IVT ###
# output_varname = 'ivt'
# domain = 'd01'
# wrf.omp_set_num_threads(8)

# ### PRECIP ###
# output_varname = 'prec'
# domain = 'd02'

# ### SR ###
# output_varname = 'sr'
# domain = 'd02'

### START PROGRAM ###
start_yr = 2015
end_yr = 2015

## Loop through all above years
for year in np.arange(start_yr, end_yr+1):
    print('Processing ...', year)
    if year == 1979:
        # get list of filenames that contain data from that year from current year folder
        filenames = []
        for name in glob.glob('/home/hasia/' + str(year) + '/wrfout_{0}_*'.format(domain)):
            filenames.append(name)
        # sort filenames so they are in chronological order
        filenames = sorted(filenames)
        
        # just get files between last 25 and last 18
        filenames = filenames[-25:-18]
        ## create new ds depending on input above
        if output_varname == 'prec':
            ds = preprocess_prec(filenames)

        elif output_varname == 'ivt':
            ds = preprocess_ivt(filenames)
            
        elif output_varname == 'sr':
            ds = preprocess_SR(filenames)

        else:
            ds = preprocess_pressure_lev_var(filenames, var1_name, var2_name, levs)

        # clip the time steps so it is just Nov 30 - Dec 31
        ds = ds.sel(time=slice('{0}-11-30T00:00'.format(str(year)), '{0}-12-31T21:00'.format(str(year))))
        
    elif year == 2015:
        # get list of filenames that contain data from that year from current year folder
        filenames = []
        for name in glob.glob('/home/hasia/' + str(year-1) + '/wrfout_{0}_*'.format(domain)):
            filenames.append(name)
        # sort filenames so they are in chronological order
        filenames = sorted(filenames)
        
        # just get last 18 files
        filenames = filenames[-18:]
        ## create new ds depending on input above
        if output_varname == 'prec':
            ds = preprocess_prec(filenames)

        elif output_varname == 'ivt':
            ds = preprocess_ivt(filenames)
        elif output_varname == 'sr':
            ds = preprocess_SR(filenames)

        else:
            ds = preprocess_pressure_lev_var(filenames, var1_name, var2_name, levs)

        # clip the time steps so it is just Jan 1 to March 31
        ds = ds.sel(time=slice('{0}-12-31T00:00'.format(str(year-1)), '{0}-03-31T21:00'.format(str(year))))
        
    else:
        # get list of filenames that contain data from that year from previous year folder
        filenames = []
        for name in glob.glob('/home/hasia/' + str(year-1) + '/wrfout_{0}_*'.format(domain)):
            filenames.append(name)
        # sort filenames so they are in chronological order
        filenames = sorted(filenames)
        # only get last 18 files
        filenames = filenames[-18:]

        # now let's get the first x files in the current year folder
        filenames2 = []
        for name in glob.glob('/home/hasia/' + str(year) + '/wrfout_{0}_*'.format(domain)):
            filenames2.append(name)
        # sort filenames so they are in chronological order
        filenames2 = sorted(filenames2)
        filenames2 = filenames2[:-17]

        print('Number of files:', len(filenames))
        print(filenames)

        print('Number of files:', len(filenames2))
        print(filenames2)

        ## create new ds depending on input above
        if output_varname == 'prec':
            ds1 = preprocess_prec(filenames)
            ds2 = preprocess_prec(filenames2)

        elif output_varname == 'ivt':
            ds1 = preprocess_ivt(filenames)
            ds2 = preprocess_ivt(filenames2)
        elif output_varname == 'sr':
            ds1 = preprocess_SR(filenames)
            ds2 = preprocess_SR(filenames2)

        else:
            ds1 = preprocess_pressure_lev_var(filenames, var1_name, var2_name, levs)
            ds2 = preprocess_pressure_lev_var(filenames2, var1_name, var2_name, levs)

        # clip the first 8 time steps and last 33 time steps off of ds2
        ds1 = ds1.sel(time=slice('{0}-01-01T00:00'.format(str(year)), '{0}-03-31T21:00'.format(str(year))))
        ds2 = ds2.sel(time=slice('{0}-04-01T00:00'.format(str(year)), '{0}-12-31T21:00'.format(str(year))))
        
        # make sure everyone has the same lats/lons
        wrflats = ds1.lat.values
        wrflons = ds1.lon.values
        # reassign lats/lons
        ds1 = ds1.assign_coords({"lon": wrflons, "lat": wrflats})
        ds2 = ds2.assign_coords({"lon": wrflons, "lat": wrflats})
        
        # merge ds1 and ds2
        ds = xr.merge([ds1, ds2])
    
    # write to netCDF
    fname = os.path.join(path_to_data, 'wrf_hasia/{0}/{1}/3hr/tmp_{2}.nc').format(domain, output_varname, str(year))
    ds.to_netcdf(path=fname, mode = 'w', format='NETCDF4')