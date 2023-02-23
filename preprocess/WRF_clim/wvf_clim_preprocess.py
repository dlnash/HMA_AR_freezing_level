"""
Filename:    wvf_clim_preprocess.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: preprocess WVF from 36 years of wrf output at selected coordinates and save as nc file
"""

import os, sys
import glob
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc
import wrf
import datetime
import xoak
import cartopy.crs as ccrs
# import personal modules
# Path to modules
sys.path.append('../../modules')

# Import my modules
from wrf_funcs_preprocess import select_single_coord_WRF
from ar_funcs import get_ar_days, duration_stats
from timeseries import select_months_ds, select_months_df
from utils import find_perpindicular_line, find_parallel_line, find_intersection_two_lines

# Set up paths

path_to_data = '/home/sbarc/students/nash/data/'
path_to_out  = '../../out/'       # output files (numerical results, intermediate datafiles) -- read & write
path_to_figs = '../../figs/'      # figures

## Get coordinates of locations needed for WVF climatology
## Get lines and points for mesoscale analysis

# hlat, hlon, tlat, tlon
Line1 = [39.5, 71., 32.25, 81.]
x1, y1, x2, y2 = Line1[3],Line1[2],Line1[1],Line1[0]

# calculate parallel lines
Line2, eq1, eq2 = find_parallel_line(x1, y1, x2, y2, -3., 69.,  80.)
Line3, eq1, eq3 = find_parallel_line(x1, y1, x2, y2, -2., 69.,  80.)

## calculate perpindicular cross section lines
newx_lst = [74.5, 76., 76.25, 78.5]
newx_lst = [74.25, 76.25]
newline = []
ptlst = []
for i, newx in enumerate(newx_lst):
    # get perpindicular line
    line, eq = find_perpindicular_line(x1, y1, x2, y2, newx)
    newline.append(line)
    # get intersecting point for each newline and Line3
    pt = find_intersection_two_lines(eq2[0], eq2[1], eq[0], eq[1])
    ptlst.append(pt)
# print(ptlst)

## add points from radiosonde locations
# ptlst.append([71.38, 42.85]) # Station 38341
# ptlst.append([77.20, 28.58]) # station 42182
ptlst.append([74.0225, 34.0876]) # Feb2010 Landslide Loc
ptlst.append([72.6625, 34.8733]) # Feb2010 Landslide Loc

ptlst = ptlst[2:] ## already processed first two points

## Get list of AR dates for specified lat/lon

def ar_dates_single_coord(start_date, end_date, lat1, lon1, mons, mone):
    filename =  'globalARcatalog_ERA-Interim_1979-2019_v3.0.nc'

    # open ds
    ds = xr.open_dataset(path_to_data + 'ar_catalog/' + filename, chunks={'time': 1460}, engine='netcdf4')
    ds = ds.squeeze()
    # remove lev and ens coords
    ds = ds.reset_coords(names=['lev', 'ens'], drop=True)

    
    # select dates within start_date, end_date and months
    ds = ds.sel(time=slice(start_date, end_date))
    ## select only months we are interested in
    ds = select_months_ds(ds, mons, mone)
    
    ## select single coordinate
    ds = ds.sel(lat=lat1, lon=lon1, method="nearest")


    
    ## now get list of dates where ds.shape != np.nan
    # convert dataset to dataframe
    df = ds.shape.to_dataframe(dim_order=['time'])
    df = df.dropna(axis='rows')
    
    dates = df.index.values
    
    return df
    
def wrf_clim_filenames(year, domain):
        
    if year == 2015:
        # get list of filenames that contain data from that year from current year folder
        filenames = []
        for name in glob.glob('/home/hasia/' + str(year-1) + '/wrfout_{0}_*'.format(domain)):
            filenames.append(name)
        # sort filenames so they are in chronological order
        filenames = sorted(filenames)
        
        # just get last 18 files
        filenames = filenames[-18:-5]
        
    else:
        # get list of filenames that contain data from that year from previous year folder
        filenames = []
        for name in glob.glob('/home/hasia/' + str(year-1) + '/wrfout_{0}_*'.format(domain)):
            filenames.append(name)
        # sort filenames so they are in chronological order
        filenames = sorted(filenames)
        # only get last 18 files
        filenames = filenames[-18:-5]
        
    return filenames

######### BEGIN INPUT ############
start_date = '1980-01-01'
end_date = '2015-02-28'
mons = 1
mone = 2
domain ='d02'

dates = []
for i, pt in enumerate(ptlst):
    times = ar_dates_single_coord(start_date, end_date, pt[1], pt[0], mons, mone)
    dates.append(times)
    
## Now get the WRF wvf values on those datetimes
varlst = ('pres', 'ua', 'va', 'QVAPOR', 'HGT')
## for each list of dates, 
yr_lst = np.arange(1980, 2016, 1)
for j, pt in enumerate(ptlst):
    lat = pt[1]
    lon = pt[0]
    ds_lst = []
    for i, yr in enumerate(yr_lst):
        times = dates[j]
        # get the list of dates and times within that year
        idx = (times.index.year == yr)
        timesteps = times.loc[idx].index.values
        # now get the wrf filenames for that year between jan/feb
        wrf_filenames = wrf_clim_filenames(yr, domain)

        if len(timesteps) == 0:
            print('No ARs during', yr)
        else:
            # now we can use this function to search the given wrf files for those timesteps and process the vertical values for the given location
            ds = select_single_coord_WRF(wrf_filenames, varlst, lat, lon, timesteps)

            # # calculate specific humidity from mixing ratio (units: kg kg-1)
            # ds = ds.assign(q=lambda ds: np.divide(ds.QVAPOR, (1+ds.QVAPOR)))

            # compute vertical moisture flux (units: m s-1*kg kg-1)
            uq = ds.ua*ds.q
            vq = ds.va*ds.q
            ds = ds.assign(wvf=lambda ds: np.sqrt(uq**2 + vq**2))
            ds_lst.append(ds)
    
    # concat/merge ds_lst to single ds file
    ds_final = xr.concat(ds_lst, dim='Time')
    
    # write to netCDF
    varname = 'wvflux'
    latstr = '{:.2f}N'.format(lat)
    lonstr = '{:.2f}E'.format(lon)
    fname = path_to_data + 'wrf_hasia/{0}/{1}_{2}_{3}.nc'.format(domain, varname, latstr, lonstr)
    ds_final.to_netcdf(path=fname, mode = 'w', format='NETCDF4')