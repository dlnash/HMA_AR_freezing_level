#!/usr/bin/python3
"""
Filename:    wrf_funcs_preprocess.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: preprocess functions for wrf 3 km output
"""

## Imports

import os, sys
import yaml
import glob
from scipy.integrate import trapz
import wrf
import netCDF4 as nc
import numpy as np
import xarray as xr
import pandas as pd

def get_wrf_filenames(casename, domain, server, arname):
    if server == 'frontera':
        # Set up paths
        path_to_data = '/home1/08540/dlnash/DATA_WRF/{0}/AnalysisData/'.format(casename)

        ## WRF filenames
        if domain == 'd01':
            ## Domain 1
            filenames = [path_to_data + "wrfout_d01.2010-02-03_00:00:00",
                             path_to_data + "wrfout_d01.2010-02-06_00:00:00",
                             path_to_data + "wrfout_d01.2010-02-09_00:00:00"]
        elif domain == 'd02':
            # Domain 2
            filenames = [path_to_data + "wrfout_d02.2010-02-03_00:00:00",
                             path_to_data + "wrfout_d02.2010-02-04_00:00:00",
                             path_to_data + "wrfout_d02.2010-02-05_00:00:00",
                             path_to_data + "wrfout_d02.2010-02-06_00:00:00",
                             path_to_data + "wrfout_d02.2010-02-07_00:00:00",
                             path_to_data + "wrfout_d02.2010-02-08_00:00:00",
                             path_to_data + "wrfout_d02.2010-02-09_00:00:00"]
    if server == 'great':
        # import configuration file for case study choice
        yaml_doc = '../data/ar_casestudy.yml'
        config = yaml.load(open(yaml_doc), Loader=yaml.SafeLoader)
        ar_dict = config[arname]
        
        if domain == 'd01':
            filenames = ar_dict['wrf_files']
        elif domain == 'd02':
            filenames = ar_dict['wrf_files2']
    
    return filenames

def preprocess_ivt(filenames):
    """preprocess_ivt
    
    Returns a ds object (xarray) including ivtu, ivtv, and iwv
    
    Parameters
    ----------
    filenames : list
        list of wrf filenames to process
  
    Returns
    -------
    ds : ds object
        includes variables ivtu, ivtv, and iwv
    
    """
    # levels to interpolate to
    interp_levs = [1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700, 650, 600, 550, 500, 450, 400, 250, 300]
    
    # arrays to append data
    ivtu_final = []
    ivtv_final = []
    iwv_final = []
    da_time = []

    for i, wrfin in enumerate(filenames):
        f = nc.Dataset(wrfin, 'r')
        # get the number of times in the current file
        ntimes = f.variables['Times'].shape[0]

        # now loop through the total number of times for this file and append data:
        for file_timeidx in range(ntimes):
            print(wrfin, 'timeidx is ', file_timeidx)
            
            # interpolate ua, va, and qvapor
            ua = wrf.getvar(f, 'ua', file_timeidx)
            interp_ua = wrf.vinterp(f,
                                   field=ua,
                                   vert_coord="pressure",
                                   interp_levels=interp_levs,
                                   extrapolate=False,
                                   field_type=None,
                                   log_p=True,
                                   timeidx=file_timeidx,
                                   meta=False)

            va = wrf.getvar(f, 'va', file_timeidx)
            interp_va = wrf.vinterp(f,
                                   field=va,
                                   vert_coord="pressure",
                                   interp_levels=interp_levs,
                                   extrapolate=False,
                                   field_type=None,
                                   log_p=True,
                                   timeidx=file_timeidx,
                                   meta=False)

            QVAPOR = wrf.getvar(f, 'QVAPOR', file_timeidx) # water vapor mixing ratio in kg kg-1
            rv = wrf.vinterp(f,
                       field=QVAPOR,
                       vert_coord="pressure",
                       interp_levels=interp_levs,
                       extrapolate=False,
                       field_type=None,
                       log_p=True,
                       timeidx=file_timeidx,
                       meta=False)

            # calculate specific humidity from qvapor
            var_q = np.divide(rv, (1+rv)) # specific humidity kg kg-1

            # integrate water vapor transport and water vapor
            pressure = np.reshape(interp_levs, (20, 1, 1))*100 # convert from hPa to Pa
            g = -9.81 # gravity constant
            ivtu = trapz(interp_ua*var_q, pressure, axis=0)/g
            ivtv = trapz(interp_va*var_q, pressure, axis=0)/g
            iwv = trapz(var_q, pressure, axis=0)/g

            # get current time step
            da_time.append(ua.Time.values)
            
            # put values into preassigned arrays
            ivtu_final.append(ivtu)
            ivtv_final.append(ivtv)
            iwv_final.append(iwv)
            
        f.close()
    
    # get lats and lons
    wrflats = ua['XLAT'].isel(west_east=0).values
    wrflons = ua['XLONG'].isel(south_north=0).values

    # put into a dataset
    var_dict = {'ivtu': (['time', 'lat', 'lon'], ivtu_final),
                'ivtv': (['time', 'lat', 'lon'], ivtv_final), 
                'iwv': (['time', 'lat', 'lon'], iwv_final)}
    ds = xr.Dataset(var_dict,
                    coords={'time': (['time'], da_time),
                            'lat': (['lat'], wrflats),
                            'lon': (['lon'], wrflons)})

    return ds

def preprocess_prec(filenames):
    """preprocess_prec
    
    Returns a ds object (xarray) including total precipitation accumulated, snow height/depth, snow accumulated
    
    Parameters
    ----------
    filenames : list
        list of wrf filenames to process
  
    Returns
    -------
    ds : ds object
        includes variables prec accumulated, snow height, snow accumulated
    
    """
    # arrays to append data
    prec_final = []
    snowh_final = []
    snow_final = []
    da_time = []

    for i, wrfin in enumerate(filenames):
        f = nc.Dataset(wrfin, 'r')
        # get the number of times in the current file
        ntimes = f.variables['Times'].shape[0]

        # now loop through the total number of times for this file and append data:
        for file_timeidx in range(ntimes):
            print(wrfin, 'timeidx is ', file_timeidx)
            rainc = wrf.getvar(f, 'RAINC', file_timeidx) # units: total mm accumulated
            rainnc = wrf.getvar(f, 'RAINNC', file_timeidx) # units: total mm accumulated
            # add rainc + rainnc for total mm accumulated since model initialization (3/31/2012 in this case)
            prec = rainc + rainnc
            snowh = wrf.getvar(f, 'SNOWH', file_timeidx) # units: m depth
            snow = wrf.getvar(f, 'SNOW', file_timeidx) # units: kg m^2 = mm

            # get current time
            da_time.append(snow.Time.values)

            # put values into preassigned arrays
            prec_final.append(prec)
            snowh_final.append(snowh)
            snow_final.append(snow)

        f.close()

    # get lats and lons
    wrflats = snow['XLAT'].isel(west_east=0).values
    wrflons = snow['XLONG'].isel(south_north=0).values
    
    # convert lats/lons to 4-byte floats (this alleviates striping issue)
    wrflats = np.float32(wrflats)
    wrflons = np.float32(wrflons)

    # put into a dataset
    var_dict = {'prec': (['time', 'lat', 'lon'], prec_final),
                'snowh': (['time', 'lat', 'lon'], snowh_final), 
                'snow': (['time', 'lat', 'lon'], snow_final)}
    ds = xr.Dataset(var_dict,
                    coords={'time': (['time'], da_time),
                            'lat': (['lat'], wrflats),
                            'lon': (['lon'], wrflons)})

    return ds

def preprocess_pressure_lev_var(filenames, var1_name, var2_name, levs):
    """preprocess_pressure_lev_var
    
    Returns a ds object (xarray) including specified variable to specified vertical coordinate
    
    Parameters
    ----------
    filenames : list
        list of wrf filenames to process
    var1_name : string
        string of variable you want to interpolate
    var2_name : string
        string of vertical variable you want to interpolate var1_name to 
    levs : vertical levels of var2_name you want to interpolate var1_name to
  
    Returns
    -------
    ds : ds object
        includes variable specified in var1_name interpolated to levels of var2_name specified
    
    Example
    -------
    # Preprocess Geopotential Heights to 250 and 500 hPa
    ds = preprocess_pressure_lev_var(filenames, 'z', 'pressure', [250., 500.])
    
    """
    # arrays to append data
    var_final = []
    da_time = []

    for i, wrfin in enumerate(filenames):
        f = nc.Dataset(wrfin, 'r')
        # get the number of times in the current file
        ntimes = f.variables['Times'].shape[0]

        # now loop through the total number of times for this file and append data:
        for file_timeidx in range(ntimes):
            print(wrfin, 'timeidx is ', file_timeidx)
            var1 = wrf.getvar(f, var1_name, file_timeidx)
            var2 = wrf.getvar(f, var2_name, file_timeidx)

            # get lats, lons, and time
            da_time.append(var1.Time.values)

            # interpolate z to pressure levels
            interp_var = wrf.interplevel(var1, var2, levs)
            # put values into preassigned arrays
            var_final.append(interp_var.values)

        f.close()

    var_final = np.stack(var_final, axis=0)
    print(var_final.shape)
    # get lats and lons
    wrflats = var1['XLAT'].isel(west_east=0).values
    wrflons = var1['XLONG'].isel(south_north=0).values
    # convert lats/lons to 4-byte floats (this alleviates striping issue)
    wrflats = np.float32(wrflats)
    wrflons = np.float32(wrflons)

    # put into a dataset
    if len(levs) < 2: 
        var_dict = {var1_name: (['time', 'lat', 'lon'], var_final)}
        ds = xr.Dataset(var_dict,
                        coords={'time': (['time'], da_time),
                                'lat': (['lat'], wrflats),
                                'lon': (['lon'], wrflons)})
    else: 
        var_dict = {var1_name: (['time', 'lev', 'lat', 'lon'], var_final)}
        ds = xr.Dataset(var_dict,
                        coords={'time': (['time'], da_time),
                                'lev': (['lev'], levs),
                                'lat': (['lat'], wrflats),
                                'lon': (['lon'], wrflons)})
    
    return ds

def preprocess_SR(filenames):
    """preprocess_prec
    
    Returns a ds object (xarray) including the fraction of frozen precipitation (unitless)
    
    Parameters
    ----------
    filenames : list
        list of wrf filenames to process
  
    Returns
    -------
    ds : ds object
        includes variables SR "fraction of frozen preciptiation"
    
    """
    # arrays to append data
    sr_final = []
    da_time = []

    for i, wrfin in enumerate(filenames):
        f = nc.Dataset(wrfin, 'r')
        # get the number of times in the current file
        ntimes = f.variables['Times'].shape[0]

        # now loop through the total number of times for this file and append data:
        for file_timeidx in range(ntimes):
            print(wrfin, 'timeidx is ', file_timeidx)
            sr = wrf.getvar(f, 'SR', file_timeidx) # units: unitless

            # get current time
            da_time.append(sr.Time.values)

            # put values into preassigned arrays
            sr_final.append(sr)

        f.close()

    # get lats and lons
    wrflats = sr['XLAT'].isel(west_east=0).values
    wrflons = sr['XLONG'].isel(south_north=0).values
    
    # convert lats/lons to 4-byte floats (this alleviates striping issue)
    wrflats = np.float32(wrflats)
    wrflons = np.float32(wrflons)

    # put into a dataset
    var_dict = {'sr': (['time', 'lat', 'lon'], sr_final)}
    ds = xr.Dataset(var_dict,
                    coords={'time': (['time'], da_time),
                            'lat': (['lat'], wrflats),
                            'lon': (['lon'], wrflons)})

    return ds

def preprocess_wrf_cfcompliant(ds):
    # clean up protocol adapted from http://gallery.pangeo.io/repos/NCAR/notebook-gallery/notebooks/Run-Anywhere/WRF/wrf_ex.html

    ## clean up time coordinate
    # da_time = ds['XTIME'].values # get array of times in datetime64 format
    # ds = ds.rename({'Time':'time'}) # rename our dimension Time to time to match conventions
    # ds = ds.assign(time=da_time) # assign time coordinate xtime values
    # ds = ds.drop('XTIME') # drop XTIME coord

    ## clean up lat lon coords
    # XLAT and XLONG – these are our Latitude and Longitude values
    # XLAT_U and XLONG_U – Lat and Long with a staggered west-east grid
    # XLAT_V and XLONG_V – Lat and Long with a staggered north-south grid
    # create lat lon coordinates from XLAT and XLONG
    lats = ds['XLAT'].isel(west_east=0).values
    lons = ds['XLONG'].isel(south_north=0).values
    ds = ds.assign_coords(lat=lats, lon=lons)
    # # add land and lake masks to coords
    # ds = ds.assign_coords(landmask=ds['LANDMASK'], lakemask=ds['LAKEMASK'])
    # rename south_north to y and west_east to x
    ds = ds.rename({'south_north':'lat', 'west_east':'lon'})
    # drop the XLAT and XLONG coordinates
    ds = ds.drop(['XLAT', 'XLONG'])
    
    return ds

def wrf_3hr_to_daily(output_varname, domain):
    '''opens preprocessed 3hrly wrf files and computes daily files'''
    start_yr = 1979
    end_yr = 2015

    path_to_data = '/scratch1/08540/dlnash/data/wrf_6km/'
    path_to_out = '/work2/08540/dlnash/frontera/data/wrf_preprocessed_data/wrf_6km/'

    ## pull wrflats and wrflons from first file
    fname = path_to_out + '{0}/{1}/3hr/tmp_1979.nc'.format(domain, output_varname)
    tmp = xr.open_dataset(fname)

    ## assign those lats to the other ds when you loop
    wrflats = tmp.lat.values
    wrflons = tmp.lon.values

    ## Loop through all above years
    for yr in np.arange(start_yr, end_yr+1):
        filename = path_to_out + '{0}/{1}/3hr/tmp_{2}.nc'.format(domain, output_varname, str(yr))

        tmp = xr.open_dataset(filename, chunks={'time': 360})
        
        if output_varname == 'prec':
            tmp = tmp.sel(time=tmp.time.dt.hour == 0)

            # # calculate mm per day rain values
            # # rain at 00:00 utc next day - rain at 00:00 UTC current day
            wrf = tmp.shift(time=-1) - tmp # if in xarray

            ## bc of spin up and how the data was generated set April 1 values to nan
            if yr > 1979:
                wrf = xr.where((wrf.time.dt.month == 4) & (wrf.time.dt.day == 1), np.nan, wrf)
                
            # # set negative snow values to 0
            # wrf.snow.where((wrf.snow < 0), 0, wrf.snow)
                
        else:
            # resample to daily
            wrf = tmp.resample(time="1D").mean('time')

        # reassign lats/lons
        wrf = wrf.assign_coords({"lon": wrflons, "lat": wrflats})
        if yr == 1979:
            wrf = wrf
        elif yr == 2015:
            wrf = wrf
        else:
            wrf = wrf.sel(time=slice('{0}-01-01'.format(yr), '{0}-12-31'.format(yr)))

        # write to netCDF
        fname = os.path.join(path_to_out, '{0}/{1}/daily/out.wrf6km.{1}.daily_{2}.nc').format(domain, output_varname, str(yr))
        wrf.to_netcdf(path=fname, mode = 'w', format='NETCDF4')
        
def select_single_coord_WRF(filenames, varlst, slat, slon, dates):
    '''function to select a all levels of variable list at a single lat/lon coord from WRF data
    
    Returns a ds object (xarray) including variables indicated
    
    Parameters
    ----------
    filenames : list
        list of wrf filenames to process
        
    varlst : list
        list of variables to extract
    
    slat : float
        latitude of single coordinate interested in selecting
    
    slon : float
        longitude of single coordinate interested in selecting
        
    dates : array
        array of dates in datetime64 format for which you want to subset WRF to
  
    Returns
    -------
    ds : ds object
        includes variables at specified coordinate for all levels
    
    '''
    # arrays to append data
    darray_final = []

    for i, wrfin in enumerate(filenames):
        f = nc.Dataset(wrfin)
        # get the number of times in the current file
        ntimes = f.variables['Times'].shape[0]
        
        # get array of all times
        var = wrf.getvar(f, 'RAINC', wrf.ALL_TIMES)
        times = var.Time.values

        d = {'wrf_times': times}
        if ntimes > 1:
            tmp = pd.DataFrame(d, index=list(np.arange(len(times))))
        elif ntimes == 1:
            tmp = pd.DataFrame(d, index=list(np.arange(len([times]))))
        tmp.index = tmp['wrf_times']

        # find idx in times where dates == times
        t = tmp.index.isin(dates)
        idx = [i for i, x in enumerate(t) if x]

        # now loop through the total number of times for this file and append data:
        for j, file_timeidx in enumerate(idx):
            print(wrfin, 'timeidx is ', file_timeidx)
            
            # make cache for current time step
            my_cache = wrf.extract_vars(f, file_timeidx, ("P", "PSFC", "PB", "PH", "PHB",
                                                          "T", "QVAPOR", "HGT", "U", "V",
                                                          "W"))
            # loop through variables
            darrays = []
            for i, var in enumerate(varlst):
                v = wrf.getvar(f, var, file_timeidx, cache=my_cache)
                # remove projection from attributes
                del v.attrs['projection']
                darrays.append(v)

            # merge list of data arrays to single dataset for current timestep
            ds = xr.merge(darrays)

            # convert ds into cf compliant so we can select given lat/lon coord
            ds = preprocess_wrf_cfcompliant(ds)
            ds = ds.sel(lat=slat, lon=slon, method="nearest")
            
            # append current time step ds to list
            darray_final.append(ds)
        
        # close current file
        f.close()
    
    # concat all time steps into single ds
    ds = xr.concat(darray_final, dim='Time')

    # calculate specific humidity from mixing ratio (units = kg kg-1)
    ds = ds.assign(q=lambda ds: np.divide(ds.QVAPOR, (1+ds.QVAPOR)))

    # calculate windspeed
    ds = ds.assign(wspd=lambda ds: np.sqrt(ds.ua**2 + ds.va**2))
    
    return ds


def calculate_WRF_vertical_cross(filenames, varlst, df, startlat, startlon, endlat, endlon):
    '''function to select a all levels of variable list from WRF data using dynamic lat/;on coord pairs
    
    Returns a ds object (xarray) including variables indicated
    
    Parameters
    ----------
    filenames : list
        list of wrf filenames to process
        
    varlst : list
        list of variables to extract
        
    df : pandas df
        pandas dataframe with dates interested in and lat/lon coord pairs
        
    startlat : string
        name of pandas df column for starting latitude coordinate
        
    startlon : string
        name of pandas df column for starting longitude coordinate
        
    endlat : string
        name of pandas df column for ending latitude coordinate
        
    endlon : string
        name of pandas df column for ending longitude coordinate
  
    Returns
    -------
    list : list of ds objects
        includes ds objects with cross sections at indicated start/end coordinate pairs
    
    '''
    # arrays to append data
    ds_lst = []
    # levels to interpolate for cross section
    crosslevs = np.arange(0, 12050, 100)

    for i, wrfin in enumerate(filenames):
        f = nc.Dataset(wrfin, 'r')
        # get the number of times in the current file
        ntimes = f.variables['Times'].shape[0]
        
        # get array of all times
        var = wrf.getvar(f, 'RAINC', wrf.ALL_TIMES)
        times = var.Time.values

        d = {'wrf_times': times}
        if ntimes > 1:
            tmp = pd.DataFrame(d, index=list(np.arange(len(times))))
        elif ntimes == 1:
            tmp = pd.DataFrame(d, index=list(np.arange(len([times]))))
        tmp.index = tmp['wrf_times']

        # find idx in times where dates == times
        dates = df.time.values
        t = tmp.index.isin(dates)
        idx = [i for i, x in enumerate(t) if x]

        # now loop through the total number of times for this file and append data:
        for j, file_timeidx in enumerate(idx):
            print(wrfin, 'timeidx is ', file_timeidx)
            
            # Set the start point and end point for the cross section
            start_point = wrf.CoordPair(lat=df[startlat].iloc[j], lon=df[startlon].iloc[j])
            end_point = wrf.CoordPair(lat=df[endlat].iloc[j], lon=df[endlon].iloc[j])
            print(start_point, end_point)
            
            # make cache for current time step
            my_cache = wrf.extract_vars(f, file_timeidx, ("P", "PSFC", "PB", "PH", "PHB",
                                                          "T", "QVAPOR", "HGT", "U", "V",
                                                          "W"))
            # loop through variables
            cross = []
            for i, var in enumerate(varlst):
                v = wrf.getvar(f, var, file_timeidx, cache=my_cache)
                z = wrf.getvar(f, 'z', file_timeidx, cache=my_cache)
                # use linear z for interpolation
                V = 10**(v/10.)
                cross_v = wrf.vertcross(V, z, levels=crosslevs, wrfin=f, start_point=start_point,
                          end_point=end_point, latlon=True, meta=True, cache=my_cache)
                
                # convert back to dBz after interpolation
                new_v = 10.0 * np.log10(cross_v)

                cross.append(new_v)
                

            # merge list of data arrays to single dataset for current timestep
            ds = xr.merge(cross)
            # calculate specific humidity from mixing ratio (units: kg kg-1)
            ds = ds.assign(q=lambda ds: np.divide(ds.QVAPOR_cross, (1+ds.QVAPOR_cross)))
            # # calculate windspeed
            # ds = ds.assign(wspd=lambda ds: np.sqrt(ds.ua**2 + ds.va**2))
            
            # compute vertical moisture flux (units: m s-1*kg kg-1)
            uq = ds.ua_cross*ds.q
            vq = ds.va_cross*ds.q
            ds = ds.assign(wvf=lambda ds: np.sqrt(uq**2 + vq**2))
            
            # Get the terrain heights along the cross section line
            terrain = wrf.getvar(f, 'ter', file_timeidx, cache=my_cache)
            ter_line = wrf.interpline(terrain, wrfin=f, start_point=start_point, end_point=end_point)
            dims = ('cross_line_idx')
            ds = ds.assign(ter=lambda ds:(dims, ter_line.data))
            
            # append current time step ds to list
            ds_lst.append(ds)
        
        # close current file
        f.close()
    
    return ds_lst


def calculate_WRF_vertical_cross_static(filenames, varlst, dates, startlat, startlon, endlat, endlon):
    '''function to select a all levels of variable list from WRF data using static lat/lon coord pairs
    
    Returns a ds object (xarray) including variables indicated
    
    Parameters
    ----------
    filenames : list
        list of wrf filenames to process
        
    varlst : list
        list of variables to extract
        
    dates : array
        array of datetimes interested in and lat/lon coord pairs
        
    startlat : float
        starting latitude coordinate
        
    startlon : float
        starting longitude coordinate
        
    endlat : float
        ending latitude coordinate
        
    endlon : float
        ending longitude coordinate
  
    Returns
    -------
    list : list of ds objects
        includes ds objects with cross sections at indicated start/end coordinate pairs
    
    '''
    # arrays to append data
    ds_lst = []
    # levels to interpolate for cross section
    crosslevs = np.arange(0, 12050, 100)

    for i, wrfin in enumerate(filenames):
        f = nc.Dataset(wrfin, 'r')
        # get the number of times in the current file
        ntimes = f.variables['Times'].shape[0]
        
        # get array of all times
        var = wrf.getvar(f, 'RAINC', wrf.ALL_TIMES)
        times = var.Time.values

        d = {'wrf_times': times}
        if ntimes > 1:
            tmp = pd.DataFrame(d, index=list(np.arange(len(times))))
        elif ntimes == 1:
            tmp = pd.DataFrame(d, index=list(np.arange(len([times]))))
        tmp.index = tmp['wrf_times']

        # find idx in times where dates == times
        t = tmp.index.isin(dates)
        idx = [i for i, x in enumerate(t) if x]

        # now loop through the total number of times for this file and append data:
        for j, file_timeidx in enumerate(idx):
            print(wrfin, 'timeidx is ', file_timeidx)
            
            # Set the start point and end point for the cross section
            start_point = wrf.CoordPair(lat=startlat, lon=startlon)
            end_point = wrf.CoordPair(lat=endlat, lon=endlon)
            print(start_point, end_point)
            
            # make cache for current time step
            my_cache = wrf.extract_vars(f, file_timeidx, ("P", "PSFC", "PB", "PH", "PHB",
                                                          "T", "QVAPOR", "HGT", "U", "V",
                                                          "W"))
            # loop through variables
            cross = []
            for i, var in enumerate(varlst):
                v = wrf.getvar(f, var, file_timeidx, cache=my_cache)
                z = wrf.getvar(f, 'z', file_timeidx, cache=my_cache)
                # use linear z for interpolation
                V = 10**(v/10.)
                cross_v = wrf.vertcross(V, z, levels=crosslevs, wrfin=f, start_point=start_point,
                          end_point=end_point, latlon=True, meta=True, cache=my_cache)
                
                # convert back to dBz after interpolation
                new_v = 10.0 * np.log10(cross_v)

                cross.append(new_v)
                

            # merge list of data arrays to single dataset for current timestep
            ds = xr.merge(cross)
            # calculate specific humidity from mixing ratio (units: kg kg-1)
            ds = ds.assign(q=lambda ds: np.divide(ds.QVAPOR_cross, (1+ds.QVAPOR_cross)))
            # # calculate windspeed
            # ds = ds.assign(wspd=lambda ds: np.sqrt(ds.ua**2 + ds.va**2))
            
            # compute vertical moisture flux (units: m s-1*kg kg-1)
            uq = ds.ua_cross*ds.q
            vq = ds.va_cross*ds.q
            ds = ds.assign(wvf=lambda ds: np.sqrt(uq**2 + vq**2))
            
            # Get the terrain heights along the cross section line
            terrain = wrf.getvar(f, 'ter', file_timeidx, cache=my_cache)
            ter_line = wrf.interpline(terrain, wrfin=f, start_point=start_point, end_point=end_point)
            dims = ('cross_line_idx')
            ds = ds.assign(ter=lambda ds:(dims, ter_line.data))
            
            # append current time step ds to list
            ds_lst.append(ds)
        
        # close current file
        f.close()
    
    return ds_lst