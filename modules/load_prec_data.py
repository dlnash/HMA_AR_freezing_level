"""
Filename:    plotter.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: Functions for plotting
"""

# Import Python modules

import os, sys


class LoadPrec:

    def __init__(self, name, age):
        self.name = name
        self.age = age
        self.path_to_data = '/home/nash/DATA/data/' 
        
    def load_ERA5_prec():
        rename_dict_prec = {'mtpr': 'prec', 
                            'latitude': 'lat',
                            'longitude': 'lon'}

        filepath_pattern = self.path_to_data + 'ERA5/prec/6hr/era5_hma_025dg_6hr_prec_*.nc'
        ds = xr.open_mfdataset(filepath_pattern, combine='by_coords')
        ds = ds.rename(rename_dict_prec)
        ds['time'] = ds.indexes['time'].normalize()
        ds = ds.assign(prec=lambda ds: ds.prec*(60*60*6)) # convert to mm accumulated per 6-hours
        
        return ds
    
    def load_WRF_6km_prec():
        names = glob.glob(self.path_to_data + 'wrf_hasia/prec/3hr/tmp_*.nc')
        filenames = sorted(names)

        tmp = xr.open_mfdataset(filenames, combine='by_coords', join='override')
        tmp = tmp.sel(time=tmp.time.dt.hour == 0)

        # # calculate mm per day rain values
        # # rain at 00:00 utc next day - rain at 00:00 UTC current day
        ds = tmp.shift(time=-1) - tmp # if in xarray

        # Trim date range
        idx = slice(start_date, end_date)
        ds = ds.sel(time=idx)

        # select only months we are interested in
        ds = select_months_ds(ds, mon_s, mon_e)

        return ds
    
    def load_GPM_prec():
        def preprocess(ds):
            '''keep only the selected lats and lons'''
            return ds.sel(lat=slice(latmin, latmax), lon=slice(lonmin, lonmax))

        filename_pattern1 = path_to_data + 'IMERGV06B/3B-DAY.MS.MRG.3IMERG.20*.nc4'
        gpm = xr.open_mfdataset(filename_pattern1, engine='netcdf4', concat_dim='time', combine='by_coords',
                                   preprocess=preprocess)

        print('ds size in GB {:0.2f}\n'.format(gpm.nbytes / 1e9))
        gpm = gpm.transpose('lat', 'lon', 'time')
        gpm = gpm.rename({'precipitationCal': 'prec'})
        gpm