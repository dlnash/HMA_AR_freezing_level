"""
Filename:    metrics.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: Metrics used to evaluate WRF output
"""

import os, sys
import xarray as xr
import pandas as pd

def calculate_weighted_ts_mean(ds, domains):
    '''
    This function takes the given xarray dataset and 
    domain bounds (lonmin, lonmax, latmin, latmax) 
    and calculates the weighted mean within each domain.
    
    Parameters
    ----------
    ds : xarray dataset
        xarray dataset object
    domains : list
        list of domain bounds - each domain contains (lonmin, lonmax, latmin, latmax) 
  
    Returns
    -------
    df_lst : list of pandas dataframes
        list of dataframes (in the order of the domains) that shows the precipitation rate and cumulative precipitation variables for the given domain
        
    '''
    df_lst = []
    for i, bnds in enumerate(domains):
        tmp = ds.sel(lat=slice(bnds[2], bnds[3]), lon=slice(bnds[0], bnds[1]))
        weights = np.cos(np.deg2rad(tmp.lat))
        tmp_weighted = tmp.weighted(weights)
        weighted_mean = tmp_weighted.mean(("lon", "lat"))

        df = weighted_mean.to_dataframe()
        df['prec_cumsum'] = df.prec.cumsum()
        df['snow_cumsum'] = df.snow.cumsum()

        df_lst.append(df)
        
    return df_lst