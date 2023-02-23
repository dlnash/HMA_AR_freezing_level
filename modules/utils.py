"""
Filename:    utils.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: Functions for loading and manipulating large datasets
"""
import os, sys
import yaml
import xarray as xr
import pandas as pd
import numpy as np
from glob import glob
import math
from datetime import timedelta, datetime

def slope_2pts(x1, y1, x2, y2):
    deltay = (y2-y1)
    deltax = (x2-x1)
    slope = deltay/deltax
    b = y1 - slope*x1
    
    return b, slope

def find_intersection_two_lines(m1, b1, m2, b2):

    x = (b2 - b1)/(m1-m2)
    y = m1*x + b1
    
    return [x, y]

def find_perpindicular_line(x1, y1, x2, y2, newx):
    
    # first calculate slope and y intercept of existing line
    b, slope = slope_2pts(x1, y1, x2, y2)
    # now get the slope of the perpindicular line
    newslope = -1./slope
    # now find the corresponding y point to newx on old line
    newy = slope*newx + b
    # now find the new line b
    newb = newy - newslope*newx
    # x3 = 65.0 # furthest longitude west in outer domain
    x3 = 68.5 # furthest longitude west in inner domain
    y3 = newslope*x3 + newb
    # hlat, hlon, tlat, tlon
    line = [y3, x3, newy, newx]
    eq = [newslope, newb]
    
    return line, eq

def find_parallel_line(x1, y1, x2, y2, d, x3, x4):
    # first calculate slope and y intercept of existing line
    b, slope = slope_2pts(x1, y1, x2, y2)
    eq1 = [slope, b]
    # using the new y-intercept (b - distance from original line)
    newb = (b+d)
    # find the associated y-value given x3
    y3 = slope*x3 + newb
    
    # find the associated y-value given x4
    y4 = slope*x4 + newb
    line = [y3, x3, y4, x4]
    eq2 = [slope, newb]
    
    return line, eq1, eq2


    
    
    

def calc_dy(slope):
    dy = math.sqrt(3**2/(slope**2+1))
    return dy

def calc_dx(slope, dy):
    dx = -slope*dy
    return dx

def check_mkdir(filename):
    '''Checks if directory exists and if not, makes that directory'''
    if not os.path.exists(os.path.dirname(filename)):
        try:
            os.makedirs(os.path.dirname(filename))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
                
class load_era5_ds:
    def __init__(self, path_to_data, pathvar, bbox, anom=True, lev=None):
        if anom == True:
            self.filepath_pattern = path_to_data + 'ERA5/{0}/anomalies/daily_filtered_anomalies_{0}_*.nc'.format(pathvar)
        elif anom == False: 
            self.filepath_pattern = path_to_data + 'ERA5/{0}/daily/out.era5_*.nc'.format(pathvar)
        
        yaml_doc = '../data/load_era5.yaml'
        config = yaml.load(open(yaml_doc), Loader=yaml.SafeLoader)
        self.rename_dict = config[pathvar]
        
        self.pathvar = pathvar
        self.bbox = bbox
        self.lev = lev
        
    def preprocess_era5(self, ds):
        lonmin, lonmax, latmin, latmax = self.bbox
        if self.pathvar == 'huvq':
            subset = ds.sel(latitude=slice(latmax, latmin), longitude=slice(lonmin, lonmax), level=self.lev)
            subset = subset.rename(self.rename_dict)
        elif self.pathvar == 'prec':
            subset = ds.sel(latitude=slice(latmax, latmin), longitude=slice(lonmin, lonmax))
            subset = subset.rename(self.rename_dict)
        elif self.pathvar == 'iwv':
            subset = ds.sel(latitude=slice(latmax, latmin), longitude=slice(lonmin, lonmax))
            subset = subset.rename(self.rename_dict)
            subset = subset.drop(['tcrw', 'tcsw', 'tcw'])
        elif self.pathvar == 'ivt':
            subset = ds.sel(latitude=slice(latmax, latmin), longitude=slice(lonmin, lonmax))
            subset = subset.rename(self.rename_dict)
        
        return subset
            
    def load_ds(self):
        f = xr.open_mfdataset(self.filepath_pattern, combine='nested', concat_dim='time', preprocess=self.preprocess_era5)
        # normalize times to 00 UTC
        f['time'] = f.indexes['time'].normalize()
        return f

    
            



def leap_year(year, calendar='standard'):
    """Determine if year is a leap year"""
    leap = False
    if ((calendar in ['standard', 'gregorian',
        'proleptic_gregorian', 'julian']) and
        (year % 4 == 0)):
        leap = True
        if ((calendar == 'proleptic_gregorian') and
            (year % 100 == 0) and
            (year % 400 != 0)):
            leap = False
        elif ((calendar in ['standard', 'gregorian']) and
                 (year % 100 == 0) and (year % 400 != 0) and
                 (year < 1583)):
            leap = False
    return leap

def get_dpm(time, calendar='standard'):
    """
    return a array of days per month corresponding to the months provided in `months`
    """
    month_length = np.zeros(len(time), dtype=np.int)
    
    dpm = {'noleap': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '365_day': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'standard': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'proleptic_gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'all_leap': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '366_day': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '360_day': [0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}
    
    cal_days = dpm[calendar]

    for i, (month, year) in enumerate(zip(time.month, time.year)):
        month_length[i] = cal_days[month]
        if leap_year(year, calendar=calendar):
            month_length[i] += 1
    return month_length

# Wrap it into a simple function
def season_mean(ds, calendar='standard'):
    # Make a DataArray of season/year groups
    year_season = xr.DataArray(ds.time.to_index().to_period(freq='Q-NOV').to_timestamp(how='E'),
                               coords=[ds.time], name='year_season')

    # Make a DataArray with the number of days in each month, size = len(time)
    month_length = xr.DataArray(get_dpm(ds.time.to_index(), calendar=calendar),
                                coords=[ds.time], name='month_length')
    # Calculate the weights by grouping by 'time.season'
    weights = month_length.groupby('time.season') / month_length.groupby('time.season').sum()

    # Test that the sum of the weights for each season is 1.0
    np.testing.assert_allclose(weights.groupby('time.season').sum().values, np.ones(4))

    # Calculate the weighted average
    return (ds * weights).groupby('time.season').sum(dim='time')


def add_days_to_date(date, days):
    """Add days to a date and return the date.
    
    Args: 
        date (string): Date string in YYYY-MM-DD format. 
        days (int): Number of days to add to date
    
    Returns: 
        date (date): Date in YYYY-MM-DD with X days added. 
    """
    
    added_date = pd.to_datetime(date) + timedelta(days=days)
    added_date = added_date.strftime("%Y-%m-%d")

    return added_date

def subtract_days_from_date(date, days, normalized=False):
    """Subtract days from a date and return the date.
    
    Args: 
        date (string): Date string in YYYY-MM-DD format. 
        days (int): Number of days to subtract from date
    
    Returns: 
        date (date): Date in YYYY-MM-DD with X days subtracted. 
    """
    
    subtracted_date = pd.to_datetime(date) - timedelta(days=days)
    
    if normalized == True:
        subtracted_date = subtracted_date.strftime("%Y-%m-%d")
    else:
        # convert the timestamp to a datetime object in the local timezone
        # subtracted_date = datetime.fromtimestamp(subtracted_date)
        subtracted_date = subtracted_date

    return subtracted_date