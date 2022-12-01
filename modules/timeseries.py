"""
Collection of functions used to analyze time series (1-dimensional) data

"""

# Imports

import os, sys
import numpy as np
import xarray as xr
import pandas as pd
from scipy.stats import ttest_1samp, t, pearsonr


# Define functions

def select_months_ds(ds, mon_s, mon_e):    
    # Select months
    if mon_s > mon_e:
        idx = (ds.time.dt.month >= mon_s) | (ds.time.dt.month <= mon_e)
    else:
        idx = (ds.time.dt.month >= mon_s) & (ds.time.dt.month <= mon_e)
    ds = ds.sel(time=idx)
    return ds

def combine_ar_ds_df(ds, df):
    # Combine AR Cat data w/ WRF data
    # Add ar time series to the WRF dataset
    ds['ar'] = ('time', df.AR_CAT)
    ds = ds.set_coords('ar')

    ds['trackID'] = ('time', df.kidmap)
    ds = ds.set_coords('trackID')
    idx = (ds.ar >= 1)
    # select AR days
    ds_ar = ds.sel(time=idx)
    return ds_ar

def persistence(x):
    """A function that tags persistent events
    
    Given a binary time series `x`, this function tags persistent events, 
    (x eq 1). A sequence of consecutive x=1 values is tagged as a single event.
    
    Parameters
    ----------
    x : array_like
        Binary time series of x=0 or x=1
        
    Returns
    -------
    event_id : array_like, int
        Array of the same length as x that contains tags for each event
    duration = array_like, int
        Values represent the number of observations for for each event. 1D array with a length 
        
        
    Example
    -------
    Given:
        x          = [0,1,1,0,0,1,0,0,1,0,0,1,1,0,1,1,1,0]
    Returns:
        event_id   = [0,1,1,0,0,2,0,0,3,0,0,4,4,0,5,5,5,0]
        duration   = [  2,      1,    1,    2,    3,     ]
        
    References ******
    ----------
    Adapted from Charles Jones' IDL subroutine

    """
    
    # number of observations in x
    ntot = len(x)

    # Loop to tag persistent events
    event_id = np.zeros(ntot, dtype=int)
    tag = 1
    test = 0    
    for k in range(ntot):       
        # test for active event
        if (x[k] == 1):
            test = 1
            event_id[k] = tag
        # test for end of event    
        elif (x[k] == 0) and (test == 1):
            tag += 1
            test = 0

    # Loop to find duration of each event
    nevents = event_id.max()       # Total number of tagged events
    duration = np.empty(nevents)
    for k in range(nevents):
        # eid = event id
        eid = k+1
        # find the event indices in x
        idx = np.where(event_id == eid)
        # find the length of the event and store in duration arr
        duration[k] = len(idx[0])
    
    return event_id, nevents, duration


#### CLEAN UP AND ADD DOC STRINGS ###
def transition_matrix(x, states):

    morder = 1                 # model order (default=1) ** make optional keyword param
    nt = len(x)-1              # number of transitions
    n = len(states)            # number of states
    transc = np.zeros((n,n), dtype=int)  # Matrix of transitions counts (initialize counts as 0)

    # Loop to count transitions
    for t in range(nt):
        i = np.where(states == x[t])[0]      # i = indice of current state
        j = np.where(states == x[t+1])[0]    # j = indice of next state    
        transc[i,j] += 1                    # add a transition count for s[i] to s[j]    

    # Compute marginal totals (sum across rows in M)
    margin = np.sum(transc, axis=1)

    # Calculate probabilities (divide each row in M by its marginal total)
    probs = transc / margin[:,np.newaxis]
    
    return transc, probs


def select_months_df(df, mon_s, mon_e):
    # Select months
    if mon_s > mon_e:
        idx = (df.index.month >= mon_s) | (df.index.month <= mon_e)
    else:
        idx = (df.index.month >= mon_s) & (df.index.month <= mon_e)

    df = df.loc[idx]
    
    return df 

def create_list_all_dates(start_date, end_date, date_lst):
    """
    From a list of dates, create an array of zeros and ones 
    where ones are where the conditions are true 
    
    Parameters
    ----------
    start_date : str
        start date of list
    end_date : str
        end date of list
    date_lst : array of datetimes
        list of datetimes where conditions are true
    
    Returns
    -------
    arr_allDays : numpy array
        numpy array where the dates where the condition is true is ones
    
    """
    # date array with all days
    dates_allDays = pd.date_range(start=start_date, end=end_date, freq='1D')
    arr_allDays = np.zeros(len(dates_allDays), dtype=np.float)

    # Loop over ar days and match to ar_full 
    for i, date in enumerate(date_lst):
        idx = np.where(dates_allDays == date)
        arr_allDays[idx] = 1
    
    return arr_allDays