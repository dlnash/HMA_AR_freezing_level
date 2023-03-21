"""
Filename:    generate_df_ivt_prec_freezing.py
Author:      Deanna Nash, dlnash@ucsb.edu
Description: Pull IVT, precip, freezing level and landslide info for each AR event during DJF
"""

# Standard Python modules
import os, sys
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta
import metpy.calc as mpcalc
from metpy.units import units

# import personal modules

# Path to modules
sys.path.append('../modules')

# Import my modules
from ar_funcs import get_topo_mask
from timeseries import select_months_ds, select_months_df

# Set up paths
server = 'great'
path_to_data = '/home/nash/DATA/data/'                                      # project data -- read only
path_to_out  = '../out/'       # output files (numerical results, intermediate datafiles) -- read & write
path_to_figs = '../figs/'      # figures

### Get list of AR dates and trackIDs when an AR crosses 1000 m elevation threshold in HMA

# identify ARs using single bound box with elevation mask during DJF
bbox = [20, 40, 65, 97] # HMA region
elev_thres = 1000.
ssn = 'DJF'

if ssn == 'DJF':
    start_date = '1979-12-01 0:00'
    end_date = '2015-02-28 18:00'
    start_mon = 12
    end_mon = 2
if ssn == 'MAM':
    start_date = '1980-03-01 0:00'
    end_date = '2014-05-31 18:00'
    start_mon = 3
    end_mon = 5

# open ds
filename =  path_to_data + 'ar_catalog/globalARcatalog_ERA-Interim_1979-2019_v3.0.nc'
ds = xr.open_dataset(filename, chunks={'time': 1460}, engine='netcdf4')
ds = ds.squeeze()
# remove lev and ens coords
ds = ds.reset_coords(names=['lev', 'ens'], drop=True)

# select lats, lons, and dates within start_date, end_date and months
lat1, lat2, lon1, lon2 = bbox
ds = ds.sel(time=slice(start_date, end_date), lat=slice(lat1,lat2), lon=slice(lon1,lon2))
ds = select_months_ds(ds, start_mon, end_mon)

# add topo mask
mask = get_topo_mask(ds.lat, ds.lon) # create a elevation dataset with same grid spacing as ds
ds = ds.where(mask.bedrock >= elev_thres) # mask ds where elevation is less than 1000 m

# convert dataset to dataframe
df = ds.kidmap.to_dataframe(dim_order=['time', 'lat', 'lon'])
df = df.dropna(axis='rows')
# keep only rows that have trackID
trackID = df.groupby('time').kidmap.unique()
# trackID # this is all trackIDs that crossed the 1000 m threshold

id_df = trackID.to_frame() # converts to a pandas dataframe
id_df = id_df.reset_index() # reset the index
id_df = id_df.rename(columns={'time': 'date'}) # rename time column into date
id_df = id_df.set_index(pd.to_datetime(id_df['date'])) # reset the index as "date"
id_df.index = id_df.index.strftime("%Y-%m-%d") # make it so the index date is normalized to daily
id_df = id_df.rename(columns={'date': 'time'}) # rename the date column back to time
id_df = id_df.reset_index() # remove the index
id_df = id_df.explode('kidmap') # explode the dataframe based on trackID
# id_df

# load AR CAT (from Nash et al. 2021)
filepath = path_to_out + 'AR-types_ALLDAYS.csv'
ar_cat = pd.read_csv(filepath)
ar_cat = ar_cat.rename(columns={'Unnamed: 0': 'date'})
ar_cat = ar_cat.set_index(pd.to_datetime(ar_cat['date']))
ar_cat = select_months_df(ar_cat, start_mon, end_mon)
ar_cat.index = ar_cat.index.strftime("%Y-%m-%d")
ar_cat = ar_cat.drop(columns=['date'])
ar_cat = ar_cat.reset_index()
idx = ar_cat['AR_CAT'] > 0
ar_cat = ar_cat.loc[idx]

# merge id_df with ar_cat
merge_ar = pd.merge(id_df, ar_cat, how='outer', on='date')
track_ids = merge_ar.kidmap.unique() # get unique list of AR track IDs
ar_dates = merge_ar.time.unique() # get unique list of AR date/times (for later)
# merge_ar

# create df with trackID, ar_cat, start date, end date, and duration of AR (how long it is within HMA region)
ar = []
data = []
for i in [1, 2, 3]:
    idx = (merge_ar.AR_CAT == i)
    ar = merge_ar.loc[idx]

    for j, ids in enumerate(track_ids):
        idx = (ar.kidmap == ids)
        tmp = ar.loc[idx]
        start = pd.to_datetime(tmp.time.min())
        stop = pd.to_datetime(tmp.time.max()) + timedelta(hours=6)
        tmp = (stop - start)
        duration = tmp.total_seconds()/(3600) # convert to number of hours

        data.append([ids, i, start, stop, duration])
    
duration_df = pd.DataFrame(data, columns=['trackID', 'ar_cat', 'start_date', 'end_date', 'duration'])
duration_df = duration_df.dropna()
# duration_df

## open landslide df
def expand_grid(lat,lon):
    '''list all combinations of lats and lons using expand_grid(lat,lon)'''
    test = [(A,B) for A in lat for B in lon]
    test = np.array(test)
    test_lat = test[:,0]
    test_lon = test[:,1]
    full_grid = pd.DataFrame({'lat': test_lat, 'lon': test_lon})
    full_grid = full_grid.sort_values(by=['lat','lon'])
    full_grid = full_grid.reset_index(drop=True)
    return full_grid


fname = path_to_data + 'CH2_generated_data/Global_Landslide_Catalog_Export.csv' #TODO check this - is it the raw downloaded data?
landslide = pd.read_csv(fname)

# Select lat/lon grid
lonmin = 65
lonmax = 100
latmin = 20
latmax = 42

## Select Landslides within Southern Asia region
idx = (landslide.latitude >= latmin) & (landslide.latitude <= latmax) & (landslide.longitude >= lonmin) & (landslide.longitude <= lonmax)
landslide = landslide.loc[idx]
# set event time as index
landslide = landslide.set_index(pd.to_datetime(landslide.event_date))
# landslide.index = landslide.index.normalize()

# select only landslide dates that are between december and may
idx = (landslide.index.month >= 12) | (landslide.index.month <= 5)
landslide = landslide[idx]

# rename and reindex
landslide = landslide.rename(columns={"latitude": "lat", "longitude": "lon", "event_date": "event_time"})
landslide = landslide.reset_index()

# round event time to the nearest 6 hours
landslide['time'] = landslide['event_date'].dt.round('6H')
landslide = landslide.set_index(pd.to_datetime(landslide.time))

# select only landslide dates that are between december and may
idx = (landslide.index.month >= 12) | (landslide.index.month <= 5)
landslide = landslide[idx]

# now we want to see if there is an AR present at the same time and location as the landslides
# open the trackID for ARs
filename =  path_to_data + 'ar_catalog/globalARcatalog_ERA-Interim_1979-2019_v3.0.nc'
ar = xr.open_dataset(filename, engine='netcdf4')
ar = ar.squeeze()

# Select months
idx = (ar.time.dt.month >= 12) | (ar.time.dt.month <= 5)
kid = ar.kidmap.sel(time=idx) # trackID for indexing

# slice the dates so both ds match
kid = kid.sel(time=slice('1979-12-01 00', '2019-05-31 00:00'))


## for each landslide_id, if the lat/lon falls within an AR, keep that AR ID and landslide ID
landslideID = []
arID = []
landslide_lat = []
landslide_lon = []
for i, row in landslide.T.iteritems():
    t = kid.sel(lat=row['lat'], lon=row['lon'], time=row['time'], method='nearest').values
    # print(t)
    if t > 0:
        landslideID.append(row['event_id'])
        arID.append(t)
        landslide_lat.append(row['lat'])
        landslide_lon.append(row['lon'])
        
d = {'landslideID': landslideID, 'trackID': arID, 
     'landslide_lat': landslide_lat, 'landslide_lon': landslide_lon}
landslide_df = pd.DataFrame(data=d)
# convert the dtype for the trackID column
landslide_df = landslide_df.astype({'trackID': 'float64'})

# merge AR duration df and landslide DF
merged_data = pd.merge(duration_df, landslide_df, how='outer', on='trackID')
# merged_data 
# note the rows that do not have a date or time 
# are landslides that are associated with a specific AR that was not considered a "HMA AR"

# drop the rows that are not a HMA AR
idx = merged_data['ar_cat'] > 0
merged_data = merged_data.loc[idx]

def ar_ivt(df, ds, domains, clim_mean, clim_std):
    '''Calculate maximum IVT for a subregion in a ds and append to dataframe.
     For each range of AR event dates, we find the maximum IVT for the duration of the AR for every grid cell. 
    '''
    # the final IVT statistic to retain
    ivtdir_vals = []
    ivt_vals = []
    freeze_vals = []
    # loop through each AR track
    for i, (arcat, track) in enumerate(zip(df.ar_cat.values, df.trackID.values)):
        start = df.start_date.values[i]
        end = df.end_date.values[i]
        # print('Getting maximum between', start, end)
        print(start, end)
        # get bbox based on ar_cat
        bnds = domains[int(arcat)-1]
        # select only the time steps for AR event and specified domain
        tmp = ds.sel(time=slice(start, end), lat=slice(bnds[2], bnds[3]), lon=slice(bnds[0], bnds[1]))
        
        ## remove climatology from freezing level 
        clim_ar = clim_mean.sel(AR_CAT = arcat, lat=slice(bnds[2], bnds[3]), lon=slice(bnds[0], bnds[1]))
        std_ar = clim_std.sel(AR_CAT = arcat, lat=slice(bnds[2], bnds[3]), lon=slice(bnds[0], bnds[1]))
        tmp['z_new'] = (tmp.z - clim_ar.z)/std_ar.z # standardized anomalies
        ## average freezing level over time, lat, lon
        freeze = tmp['z_new'].max(['time'], skipna=True).mean(['lat', 'lon']).values
        freeze_vals.append(freeze.tolist())
        
        ### localized IVT maxima during event
        # event_max = tmp.where(tmp.ivt==tmp.ivt.max(), drop=True).squeeze()
        event_max = tmp.where(tmp.ivt==tmp.ivt.max(), drop=True).squeeze().load() # this was taking too long, decided to load earlier
        ## pull IVT and IVTDIR where ivt is max
        uvec = event_max.ivtu.values
        uvec = units.Quantity(uvec, "m/s")
        vvec = event_max.ivtv.values
        vvec = units.Quantity(vvec, "m/s")
        ivtdir = mpcalc.wind_direction(uvec, vvec)
        ivtdir_vals.append(ivtdir.item())
        ivt_vals.append(event_max.ivt.values.tolist())
        
        # # pull freezing level anomaly where ivt is max
        # freeze_vals.append(event_max.z_new.values)
        
    final = [ivtdir_vals, ivt_vals, freeze_vals]
        
    return final

def ar_precip(df, ds, domains, mode):
    '''Calculate precipitation statistics for a subregion in a ds and append to dataframe.
     Mode is chosen based on calculation. For each range of AR event dates, we calculate the total accumulated precip for every grid cell. 
     Then we remove all gridcells that had less than 1 mm of rain per event (these are not included in any calc)
     Then we weight the gridcells by the cosine of the latitude.
     Then based on mode selected, different statistics are retained:
         'mean-total' averages all viable gridcells within the subregion and retains this number
         'max-total' selects the maximum gridcell value to append
         'percentile-total' calcuates the 95th percentile and then averages all the grid cells that exceed this threshold
    '''
    # the final precip statistic to retain
    m1_vals = []

    for i, (arcat, track) in enumerate(zip(df.ar_cat.values, df.trackID.values)):
        start = df.start_date.values[i]
        end = df.end_date.values[i]
        # print('Getting maximum between', start, end)
        print(start, end)
        # get bbox based on ar_cat
        bnds = domains[int(arcat)-1]
        # select only the time steps for AR event and specified domain
        tmp = ds.sel(time=slice(start, end), lat=slice(bnds[2], bnds[3]), lon=slice(bnds[0], bnds[1]))

        ### event-total precipitation per event for every grid cell
        tmp = tmp.sum('time')
        ### mask out grid cells with less than 1 mm per event
        tmp2 = xr.where(cond=(tmp.prec > 1), x=tmp.prec, y=np.nan)

        ### area weighted
        # tmp = tmp2.weighted(tmp.weights)

        if mode == 'mean-total':
            ## mode 1: mean-total
            # average over gridcells in weighted subregion
            mean_tot = tmp.mean(['lat', 'lon'], skipna=True)
            # append to list
            m1_vals.append(mean_tot.values.tolist())
        elif mode == 'max-total':
            ## mode 2: max-total
            ### localized precip maxima during event
            event_max = tmp2.max(['lat', 'lon'])
            m1_vals.append(event_max.values.tolist())
        elif mode == 'percentile-total':
            ## mode 3: percentile-total
            ###  get 95th percentile thres
            q_thres = tmp2.quantile(0.95, dim=['lat', 'lon'], interpolation='linear')
            ## mask out grid cells below threshold
            perc_prec = xr.where(cond=(tmp2 > q_thres), x=tmp2, y=np.nan)
            # average over all grid cells skipping nans
            mean = perc_prec.mean(['lat', 'lon'], skipna=True)
            m1_vals.append(mean.values.tolist())

        
    return m1_vals

## pull wrflats and wrflons from first file
fname = path_to_data + 'wrf_hasia/d01/ivt/3hr/tmp_2015.nc'
tmp = xr.open_dataset(fname)
# print(tmp.time[:100])
# print(tmp.time[-100:])

## assign those lats to the other ds when you loop
wrflats = tmp.lat.values
wrflons = tmp.lon.values

fname = path_to_data + 'wrf_hasia/d02/prec/3hr/tmp_2014.nc'
tmp = xr.open_dataset(fname)
# print(tmp.time[:100])
# print(tmp.time[-100:])

## assign those lats to the other ds when you loop
wrflats2 = tmp.lat.values
wrflons2 = tmp.lon.values

def preprocess_ivt(ds):
    '''keep only the current year'''
    year = ds.time.dt.year.max().values
    ds = ds.assign_coords({"lon": wrflons, "lat": wrflats})
    if year == 1980:
        ds = ds
    else:
        ds = ds.sel(time=slice('{0}-01-01 00:00'.format(year), '{0}-12-31 21:00'.format(year)))
    return ds

def preprocess_prec(ds):
    '''keep only the current year'''
    year = ds.time.dt.year.max().values
    ds = ds.assign_coords({"lon": wrflons2, "lat": wrflats2})
    if year == 1980:
        ds = ds
    else:
        ds = ds.sel(time=slice('{0}-01-01 00:00'.format(year), '{0}-12-31 21:00'.format(year)))
    return ds

domains = ['d01', 'd02', 'd01']
varname_lst = ['ivt', 'prec', 'zerodegisotherm']

## loop through each ds
ds_lst = []
for i, (dom, varname) in enumerate(zip(domains, varname_lst)):
    print('Reading...', varname)
    if server == 'great':
        data_path = path_to_data + 'wrf_hasia/'
    else:
        data_path = path_to_data + 'wrf_preprocessed_data/wrf_6km/'
        
    filename_pattern = '{0}/{1}/3hr/tmp_*.nc'.format(dom, varname)
    fname = data_path + filename_pattern
    
    if varname == 'ivt':
        ds = xr.open_mfdataset(fname, preprocess=preprocess_ivt)
        ds = ds.assign(ivt=lambda ds: np.sqrt(ds.ivtu**2 + ds.ivtv**2))
    elif varname == 'prec':
        ds = xr.open_mfdataset(fname, preprocess=preprocess_prec)
        ## shift subtraction to get mm per hour 
        # # rain at next time step - rain at current time step
        ds = ds.shift(time=-1) - ds # if in xarray
    elif varname == 'geopotential':
        ds = ds.sel(lev=250.)
    elif varname == 'zerodegisotherm':
        ds = xr.open_mfdataset(fname, preprocess=preprocess_ivt)
    
    # subset to just ar days
    # ds = ds.sel(time = slice(start_date, end_date))
    # ds = select_months_ds(ds, start_mon, end_mon)
    ds = ds.sel(time = ar_dates[:-1])
    
    ds_lst.append(ds)
    
prec = ds_lst[1]

## combine IVT and freezing level into 1 ds
wrf_d01 = xr.combine_by_coords([ds_lst[0], ds_lst[2]]) 

# latmin, latmax, lonmin, lonmax
ext1 = [71, 79, 32, 37] # Western precip anomalies
ext2 = [69, 74, 37, 40] # Northwestern precip anomalies
ext3 = [90, 99, 24, 30] # Eastern precip anomalies

region_name = ['western', 'northwestern', 'eastern']
domains = [ext1, ext2, ext3]
clim_mean = xr.open_dataset('/home/sbarc/students/nash/data/wrf_hasia/d01/zerodegisotherm_ar_clim_new.nc') # freezing level climatology
clim_std = xr.open_dataset('/home/sbarc/students/nash/data/wrf_hasia/d01/zerodegisotherm_ar_std_new.nc') # freezing level standard deviation

print('Processing IVT and freezing level...')
## For each row, calculate the maximum IVT within the region between start and end
ivt_final = ar_ivt(merged_data, wrf_d01, domains, clim_mean, clim_std)

print('Processing precipitation...')
## For each row, calculate the maximum precip within the region between start and end
prec_final = ar_precip(merged_data, prec, domains, 'max-total')

## attach data to existing df
merged_data['ivt'] = ivt_final[1]
merged_data['ivtdir'] = ivt_final[0]
merged_data['freeze'] = ivt_final[2]
merged_data['prec'] = prec_final

# Export dataframes as csv
merged_data.to_csv(path_to_out + '{0}_ivt_ar_types_freezing_level_max_prec_20230317.csv'.format(ssn))