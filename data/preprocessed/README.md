## File structure

Once data is downloaded for both `/downloads` and `/dryad`, and preprocessing scripts are complete, the file structure of this directory should be as follows:

```
.
├── ivt
│   ├── daily_mean_clim_ivt.nc
│   ├── daily_std_clim_ivt.nc
│   ├── filtered_daily_mean_clim_ivt.nc
│   ├── out.wrf.d01.ivt.3hr_1979.nc
│   └── out.wrf.d01.ivt.3hr_*.nc
├── sr
│   ├── daily_mean_clim_sr.nc
│   ├── daily_std_clim_sr.nc
│   ├── filtered_daily_mean_clim_sr.nc
│   ├── out.wrf.d02.sr.daily_1979.nc
│   └── out.wrf.d02.sr.daily_*.nc
└── zerodegisotherm
    ├── ar_mean_clim_zerodegisotherm.nc
    ├── ar_std_clim_zerodegisotherm.nc
    ├── daily_mean_clim_zerodegisotherm.nc
    ├── daily_std_clim_zerodegisotherm.nc
    ├── filtered_daily_mean_clim_zerodegisotherm.nc
    ├── out.wrf.d01.zerodegisotherm.daily_1979.nc
    └── out.wrf.d01.zerodegisotherm.daily_*.nc
```
