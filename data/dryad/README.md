# WRF High Mountain Asia Atmospheric River Freezing Level (1979-2015)
---

This study uses 36 years of Climate Forecast System Reanalysis (CFSR)[^1] dynamically downscaled over HMA to 20 km and 6.7 km spatial resolution and 3-hourly temporal resolution using the Advanced Weather Research and Forecasting (ARW-WRF; hereafter WRF) model version 3.7.1[^2][^3]. A subset of preprocessed data used in this analysis is provided here. This includes 35 years of 3-hourly 20 km integrated water vapor transport (IVT) and the height of the 0&deg; isotherm (hereafter, the freezing level), and 6.7 km precipitation, and fraction of frozen precipitation.


## Description of the data and file structure

ivt.zip contains annual netCDF (.nc) files (from 1979 to 2015) with WRF outer domain (20 km) 3-hourly integrated water vapor transport magnitude (ivt, kg m<sup>-1</sup> s<sup>-1</sup>), zonal IVT (ivtu, kg m<sup>-1</sup> s<sup>-1</sup>), meridional IVT (ivtv, kg m<sup>-1</sup> s<sup>-1</sup>), and integrated water vapor (iwv, mm). 

prec.zip contains annual .nc files (from 1979 to 2015) with WRF inner domain (6.7 km) 3-hourly accumulated precipitation (prec, mm) and fraction of frozen precipitation (sr, unitless). 

WRF_norris.zip contains the raw inner and outer domain WRF simulations for the two case study events. 

wvf.zip contains two .nc files (location indicated by the filename) with WRF inner domain (6.7 km) vertical water vapor flux (wvf, m s<sup>-1</sup>) for all time steps in which an Atmospheric River was identified at that location.

zerodegisotherm.zip contains annual .nc files (from 1979 to 2015) with WRF outer domain (20 km) 3-hourly freezing level (z, m). 

### File Structure

```
.
├── README.md
├── ivt.zip
│   ├── out.wrf.d01.ivt.3hr_1979.nc
│   └── out.wrf.d01.ivt.3hr_*.nc
├── prec.zip
│   ├── out.wrf.d02.prec.3hr_1979.nc
│   └── out.wrf.d02.prec.3hr_*.nc
├── WRF_norris.zip
│   ├── wrfout_d01_1988-12-31_03:00:00
│   ├── wrfout_d01_1989-01-05_03:00:00
│   ├── wrfout_d01_2010-02-04_03:00:00
│   ├── wrfout_d02_1988-12-31_03:00:00
│   ├── wrfout_d02_1989-01-05_03:00:00
│   └── wrfout_d02_2010-02-04_03:00:00
├── wvf.zip
│   ├── wvflux_34.09N_74.02E.nc
│   └── wvflux_34.87N_72.66E.nc
└── zerodegisotherm.zip
    ├── out.wrf.d01.zerodegisotherm.3hr_1979.nc
    └── out.wrf.d01.zerodegisotherm.3hr_*.nc
```


## Sharing/Access information

For detailed information on how WRF simulations were completed, see Norris et al.[^3][^4].

Data was derived from the following sources:
  * Climate Forecast System Reanalysis (CFSR)[^1] is available for download from [https://rda.ucar.edu/datasets/ds093.0/](https://rda.ucar.edu/datasets/ds093.0/)
  * Advanced Weather Research and Forecasting (ARW-WRF) version 3.7.1[^2] is available for download from [https://www2.mmm.ucar.edu/wrf/users/wrf_files/wrfv3.7/updates-3.7.1.html](https://www2.mmm.ucar.edu/wrf/users/wrf_files/wrfv3.7/updates-3.7.1.html)



## Code/Software

The code, including the conda environment, preprocessing scripts, figure scripts, and information on where to download other publically available data used for this analysis can be found at <https://github.com/dlnash/HMA_AR_freezing_level>.



## References

[^1]: Saha S, Moorthi S, Pan HL, Wu X, Wang JJ, Nadiga S, ... Goldberg M (2010) The NCEP Climate Forecast System Reanalysis. Bulletin of the American Meteorological Society 91(8):1015–1058, DOI <https://doi.org/10.1175/2010BAMS3001.1>

[^2]: Skamarock WC, Klemp JB, Dudhia J, Gill DO, Barker DM, Wang W, Powers JG (2008) A description of the advanced research WRF Version 3, NCAR technical note, Mesoscale and Microscale Meteorology Division. Tech. rep., National Center for Atmospheric Research, Boulder, CO, USA, DOI <http://dx.doi.org/10.5065/D68S4MVH>

[^4]: Norris J, Carvalho LMV, Jones C, Cannon F, Bookhagen B, Palazzi E, Tahir AA (2017) The spatiotemporal variability of precipitation over the Himalaya: evaluation of one-year WRF model simulation. Climate Dynamics 49(5-6):2179–2204, DOI <https://doi.org/10.1007/s00382-016-3414-y>

[^3]: Norris J, Carvalho LMV, Jones C, Cannon F (2019) Deciphering the contrasting climatic trends between the central Himalaya and Karakorum with 36 years of WRF simulations. Climate Dynamics 52(1-2):159–180, DOI <https://doi.org/10.1007/s00382-018-4133-3>