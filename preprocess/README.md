### Data and scripts used for preprocessing

1. WRF_clim/preprocess_wrf_climatology.py

Preprocesses chosen variables from 36 years of wrf output and save as yearly nc files (3-hourly resolution). Generates yearly files for the following variables: D01 Zero degree isotherm, D01 IVT, D02 Precipitation, and D02 Fraction of Frozen Precipitation. This preprocessed data is available for download on dryad. The data needed for this script includes: 

```
../home/hasia/
└── wrfout*
```

2. WRF_clim/calculate_filtered_climatology.py

Generates and saves daily variable files from 3-hourly, computes and saves daily standard deviation, and the filtered and unfiltered annual climatology for all variables processed in "preprocess_wrf_climatology.py". The data needed for this script includes:

```
../data/dryad/
├── ivt
│   ├── out.wrf.d01.ivt.3hr_1979.nc
│   └── out.wrf.d01.ivt.3hr_*.nc
├── prec
│   ├── out.wrf.d02.prec.3hr_1979.nc
│   └── out.wrf.d02.prec.3hr_*.nc
├── WRF_norris
│   ├── wrfout_d01_1988-12-31_03:00:00
│   ├── wrfout_d01_1989-01-05_03:00:00
│   ├── wrfout_d01_2010-02-04_03:00:00
│   ├── wrfout_d02_1988-12-31_03:00:00
│   ├── wrfout_d02_1989-01-05_03:00:00
│   └── wrfout_d02_2010-02-04_03:00:00
├── wvf
│   ├── wvflux_34.09N_74.02E.nc
│   └── wvflux_34.87N_72.66E.nc
└── zerodegisotherm
    ├── out.wrf.d01.zerodegisotherm.3hr_1979.nc
    └── out.wrf.d01.zerodegisotherm.3hr_*.nc
```

3. generate_df_ivt_prec_freezing.py

Generates data frame with IVT, precipitation, freezing level, and landslide information for each AR event. Saves to "../out/DJF_ivt_ar_types_freezing_level_max_prec.csv". The data needed for this script includes:

```
../data/dryad/
├── ivt
│   ├── out.wrf.d01.ivt.3hr_1979.nc
│   └── out.wrf.d01.ivt.3hr_*.nc
├── prec
│   ├── out.wrf.d02.prec.3hr_1979.nc
│   └── out.wrf.d02.prec.3hr_*.nc
├── WRF_norris
│   ├── wrfout_d01_2010-02-04_03:00:00
│   └── wrfout_d02_2010-02-04_03:00:00
└── zerodegisotherm
    ├── out.wrf.d01.zerodegisotherm.3hr_1979.nc
    └── out.wrf.d01.zerodegisotherm.3hr_*.nc
../data/downloaded
├── globalARcatalog_ERA-Interim_1979-2019_v3.0.nc
└── Global_Landslide_Catalog_Export.csv
```

4. WRF_clim/wvf_climatology.py

Calculates the vertical water vapor flux during AR days for the identified locations. Creates the files "../data/dryad/wvf/wvflux_*.nc", which are available for download on dryad. The data needed for this script includes: 

```
../home/hasia/
└── wrfout*

../data/downloaded
└── globalARcatalog_ERA-Interim_1979-2019_v3.0.nc
```