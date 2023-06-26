## Data paths for each notebook

### Downloads

Publically available data.

```
.
├── ERA5
│   ├── era5_hma_025deg_1hr_250z_198901.nc
│   ├── era5_hma_025deg_1hr_250z_201002.nc
│   ├── era5_hma_025deg_1hr_ivt_198901.nc
│   └── era5_hma_025deg_1hr_ivt_201002.nc
├── globalARcatalog_ERA-Interim_1979-2019_v3.0.nc
└── Global_Landslide_Catalog_Export.csv
```


#### Dryad

A subset of preprocessed WRF data available on Dryad for download.

```
.
├── ivt_split.zip
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

### Preprocessed 

Create these files after downloading all files from /data/downloads/ and /data/dryad/.

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

### Other (provided via github)

```
.
├── ar_casestudy.yml
├── AR-types_ALLDAYS.csv
├── namelist_hasia.input
└── namelist_hasia.wps
```


### Data used to generate figures

The following describes the data used for each of the notebooks in `../analysis/`.

1. fig01_norris_wrf_geogrid.ipynb

```
../data/dryad
└── WRF_norris
    ├── wrfout_d01_2010-02-04_03:00:00
    └── wrfout_d02_2010-02-04_03:00:00
../data/other
└── ar_casestudy.yml
```

2. fig02_trends_WRF_ivt.ipynb

```
../data/dryad
└── WRF_norris
    └── wrfout_d01_2010-02-04_03:00:00
../data/preprocessed
└── ivt
    ├── filtered_daily_mean_clim_zerodegisotherm.nc
    ├── out.wrf.d01.zerodegisotherm.daily_1979.nc
    └── out.wrf.d01.zerodegisotherm.daily_*.nc
../data/other
└── AR-types_ALLDAYS.csv
```

3. fig03_trends_WRF_0deg_isotherm.ipynb

```
../data/dryad
└── WRF_norris
    └── wrfout_d01_2010-02-04_03:00:00
../data/preprocessed
└── zerodegisotherm
    ├── filtered_daily_mean_clim_ivt.nc
    ├── out.wrf.d01.ivt.3hr_1979.nc
    └── out.wrf.d01.ivt.3hr_*.nc
../data/other
└── AR-types_ALLDAYS.csv
```

4. fig04_freezing_composites.ipynb

```
../data/dryad
└── WRF_norris
    └── wrfout_d01_2010-02-04_03:00:00
../data/preprocessed
└── prec
    ├── out.wrf.d02.prec.daily_1979.nc
    └── out.wrf.d02.prec.daily_*.nc
../out
└── DJF_ivt_ar_types_freezing_level_max_prec.csv
```

5. fig05_IVT_Precip_Scatter.ipynb

```
../out
└── DJF_ivt_ar_types_freezing_level_max_prec.csv
```

6. fig06_WRF-case_prec_summary.ipynb
```
../data/dryad
├── prec.zip
│   ├── out.wrf.d02.prec.3hr_1989.nc
│   └── out.wrf.d02.prec.3hr_2010.nc
├── WRF_norris.zip
│   └── wrfout_d01_2010-02-04_03:00:00
└── zerodegisotherm.zip
    ├── out.wrf.d01.zerodegisotherm.3hr_1989.nc
    └── out.wrf.d01.zerodegisotherm.3hr_2010.nc
../out
└── DJF_ivt_ar_types_freezing_level_max_prec.csv
```

7. fig07-8_synoptic_plots_all.ipynb

```
../data/dryad
├── WRF_norris.zip
│   └── wrfout_d01_2010-02-04_03:00:00
└── prec.zip
    ├── out.wrf.d02.prec.3hr_1989.nc
    └── out.wrf.d02.prec.3hr_2010.nc
../data/downloads
└── ERA5
    ├── era5_hma_025deg_1hr_250z_198901.nc
    ├── era5_hma_025deg_1hr_250z_201002.nc
    ├── era5_hma_025deg_1hr_ivt_198901.nc
    └── era5_hma_025deg_1hr_ivt_201002.nc

```

8. fig09-10_cross_section_static.ipynb

```
../data/dryad
├── ivt.zip
│   ├── out.wrf.d01.ivt.3hr_1989.nc
│   └── out.wrf.d01.ivt.3hr_2010.nc
├── prec.zip
│   ├── out.wrf.d02.prec.3hr_1989.nc
│   └── out.wrf.d02.prec.3hr_2010.nc
├── WRF_norris.zip
│   ├── wrfout_d01_1988-12-31_03:00:00
│   ├── wrfout_d01_1989-01-05_03:00:00
│   ├── wrfout_d01_2010-02-04_03:00:00
│   ├── wrfout_d02_1988-12-31_03:00:00
│   ├── wrfout_d02_1989-01-05_03:00:00
│   └── wrfout_d02_2010-02-04_03:00:00
└── zerodegisotherm.zip
    ├── out.wrf.d01.zerodegisotherm.3hr_1989.nc
    └── out.wrf.d01.zerodegisotherm.3hr_2010.nc
../data/other
└── ar_casestudy.yml
```

9. fig10_wvf_climatology.ipynb

```
../data/dryad
├── WRF_norris.zip
│   ├── wrfout_d02_1988-12-31_03:00:00
│   ├── wrfout_d02_1989-01-05_03:00:00
│   └── wrfout_d02_2010-02-04_03:00:00
└── wvf.zip
    ├── wvflux_34.09N_74.02E.nc
    └── wvflux_34.87N_72.66E.nc
```