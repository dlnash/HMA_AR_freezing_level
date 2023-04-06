## File structure

Once data is downloaded, the file structure of this directory should be as follows:

```
├── ERA5
│   ├── era5_hma_025deg_1hr_250z_198901.nc
│   ├── era5_hma_025deg_1hr_250z_201002.nc
│   ├── era5_hma_025deg_1hr_ivt_198901.nc
│   └── era5_hma_025deg_1hr_ivt_201002.nc
├── globalARcatalog_ERA-Interim_1979-2019_v3.0.nc
└── Global_Landslide_Catalog_Export.csv
```

## Directions for downloading data

1. Download ERA5 data using cds-api and ERA5_config.yml

```
# in download directory
conda activate cds
bash download_ERA5.sh
```

3. Download the AR detection result

Navigate to [https://dataverse.ucla.edu/dataverse/ar](https://dataverse.ucla.edu/dataverse/ar) and select 'globalARcatalog_ERA-Interim_1979-2019_v3.0.nc' to download. You will need to provide information regarding the use of this dataset.

2. Download Global Landslide Catalog 

Navigate to [https://data.nasa.gov/Earth-Science/Global-Landslide-Catalog-Export/dd9e-wu2v](https://data.nasa.gov/Earth-Science/Global-Landslide-Catalog-Export/dd9e-wu2v) and select Export, then CSV to download. 