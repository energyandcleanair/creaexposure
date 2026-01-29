# CREA EXPOSURE
Creating adjusted global & country maps based on observations and literature basemaps




## Generate MERRA2 PM2.5 and OMI NO2 rasters
Predictors for PM2.5 and NO2 exposure maps include MERRA2 PM2.5 and OMI NO2 respectively. Namely, we consider the difference between:
- the year we're aiming to build an exposure map for
- the latest/closest year available from PM2.5 and NO2 priors (i.e. van Donkelaar and Larkin respectively)

To build these yearly maps, if not already available in the GIS folder (`concentration` subfolder), visit https://giovanni.gsfc.nasa.gov/giovanni/

Here are the datasets you need to build yearly-averaged maps of:
- PM2.5: Total Surface Mass Concentration - PM 2.5 M2TMNXAER v5.12.4
- NO2: NO2 Tropospheric Column (30% Cloud Screened) (OMNO2d v003)

Once built, download in `.tif` format and copy the files in `gis/concentration` with the names:
- PM2.5: pm25_merra2_{year}.tif
- NO2: no2_omi_{year}.tif

These files will be used by `data.pm25_merra2_diff` and `data.no2_omi_diff` respectively.