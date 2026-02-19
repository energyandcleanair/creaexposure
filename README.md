# creaexposure

R package for building observation-adjusted air pollution exposure maps (PM2.5 and NO2). Combines satellite-derived prior maps with ground-based observations using GAM or Random Forest models to produce locally-corrected exposure estimates at multiple resolutions.

## What it does

- Loads prior basemaps (van Donkelaar PM2.5, Larkin NO2) as spatial baselines
- Retrieves ground observations from CREADB via `rcrea`
- Trains GAM or Random Forest models to predict the spatial adjustment (observation minus prior) across unobserved locations
- Applies masking to limit extrapolation far from stations or urban areas
- Outputs adjusted raster maps with diagnostics and validation plots

**Resolutions:** 30 arcsec (~1 km), 2.5 min (~5 km), 15 min (~27 km), 1 degree

## Setup

Environment variable (`.Renviron`):
```bash
GIS_DIR=~/gis   # path to GIS data directory
```

Expected GIS directory structure:
```
$GIS_DIR/
├── concentration/      # Prior maps + MERRA2/OMI yearly rasters
├── population/         # Population density rasters
├── elevation/          # SRTM elevation
├── boundaries/         # Distance-to-coast
└── landcover/          # Urban distance, GRUMP
```

## Quick start

```r
maps <- build_maps(
  res       = "2pt5_min",
  regions   = list("India" = "IN"),
  polls     = c("pm25", "no2"),
  year      = 2022,
  model     = "rf"
)

maps$maps$pm25        # adjusted raster
maps$diagnostics      # validation plots
```

## Predictors

Below are currently available predictors:

| Name | Description | Source |
|------|-------------|--------|
| `pm25_prior` | PM2.5 satellite-derived baseline | van Donkelaar Global HybridPM2.5 |
| `no2_prior` | NO2 baseline | Larkin Global LUR NO2 (2011) |
| `pm25_merra2_diff` | PM2.5 temporal change: prior year → target year | MERRA2 (`gis/concentration/pm25_merra2_{year}.tif`) |
| `no2_omi_diff` | NO2 temporal change: 2011 → target year | OMI (`gis/concentration/no2_omi_{year}.tif`) |
| `pop` | Population density | GPW v4 2020 |
| `pop_05deg` | Population density smoothed at 0.5° | GPW v4 2020 |
| `pop_ratio_log` | Log population growth ratio 2010→2020 | GPW v4 |
| `srtm` | Elevation | SRTM 1km |
| `srtm_05deg` | Elevation smoothed at 0.5° | SRTM |
| `srtm_diff05deg` | Local elevation minus 0.5° background (topographic relief) | SRTM |
| `distance_urban` | Distance to nearest urban area | GRUMP v1 |
| `grump` | Urban/rural classification | GRUMP v1 |
| `distance_coast` | Distance to coastline | NASA Ocean Color |
| `gadm0` | Country-level admin boundary (rasterized) | GADM |
| `gadm1` | Sub-national admin boundary (rasterized) | GADM |
| `lon` / `lat` | Coordinates | — |


## Generate MERRA2 PM2.5 (`pm25_merra2_{year}`) and OMI NO2 (`no2_omi_{year}`) rasters
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