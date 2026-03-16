# creaexposure

`creaexposure` is a CREA R package responsible for two things:

1. **Concentration retrieval** — unified API to access baseline concentration maps (PM2.5, NO2, O3) from multiple sources. It is used for instance in `creahia`.
2. **Concentration map updates** — builds either:
   - observation-adjusted concentration maps by existing priors with ground observations using GAM/RF models
   - standard-adjusted concentration maps by devising realistic concentration maps matching a regulatory target given historical conditions.

```r
# Get PM2.5 concentration map for 2022 (van Donkelaar v5)
r <- get_concentration("pm25", source = "vandonkelaar", version = "v5", year = 2022)

# Get NO2 (Larkin) scaled to 2023 via OMI ratio
r <- get_concentration("no2", source = "larkin", scale_year = 2023)

# Check unit
terra::units(r)  # "µg/m3"
```

## Setup
The concentration maps are stored in CREA's GIS folder on Google Cloud Storage and typically synced in a local folder for the user to access them.

### Sync local folder with online bucket
You first need to set up your local GIS folder root, using: 
```bash
# .Renviron
GIS_DIR=~/gis
```

To synchronise your whole gis folder with CREA's GCS, you should run:
```bash
source .Renviron
gsutil -m rsync -r $GIS_DIR gs://crea-data/gis
gsutil -m rsync -r gs://crea-data/gis $GIS_DIR
```

If you'd like to synchronise your concentration folder only, use:
```bash
source .Renviron
gsutil -m rsync -r $GIS_DIR/concentration gs://crea-data/gis/concentration
gsutil -m rsync -r gs://crea-data/gis/concentration $GIS_DIR/concentration
```

### Instal creaexposure package
```r
remotes::install_github("energyandcleanair/creaexposure")
```

## Concentration API

`creaexposure` offers a single entry point to retrieve any baseline concentration map:

```r
# PM2.5 (van Donkelaar v5, latest year)
r <- get_concentration("pm25")

# PM2.5 for a specific year/version
r <- get_concentration("pm25", source = "vandonkelaar", version = "v5", year = 2022)

# NO2 (Larkin, temporally scaled to 2023 via OMI ratio)
r <- get_concentration("no2", source = "larkin", scale_year = 2023)

# O3 (GEOSChem, specific layer)
r <- get_concentration("o3", source = "geoschem", variant = "sm8h")

# PM2.5 without sea salt/dust (ACAG = vandonkelaar v4)
r <- get_concentration("pm25", source = "vandonkelaar", version = "v4", variant = "no_ssdust")

# Available years
get_concentration_available_years("pm25", source = "vandonkelaar")
get_concentration_closest_year("pm25", source = "vandonkelaar", year = 2023)
```

### Sources

| Pollutant | Source | Versions | Unit | Notes |
|-----------|--------|----------|------|-------|
| PM2.5 | `vandonkelaar` (default) | `v5` (default), `v6`, `v4` | µg/m3 | v4 = ACAG, supports `no_ssdust` variant |
| PM2.5 | `merra2` | `default` | kg/m3 | Used as relative predictor (temporal diff) |
| PM2.5 | `tap` | `china` | µg/m3 | China 1km from tapdata.org.cn |
| NO2 | `larkin` (default) | `default` | µg/m3 | Single year (2011), use `scale_year` for temporal adjustment |
| NO2 | `omi` | `default` | molecules/cm2 | Used as ratio for temporal scaling |
| NO2 | `tap` | `china` | µg/m3 | China 1km from tapdata.org.cn |
| O3 | `geoschem` (default) | `default` | unknown | Variants: `m3m` (default), `sm8h` |

Units are stored in the registry and set on the returned `SpatRaster` via `terra::units(r)`.

### Parameters

- **`source`** — data provider. If NULL, uses pollutant default.
- **`version`** — dataset version (e.g. `v5`, `v6`, `v4`). If NULL, uses source default.
- **`variant`** — dataset variant for sub-products (e.g. `no_ssdust`, `m3m`/`sm8h`).
- **`scale_year`** — applies temporal scaling (NO2 larkin: multiplicative OMI ratio).
- **`grid_raster`** — optional SpatRaster to resample the result to.

## Folder structure

All concentration data lives under `$GIS_DIR/concentration/`:

```
$GIS_DIR/concentration/
├── pm25/
│   ├── vandonkelaar/
│   │   ├── v5/    V5GL0502.HybridPM25.Global.{YYYY}01-{YYYY}12.tif, SOURCE.md
│   │   ├── v6/    V6GL02.02.CNNPM25.Global.{YYYY}01-{YYYY}12.tif, SOURCE.md
│   │   └── v4/    ACAG_PM25_GWR_V4GL03_{YYYY}01_{YYYY}12_0p01.tif, SOURCE.md
│   ├── merra2/
│   │   └── default/   pm25_merra2_{YYYY}.tif, SOURCE.md
│   └── tap/
│       └── china/     {YYYY}.tif, SOURCE.md
├── no2/
│   ├── larkin/
│   │   └── default/   no2_agg10_ugm3_wsg84.tif, SOURCE.md
│   ├── omi/
│   │   └── default/   no2_omi_{YYYY}.tif, SOURCE.md
│   └── tap/
│       └── china/     {YYYY}.tif, SOURCE.md
├── o3/
│   └── geoschem/
│       └── default/   o3_m3m.tif, o3_sm8h.tif, SOURCE.md
├── [legacy flat files kept for backward compat]
└── tap/
    └── urls.txt       (TAP download links)
```

Files keep their **original filenames** for traceability. The folder hierarchy provides the structure.

Each version folder contains a `SOURCE.md` file documenting the data provenance (URL, coverage, resolution, unit, processing). When adding a new dataset, always create a `SOURCE.md` in its version folder.

### Processing raw data

Each source has a `process_*()` function that converts raw formats (.nc, .grd) to .tif:

```r
process_vandonkelaar(year = 2023, version = "v5")  # .nc → .tif
process_larkin()                                     # aggregate, reproject, ppb→ug/m3
process_geoschem()                                   # extract M3M/SM8h layers from .nc
process_omi(year = 2023)                             # copy from legacy location
process_merra2(year = 2023)                          # copy from legacy location
process_tap(year = 2023)                             # download, rasterize, gap-fill
```

## Exposure maps

Builds observation-adjusted maps by training models on the residual (observed − prior):

```r
maps <- build_maps(
  res     = "2pt5_min",
  regions = list("India" = "IN"),
  polls   = c("pm25", "no2"),
  year    = 2022,
  model   = "rf"
)
```

**Resolutions:** 30 arcsec (~1 km), 2.5 min (~5 km), 15 min (~27 km), 1 degree

### Predictors

| Name | Description | Source |
|------|-------------|--------|
| `pm25_prior` | PM2.5 baseline | van Donkelaar |
| `no2_prior` | NO2 baseline | Larkin (scaled via OMI) |
| `pm25_merra2_diff` | PM2.5 temporal change | MERRA2 |
| `no2_omi_diff` | NO2 temporal change | OMI |
| `pop`, `pop_05deg` | Population density | GPW v4 |
| `pop_ratio_log` | Population growth 2010→2020 | GPW v4 |
| `srtm`, `srtm_05deg`, `srtm_diff05deg` | Elevation / relief | SRTM |
| `distance_urban`, `grump` | Urban proximity | GRUMP v1 |
| `distance_coast` | Coastline distance | NASA |
| `gadm0`, `gadm1` | Admin boundaries | GADM |
| `lon`, `lat` | Coordinates | — |

## Adding a new concentration source

Checklist for adding a new dataset to the concentration API:

1. **Create `R/process_<source>.R`** — function that converts raw data to `.tif`
   - Output to `$GIS_DIR/concentration/{pollutant}/{source}/{version}/`
   - Use `.concentration_dir()` to build the output path
   - Skip if output already exists
   - Clean up raw downloads after processing

2. **Register in `.concentration_sources`** (in `R/concentration.R`)
   - Add entry under the pollutant with `default_version`, and for each version:
     - `unit` — the unit of the raster values (e.g. `"µg/m3"`)
     - `year_regex` — regex to discover years from filenames (or `NULL` if yearless)
     - `file_template` — function(year, variant) returning the `.tif` filename
     - `fixed_years` — if the source only covers specific years
   - If this should be the default source for a pollutant, update `.default_sources`

3. **Create `SOURCE.md`** in the version folder with:
   - URL / download instructions
   - Coverage (global, China, etc.)
   - Resolution
   - Unit
   - Processing function reference

4. **Update the sources table** in this README

5. **Test**
   ```r
   process_newsource(year = 2023)
   r <- get_concentration("pm25", source = "newsource")
   terra::units(r)  # should return the registered unit
   get_concentration_available_years("pm25", source = "newsource")
   ```

```
