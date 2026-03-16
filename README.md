# creaexposure

`creaexposure` is a CREA R package responsible for three things:

1. **Concentration retrieval** — unified API to access baseline concentration maps (PM2.5, NO2, O3) from multiple sources. It is used for instance in `creahia`.
2. **[Observation-adjusted concentration maps](docs/observation-adjusted-maps.md)** — fuse baseline priors with ground observations using GAM/RF models.
3. **[Standard-adjusted concentration maps](docs/standard-adjusted-maps.md)** — devise realistic concentration maps matching a regulatory target given historical trends.

```r
r <- get_concentration("pm25", source = "vandonkelaar", version = "v5", year = 2022)
```

## Setup

### Install creaexposure package
```r
remotes::install_github("energyandcleanair/creaexposure")
```

### Sync local folder with online bucket

The concentration maps are stored in CREA's GIS folder on Google Cloud Storage and typically synced in a local folder for the user to access them.

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

## Usage

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
- **`scale_year`** — applies temporal scaling (NO2 larkin: multiplicative OMI ratio). Only available for certain datasets.
- **`grid_raster`** — optional SpatRaster to resample the result to.

## Adding a new concentration source

See [docs/adding-a-concentration-source.md](docs/adding-a-concentration-source.md).
