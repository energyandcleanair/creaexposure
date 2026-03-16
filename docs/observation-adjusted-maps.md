# Building observation-adjusted concentration maps

`creaexposure` offers a set of functions to build observation-adjusted maps by training models.

```r
maps <- build_maps(
  res     = "2pt5_min",
  regions = list("India" = "IN"),
  polls   = c("pm25", "no2"),
  year    = 2022,
  model   = "rf"
)
```


## Predictors

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
