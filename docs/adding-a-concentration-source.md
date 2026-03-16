# Adding a new concentration source

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

## Processing raw data

Each source has a `process_*()` function that converts raw formats (.nc, .grd) to .tif:

```r
process_vandonkelaar(year = 2023, version = "v5")  # .nc → .tif
process_larkin()                                     # aggregate, reproject, ppb→ug/m3
process_geoschem()                                   # extract M3M/SM8h layers from .nc
process_omi(year = 2023)                             # copy from legacy location
process_merra2(year = 2023)                          # copy from legacy location
process_tap(year = 2023)                             # download, rasterize, gap-fill
```

## Checklist

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

4. **Update the sources table** in the README

5. **Test**
   ```r
   process_newsource(year = 2023)
   r <- get_concentration("pm25", source = "newsource")
   terra::units(r)  # should return the registered unit
   get_concentration_available_years("pm25", source = "newsource")
   ```
