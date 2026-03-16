#' Process GEOSChem O3 concentration file
#'
#' Converts the raw .nc file to separate .tif layers (M3M and SM8h),
#' and places them in o3/geoschem/default/
#'
#' Raw file: O3_77e3b7-xmessy_mmd_kk.nc
#'   Place in: $GIS_DIR/concentration/o3/geoschem/default/
#'
#' Uses logic ported from creahia's rasterise_geoschem_nc().
#'
#' @return Paths to the processed .tif files (invisibly).
#' @export
process_geoschem <- function() {

  out_dir <- .concentration_dir("o3", "geoschem", "default")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_m3m <- file.path(out_dir, "o3_m3m.tif")
  out_sm8h <- file.path(out_dir, "o3_sm8h.tif")

  if (file.exists(out_m3m) && file.exists(out_sm8h)) {
    message("O3 geoschem .tif files already exist.")
    return(invisible(c(out_m3m, out_sm8h)))
  }

  raw_filename <- "O3_77e3b7-xmessy_mmd_kk.nc"
  raw_path <- file.path(out_dir, raw_filename)

  if (!file.exists(raw_path)) {
    # Check legacy flat location
    raw_path_legacy <- creahelpers::get_concentration_path(raw_filename)
    if (file.exists(raw_path_legacy)) {
      raw_path <- raw_path_legacy
    } else {
      stop(glue::glue(
        "Raw GEOSChem O3 file not found.\n",
        "Expected at: {file.path(out_dir, raw_filename)}\n",
        "Or legacy: {creahelpers::get_concentration_path(raw_filename)}"
      ))
    }
  }

  # Load the netCDF and extract layers
  r <- terra::rast(raw_path)
  layer_names <- names(r)

  # Find M3M and SM8h layers
  m3m_idx <- grep("M3M_lev31", layer_names)
  sm8h_idx <- grep("SM8h_lev31", layer_names)

  if (length(m3m_idx) == 0 || length(sm8h_idx) == 0) {
    stop(glue::glue(
      "Could not find expected layers in {raw_filename}.\n",
      "Available layers: {paste(layer_names, collapse=', ')}\n",
      "Expected layers containing 'M3M_lev31' and 'SM8h_lev31'."
    ))
  }

  message("Extracting O3 M3M layer...")
  r_m3m <- r[[m3m_idx[1]]]
  terra::crs(r_m3m) <- "epsg:4326"
  terra::writeRaster(r_m3m, out_m3m, overwrite = TRUE)
  message(glue::glue("Written: {out_m3m}"))

  message("Extracting O3 SM8h layer...")
  r_sm8h <- r[[sm8h_idx[1]]]
  terra::crs(r_sm8h) <- "epsg:4326"
  terra::writeRaster(r_sm8h, out_sm8h, overwrite = TRUE)
  message(glue::glue("Written: {out_sm8h}"))

  return(invisible(c(out_m3m, out_sm8h)))
}
