#' Process van Donkelaar PM2.5 concentration files
#'
#' Converts raw .nc files to .tif and places them in the standardized
#' folder structure: pm25/vandonkelaar/{version}/{filename}.tif
#'
#' Raw .nc files can be downloaded from:
#' - v5: https://sites.wustl.edu/acag/datasets/surface-pm2-5/
#'   Files: V5GL0502.HybridPM25.Global.{YYYY}01-{YYYY}12.nc
#' - v6: https://sites.wustl.edu/acag/datasets/surface-pm2-5/
#'   Files: V6GL02.02.CNNPM25.Global.{YYYY}01-{YYYY}12.nc
#' - v4 (ACAG): already .tif, just needs placement in folder structure
#'
#' @param year Integer. Year to process.
#' @param version Character. "v5" (default), "v6", or "v4".
#'
#' @return Path to the processed .tif file (invisibly).
#' @export
process_vandonkelaar <- function(year, version = "v5") {

  version_config <- .get_version_config("pm25", "vandonkelaar", version)
  out_dir <- .concentration_dir("pm25", "vandonkelaar", version)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_filename <- version_config$file_template(year)
  out_path <- file.path(out_dir, out_filename)

  if (file.exists(out_path)) {
    message(glue::glue("Already exists: {out_path}"))
    return(invisible(out_path))
  }

  if (version %in% c("v5", "v6")) {
    # Convert .nc to .tif
    basename <- switch(version,
      "v5" = "V5GL0502.HybridPM25.Global",
      "v6" = "V6GL02.02.CNNPM25.Global"
    )
    nc_filename <- glue::glue("{basename}.{year}01-{year}12.nc")
    nc_path <- file.path(out_dir, nc_filename)

    if (!file.exists(nc_path)) {
      # Also check legacy flat location
      nc_path_legacy <- creahelpers::get_concentration_path(nc_filename)
      if (file.exists(nc_path_legacy)) {
        nc_path <- nc_path_legacy
      } else {
        stop(glue::glue(
          "Raw .nc file not found.\n",
          "Expected at: {file.path(out_dir, nc_filename)}\n",
          "Download from https://sites.wustl.edu/acag/datasets/surface-pm2-5/"
        ))
      }
    }

    message(glue::glue("Converting {nc_filename} -> {out_filename}"))
    r <- raster::raster(nc_path)
    raster::writeRaster(r, out_path, overwrite = TRUE)
    message(glue::glue("Written: {out_path}"))

  } else if (version == "v4") {
    # ACAG v4 files are already .tif, check if they exist in the right place
    # They may need to be moved from the legacy flat folder
    legacy_path <- creahelpers::get_concentration_path(out_filename)
    if (file.exists(legacy_path) && !file.exists(out_path)) {
      file.copy(legacy_path, out_path)
      message(glue::glue("Copied from legacy: {out_path}"))
    } else if (!file.exists(out_path)) {
      stop(glue::glue("ACAG v4 file not found: {out_filename}"))
    }
  }

  return(invisible(out_path))
}
