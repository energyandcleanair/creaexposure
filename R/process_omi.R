#' Process OMI NO2 tropospheric column files
#'
#' Standardizes OMI NO2 files to .tif in no2/omi/default/
#' Files are expected to already be in .tif format with naming no2_omi_{year}.tif
#'
#' Raw data download:
#'   OMI NO2 tropospheric column from NASA Giovanni or similar
#'   Expected naming: no2_omi_{year}.tif
#'
#' @param year Integer. Year to process.
#'
#' @return Path to the processed .tif file (invisibly).
#' @export
process_omi <- function(year) {

  out_dir <- .concentration_dir("no2", "omi", "default")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_filename <- glue::glue("no2_omi_{year}.tif")
  out_path <- file.path(out_dir, out_filename)

  if (file.exists(out_path)) {
    message(glue::glue("Already exists: {out_path}"))
    return(invisible(out_path))
  }

  # Check legacy flat location
  legacy_path <- creahelpers::get_concentration_path(out_filename)
  if (file.exists(legacy_path)) {
    file.copy(legacy_path, out_path)
    message(glue::glue("Copied from legacy: {out_path}"))
    return(invisible(out_path))
  }

  stop(glue::glue(
    "OMI NO2 file not found: {out_filename}\n",
    "Expected at: {out_path}\n",
    "Or legacy: {legacy_path}"
  ))
}
