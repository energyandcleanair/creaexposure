#' Process MERRA2 PM2.5 files
#'
#' Standardizes MERRA2 PM2.5 files to .tif in pm25/merra2/default/
#' Files are expected to already be in .tif format with naming pm25_merra2_{year}.tif
#'
#' Raw data download:
#'   MERRA2 PM2.5 from NASA Giovanni
#'   Expected naming: pm25_merra2_{year}.tif
#'
#' @param year Integer. Year to process.
#'
#' @return Path to the processed .tif file (invisibly).
#' @export
process_merra2 <- function(year) {

  out_dir <- .concentration_dir("pm25", "merra2", "default")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_filename <- glue::glue("pm25_merra2_{year}.tif")
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
    "MERRA2 PM2.5 file not found: {out_filename}\n",
    "Expected at: {out_path}\n",
    "Or legacy: {legacy_path}"
  ))
}
