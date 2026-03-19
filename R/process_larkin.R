#' Process Larkin Global LUR NO2 data
#'
#' Takes the raw Global_LUR_NO2_2011_16b_2.tif (full resolution),
#' aggregates 10x, reprojects to EPSG:4326, converts ppb to ug/m3,
#' and writes the processed file to no2/larkin/default/
#'
#' Raw file download:
#'   Manual download from Larkin et al. (2017) supplementary data
#'   File: Global_LUR_NO2_2011_16b_2.tif
#'   Place in: $GIS_DIR/concentration/no2/larkin/default/
#'
#' @return Path to the processed .tif file (invisibly).
#' @export
process_larkin <- function() {

  out_dir <- .concentration_dir("no2", "larkin", "default")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_filename <- "no2_agg10_ugm3_wsg84.tif"
  out_path <- file.path(out_dir, out_filename)

  if (file.exists(out_path)) {
    message(glue::glue("Already exists: {out_path}"))
    return(invisible(out_path))
  }

  raw_filename <- "Global_LUR_NO2_2011_16b_2.tif"
  raw_path <- file.path(out_dir, raw_filename)

  if (!file.exists(raw_path)) {
    # Check legacy flat location
    raw_path_legacy <- creahelpers::get_concentration_path(raw_filename)
    if (file.exists(raw_path_legacy)) {
      raw_path <- raw_path_legacy
    } else {
      stop(glue::glue(
        "Raw Larkin NO2 file not found.\n",
        "Expected at: {file.path(out_dir, raw_filename)}\n",
        "Download from Larkin et al. (2017) supplementary data."
      ))
    }
  }

  message("Processing Larkin NO2: aggregate 10x, reproject to WGS84, convert ppb -> ug/m3")
  no2 <- terra::rast(raw_path)
  terra::NAflag(no2) <- 128
  no2_agg10 <- terra::aggregate(no2, fact = 10, cores = parallel::detectCores() - 1)
  no2_wsg84 <- no2_agg10 %>% terra::project("epsg:4326")
  no2_wsg84 <- no2_wsg84 * 1.88  # ppb to ug/m3
  terra::writeRaster(no2_wsg84, out_path, overwrite = TRUE)
  message(glue::glue("Written: {out_path}"))

  return(invisible(out_path))
}
