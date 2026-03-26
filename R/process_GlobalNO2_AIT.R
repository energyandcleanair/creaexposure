#' Process GlobalNO2_AIT (https://zenodo.org/records/13842191/) concentration files
#'
#' Downloads the Annually.nc file from GlobalNO2_AIT, extracts the specified year,
#' converts from ppb to µg/m³, and writes to .tif.
#'
#' Data is organized as: no2/ait/default/no2_ait_{year}.tif
#'
#' @param year Integer. Year to process (2005–2023).
#' @param cache_dir Character. Directory to cache the downloaded .nc file.
#'   Defaults to "cache".
#'
#' @return Path to the processed .tif file (invisibly).
#' @export
process_GlobalNO2_AIT <- function(year) {

  out_dir <- .concentration_dir("no2", "ait", "default")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_filename <- glue::glue("no2_ait_{year}.tif")
  out_path <- file.path(out_dir, out_filename)

  if (file.exists(out_path)) {
    message(glue::glue("Already exists: {out_path}"))
    return(invisible(out_path))
  }

  # Raw data directory for this year
  raw_dir <- file.path(out_dir, "raw", year)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  nc_path <- file.path(raw_dir, "Annually.nc")

  if (!file.exists(nc_path)) {
    download_link <- "https://zenodo.org/records/13842191/files/Annually.nc?download=1"
    message("Downloading GlobalNO2_AIT Annually.nc...")
    httr::GET(download_link, httr::write_disk(nc_path, overwrite = TRUE))
  }

  nc_grid <- ncdf4::nc_open(nc_path)

  time_vals <- ncdf4::ncvar_get(nc_grid, "time")
  time_idx <- which(as.Date("2005-12-31") + time_vals == as.Date(paste0(year, "-12-31")))

  if (length(time_idx) == 0) {
    available <- format(as.Date("2005-12-31") + time_vals, "%Y")
    stop(glue::glue(
      "Year {year} not found in {nc_path}.\n",
      "Available years: {paste(available, collapse = ', ')}"
    ))
  }

  message(glue::glue("Extracting NO2_AiT for {year} (time index {time_idx})..."))
  r <- .nc_year_to_rast(nc_grid, time_idx)

  message(glue::glue("Saving to: {out_path}"))
  terra::writeRaster(r, out_path, overwrite = TRUE)

  # Clean up raw downloads
  ncdf4::nc_close(nc_grid)
  if (dir.exists(raw_dir)) {
    unlink(raw_dir, recursive = TRUE)
    message("Cleaned up raw downloads: ", raw_dir)
  }

  return(invisible(out_path))
}


#' Extract a single year slice from GlobalNO2_AIT nc and return a SpatRaster
#'
#' @param nc_grid An open ncdf4 connection.
#' @param time_idx Integer. Index along the time dimension to extract.
#' @return A SpatRaster in EPSG:4326 with values in µg/m³.
.nc_year_to_rast <- function(nc_grid, time_idx) {
  no2_slice <- ncdf4::ncvar_get(nc_grid, "NO2_AiT",
    start = c(1, 1, time_idx),
    count = c(-1, -1, 1)
  )

  lon <- ncdf4::ncvar_get(nc_grid, "longitude")
  lat <- ncdf4::ncvar_get(nc_grid, "latitude")

  r <- terra::rast(t(no2_slice))           # transpose [lon, lat] → [nlat, nlon]
  r <- terra::flip(r, direction = "vertical")
  terra::ext(r) <- terra::ext(min(lon), max(lon), min(lat), max(lat))
  terra::crs(r) <- "EPSG:4326"
  ###### Not sure about the unit. Maybe ppb? #####
  r <- r * 1.88 # ppb → µg/m³

  return(r)
}

