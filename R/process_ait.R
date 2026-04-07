#' Process GlobalNO2_AIT (https://zenodo.org/records/13842191) concentration files
#'
#' Downloads the Annually.nc file from GlobalNO2_AIT (which contains all available
#' years in a single file), extracts the requested year(s), converts from ppb to
#' µg/m³, and writes to .tif.
#'
#' Data is organized as: no2/ait/default/no2_ait_{year}.tif
#'
#' @param year Integer or integer vector. Year(s) to process (2005–2023).
#'   If NULL, all years available in the source file are processed.
#'
#' @return Character vector of paths to the processed .tif files (invisibly).
#' @export
process_ait <- function(year = NULL) {

  out_dir <- .concentration_dir("no2", "ait", "default")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Annually.nc contains all years — download once to a shared raw location
  raw_dir <- file.path(out_dir, "raw")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  nc_path <- file.path(raw_dir, "Annually.nc")

  if (!file.exists(nc_path)) {
    download_link <- "https://zenodo.org/records/13842191/files/Annually.nc?download=1"
    message("Downloading GlobalNO2_AIT Annually.nc (~100 MB, this may take a while)...")
    old_timeout <- getOption("timeout")
    options(timeout = 600)
    on.exit(options(timeout = old_timeout), add = TRUE)
    status <- utils::download.file(download_link, destfile = nc_path, mode = "wb", quiet = TRUE)
    if (!file.exists(nc_path) || status != 0 || file.size(nc_path) == 0) {
      stop("Failed to download GlobalNO2_AIT Annually.nc from ", download_link)
    }
  }

  nc_grid <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc_grid), add = TRUE)

  time_vals <- ncdf4::ncvar_get(nc_grid, "time")
  available_years <- as.integer(format(as.Date("2005-12-31") + time_vals, "%Y"))

  if (is.null(year)) {
    years_to_process <- available_years
  } else {
    years_to_process <- as.integer(year)
    missing_years <- setdiff(years_to_process, available_years)
    if (length(missing_years) > 0) {
      stop(glue::glue(
        "Year(s) {paste(missing_years, collapse = ', ')} not found in {nc_path}.\n",
        "Available years: {paste(available_years, collapse = ', ')}"
      ))
    }
  }

  out_paths <- vapply(years_to_process, function(yr) {
    out_filename <- glue::glue("no2_ait_{yr}.tif")
    out_path <- file.path(out_dir, out_filename)

    if (file.exists(out_path)) {
      message(glue::glue("Already exists: {out_path}"))
      return(out_path)
    }

    time_idx <- which(available_years == yr)
    message(glue::glue("Extracting NO2_AiT for {yr} (time index {time_idx})..."))
    r <- .ait_nc_year_to_rast(nc_grid, time_idx)

    message(glue::glue("Saving to: {out_path}"))
    terra::writeRaster(r, out_path, overwrite = TRUE)
    out_path
  }, character(1))

  # Clean up raw downloads
  unlink(raw_dir, recursive = TRUE)
  message("Cleaned up raw downloads: ", raw_dir)

  return(invisible(out_paths))
}


#' Extract a single year slice from GlobalNO2_AIT nc and return a SpatRaster
#'
#' @param nc_grid An open ncdf4 connection.
#' @param time_idx Integer. Index along the time dimension to extract.
#' @return A SpatRaster in EPSG:4326 with values in µg/m³.
.ait_nc_year_to_rast <- function(nc_grid, time_idx) {
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
  # Unit is ppbv (confirmed from Mu et al. 2026, https://doi.org/10.5194/essd-2025-821)
  r <- r * 1.88 # ppb → µg/m³ (standard NO2 conversion at 25°C, 1 atm)

  return(r)
}

