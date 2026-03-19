#' Process TAP (http://tapdata.org.cn/) concentration files
#'
#' Downloads annual .nc files from TAP, rasterizes the point-based data,
#' fills inland gaps, and writes to .tif.
#'
#' Data is organized as: pm25/tap/china/{year}.tif
#'
#' A shared urls.txt file is expected at $GIS_DIR/concentration/tap/urls.txt
#' with download links for all years.
#' Lines containing "Grid_lonlat" are grid coordinate files.
#' Lines containing "China_PM25_1km_{year}" are annual PM2.5 data files.
#'
#' @param year Integer. Year to process.
#' @param urls_file Character. Path to urls.txt. If NULL, uses
#'   $GIS_DIR/concentration/tap/urls.txt
#' @param pop Optional SpatRaster to use as template grid.
#'
#' @return Path to the processed .tif file (invisibly).
#' @export
process_tap <- function(year, urls_file = NULL, pop = NULL) {

  out_dir <- .concentration_dir("pm25", "tap", "china")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_filename <- paste0(year, ".tif")
  out_path <- file.path(out_dir, out_filename)

  if (file.exists(out_path)) {
    message(glue::glue("Already exists: {out_path}"))
    return(invisible(out_path))
  }

  # Raw data directory for this year
  raw_dir <- file.path(out_dir, "raw", year)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  # Locate urls.txt (shared file at pm25/tap/)
  if (is.null(urls_file)) {
    urls_file <- creahelpers::get_concentration_path(file.path("pm25", "tap", "urls.txt"))
  }
  if (!file.exists(urls_file)) {
    stop(glue::glue(
      "URLs file not found: {urls_file}\n",
      "Create a urls.txt with download links from http://tapdata.org.cn/"
    ))
  }

  urls <- readLines(urls_file, warn = FALSE)
  urls <- urls[urls != ""]

  grid_url <- urls[grep("Grid_lonlat", urls)]
  # Filter for this specific year
  year_pattern <- paste0("China_PM25_1km_", year)
  data_url <- urls[grep(year_pattern, urls)]

  if (length(grid_url) == 0) stop("Grid coordinates URL not found in urls.txt")
  if (length(data_url) == 0) stop(glue::glue("No PM2.5 data URL found for year {year} in urls.txt"))
  data_url <- data_url[1]

  # Download and extract grid coordinates
  grid_zip <- file.path(raw_dir, "Grid_lonlat.nc.zip")
  grid_nc <- file.path(raw_dir, "Grid_lonlat.nc")

  if (!file.exists(grid_nc)) {
    message("Downloading grid coordinates...")
    utils::download.file(grid_url, grid_zip, mode = "wb", quiet = FALSE)
    utils::unzip(grid_zip, exdir = raw_dir)
  }

  message("Reading grid coordinates...")
  nc_grid <- ncdf4::nc_open(grid_nc)
  lon <- ncdf4::ncvar_get(nc_grid, "Longitude")
  lat <- ncdf4::ncvar_get(nc_grid, "Latitude")
  ncdf4::nc_close(nc_grid)

  lon[lon == -999] <- NA
  lat[lat == -999] <- NA

  # Create template raster
  if (!is.null(pop)) {
    message("Using population grid as template")
    template <- pop[[1]]
    terra::values(template) <- NA
  } else {
    lon_range <- range(lon, na.rm = TRUE)
    lat_range <- range(lat, na.rm = TRUE)
    template <- terra::rast(
      xmin = lon_range[1], xmax = lon_range[2],
      ymin = lat_range[1], ymax = lat_range[2],
      resolution = 0.01, crs = "EPSG:4326"
    )
  }

  # Download and extract annual PM2.5 file
  filename <- paste0("China_PM25_1km_", year, ".nc.zip")
  zip_path <- file.path(raw_dir, filename)
  nc_path <- sub("\\.zip$", "", zip_path)

  if (!file.exists(nc_path)) {
    message("Downloading: ", filename)
    utils::download.file(data_url, zip_path, mode = "wb", quiet = FALSE)
    utils::unzip(zip_path, exdir = raw_dir)
  }

  message("Processing: ", basename(nc_path))
  nc <- ncdf4::nc_open(nc_path)
  pm25 <- ncdf4::ncvar_get(nc, "PM25")
  ncdf4::nc_close(nc)

  pm25[pm25 == -999] <- NA

  df <- data.frame(
    x = as.vector(lon),
    y = as.vector(lat),
    pm25 = as.vector(pm25)
  )
  df <- df[stats::complete.cases(df), ]

  pts <- terra::vect(df, geom = c("x", "y"), crs = "EPSG:4326")
  pm25_rast <- terra::rasterize(pts, template, field = "pm25", fun = mean)
  names(pm25_rast) <- paste0("PM25_", year)

  # Fill small inland gaps
  message("Filling inland gaps...")
  pm25_rast <- .fill_inland_gaps(pm25_rast)
  names(pm25_rast) <- paste0("PM25_", year)

  # Save
  message(glue::glue("Saving to: {out_path}"))
  terra::writeRaster(pm25_rast, out_path, overwrite = TRUE)

  # Clean up raw downloads
  if (dir.exists(raw_dir)) {
    unlink(raw_dir, recursive = TRUE)
    message("Cleaned up raw downloads: ", raw_dir)
  }

  return(invisible(out_path))
}


#' Fill small inland gaps using focal interpolation
#' Preserves ocean (large contiguous NA areas) by only filling gaps
#' that are surrounded by valid data
.fill_inland_gaps <- function(r, max_iterations = 10) {
  filled <- r
  total_filled <- 0

  for (i in seq_len(max_iterations)) {
    na_count_before <- sum(is.na(terra::values(filled)))

    filled <- terra::focal(
      filled, w = 3, fun = "mean",
      na.policy = "only", na.rm = TRUE
    )

    na_count_after <- sum(is.na(terra::values(filled)))
    cells_filled <- na_count_before - na_count_after
    total_filled <- total_filled + cells_filled

    if (cells_filled == 0) break
  }

  message("  Filled ", format(total_filled, big.mark = ","), " gap cells in ", i, " iterations")

  # Mask to original land extent to prevent ocean filling
  land_mask <- terra::focal(
    !is.na(r), w = 11, fun = "max", na.rm = TRUE
  )
  filled <- terra::mask(filled, land_mask, maskvalues = 0)

  return(filled)
}
