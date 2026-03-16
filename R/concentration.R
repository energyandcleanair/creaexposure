# ==============================================================================
# Unified concentration retrieval API
#
# Reads processed .tif files from the standardized GIS folder structure:
#   $GIS_DIR/concentration/{pollutant}/{source}/{version}/
#
# Files keep their original filenames for traceability. Each source/version
# defines a year_regex and file_template to discover and construct filenames.
#
# Processing (raw -> .tif) is handled by process_*.R files.
# ==============================================================================


# --- Source Registry ----------------------------------------------------------
# Each source defines:
#   default_version: version to use when none specified
#   versions: list of version configs, each with:
#     unit: unit of the raster values (e.g. "ug/m3", "kg/m3", "molecules/cm2")
#     year_regex: regex to extract year (capture group 1) from .tif filenames
#     file_template: function(year, variant) -> filename
#     fixed_years: (optional) hardcoded years when discovery isn't possible
#     default_variant: (optional) default variant name
#     variants: (optional) list of variant names with their file templates

.concentration_sources <- list(
  pm25 = list(
    vandonkelaar = list(
      default_version = "v5",
      versions = list(
        v5 = list(
          unit = "µg/m3",
          year_regex = "^V5GL0502\\.HybridPM25\\.Global\\.(\\d{4})\\d{2}-\\d{4}\\d{2}\\.tif$",
          file_template = function(year, variant = NULL) {
            glue::glue("V5GL0502.HybridPM25.Global.{year}01-{year}12.tif")
          }
        ),
        v6 = list(
          unit = "µg/m3",
          year_regex = "^V6GL02\\.02\\.CNNPM25\\.Global\\.(\\d{4})\\d{2}-\\d{4}\\d{2}\\.tif$",
          file_template = function(year, variant = NULL) {
            glue::glue("V6GL02.02.CNNPM25.Global.{year}01-{year}12.tif")
          }
        ),
        v4 = list(
          # ACAG = vandonkelaar v4
          unit = "µg/m3",
          year_regex = "^ACAG_PM25_GWR_V4GL03_(\\d{4})\\d{2}_\\d{4}\\d{2}_0p01\\.tif$",
          default_variant = NULL,
          file_template = function(year, variant = NULL) {
            if (!is.null(variant) && variant == "no_ssdust") {
              glue::glue("ACAG_PM25_noDUSTnoSEASALT_GWR_V4GL03_{year}01_{year}12_0p01.tif")
            } else {
              glue::glue("ACAG_PM25_GWR_V4GL03_{year}01_{year}12_0p01.tif")
            }
          }
        )
      )
    ),
    tap = list(
      default_version = "china",
      versions = list(
        china = list(
          unit = "µg/m3",
          year_regex = "^(\\d{4})\\.tif$",
          file_template = function(year, variant = NULL) paste0(year, ".tif")
        )
      )
    ),
    merra2 = list(
      default_version = "default",
      versions = list(
        default = list(
          unit = "kg/m3",
          year_regex = "^pm25_merra2_(\\d{4})\\.tif$",
          file_template = function(year, variant = NULL) glue::glue("pm25_merra2_{year}.tif")
        )
      )
    )
  ),
  no2 = list(
    larkin = list(
      default_version = "default",
      versions = list(
        default = list(
          unit = "µg/m3",
          year_regex = NULL,
          file_template = function(year, variant = NULL) "no2_agg10_ugm3_wsg84.tif",
          fixed_years = c(2011L)
        )
      )
    ),
    omi = list(
      default_version = "default",
      versions = list(
        default = list(
          unit = "molecules/cm2",
          year_regex = "^no2_omi_(\\d{4})\\.tif$",
          file_template = function(year, variant = NULL) glue::glue("no2_omi_{year}.tif")
        )
      )
    ),
    tap = list(
      default_version = "china",
      versions = list(
        china = list(
          unit = "µg/m3",
          year_regex = "^(\\d{4})\\.tif$",
          file_template = function(year, variant = NULL) paste0(year, ".tif")
        )
      )
    )
  ),
  o3 = list(
    geoschem = list(
      default_version = "default",
      versions = list(
        default = list(
          unit = "unknown",
          year_regex = NULL,
          fixed_years = NULL,
          default_variant = "m3m",
          file_template = function(year = NULL, variant = "m3m") {
            switch(variant,
              "m3m" = "o3_m3m.tif",
              "sm8h" = "o3_sm8h.tif",
              stop("Unknown O3 variant: ", variant, ". Expected 'm3m' or 'sm8h'.")
            )
          }
        )
      )
    )
  )
)

.default_sources <- list(
  pm25 = "vandonkelaar",
  no2 = "larkin",
  o3 = "geoschem"
)


# --- Public API ---------------------------------------------------------------

#' Get concentration map for a given pollutant
#'
#' Retrieves a processed .tif concentration raster from the standardized
#' folder structure.
#'
#' @param pollutant Character. One of "pm25", "no2", "o3".
#' @param source Character. Source name (e.g., "vandonkelaar", "larkin", "tap").
#'   If NULL, uses the default source for the pollutant.
#' @param year Numeric or character. Year of the concentration map.
#'   If NULL, uses the latest available year. Supports "midYYYYmidYYYY" format.
#' @param version Character. Version of the source dataset (e.g., "v5", "v6").
#'   If NULL, uses the default version for the source.
#' @param variant Character. Variant of the dataset (e.g., "no_ssdust" for
#'   vandonkelaar v4, "m3m"/"sm8h" for geoschem O3). If NULL, uses default.
#' @param grid_raster Optional SpatRaster to resample the result to.
#' @param scale_year Optional numeric. If provided, applies temporal scaling
#'   (e.g., for NO2 larkin: applies OMI ratio from base year to scale_year).
#'
#' @return A terra::SpatRaster
#' @export
get_concentration <- function(pollutant,
                              source = NULL,
                              year = NULL,
                              version = NULL,
                              variant = NULL,
                              grid_raster = NULL,
                              scale_year = NULL) {

  pollutant <- tolower(pollutant)
  if (is.null(source)) source <- get_concentration_default_source(pollutant)
  source_config <- .get_source_config(pollutant, source)
  if (is.null(version)) version <- source_config$default_version
  version_config <- .get_version_config(pollutant, source, version)
  if (is.null(variant)) variant <- version_config$default_variant

  # Resolve year (skip for non-year-based sources like geoschem O3)
  has_years <- !is.null(version_config$year_regex) || !is.null(version_config$fixed_years)
  if (has_years) {
    year <- .parse_year(year)
    if (is.null(year)) {
      available <- get_concentration_available_years(pollutant, source, version)
      if (length(available) == 0) {
        stop(glue::glue("No years available for {pollutant}/{source}/{version}"))
      }
      year <- max(available)
    } else {
      year <- get_concentration_closest_year(pollutant, source, year, version)
    }
  }

  # Build filename and load
  filename <- version_config$file_template(year, variant)
  path <- .concentration_path(pollutant, source, version, filename)
  if (!file.exists(path)) {
    stop(glue::glue(
      "Concentration file not found: {path}\n",
      "Run the appropriate process_*() function to generate it."
    ))
  }
  r <- terra::rast(path)

  # Apply temporal scaling if requested
  if (!is.null(scale_year) && !is.null(year) && scale_year != year) {
    r <- .apply_temporal_scaling(r, pollutant, source, year, scale_year, grid_raster)
  }

  # Resample to grid if provided
  if (!is.null(grid_raster)) {
    r <- r %>% terra::resample(terra::rast(grid_raster), method = "bilinear")
  }

  names(r) <- pollutant
  if (!is.null(version_config$unit)) {
    terra::units(r) <- version_config$unit
  }
  return(r)
}


#' List available years for a pollutant/source combination
#'
#' Scans the folder for .tif files matching the source's year_regex.
#'
#' @param pollutant Character.
#' @param source Character. If NULL, uses default source.
#' @param version Character. If NULL, uses default version.
#'
#' @return Integer vector of available years, sorted.
#' @export
get_concentration_available_years <- function(pollutant,
                                              source = NULL,
                                              version = NULL) {

  pollutant <- tolower(pollutant)
  if (is.null(source)) source <- get_concentration_default_source(pollutant)
  source_config <- .get_source_config(pollutant, source)
  if (is.null(version)) version <- source_config$default_version
  version_config <- .get_version_config(pollutant, source, version)

  # Some sources have fixed years (e.g., larkin = 2011 only)
  if (!is.null(version_config$fixed_years)) {
    return(version_config$fixed_years)
  }

  year_regex <- version_config$year_regex
  if (is.null(year_regex)) return(integer(0))

  dir_path <- .concentration_dir(pollutant, source, version)
  if (!dir.exists(dir_path)) return(integer(0))

  files <- list.files(dir_path, pattern = "\\.tif$", full.names = FALSE)
  m <- stringr::str_match(files, year_regex)
  years <- as.integer(m[, 2])
  sort(years[!is.na(years)])
}


#' Find the closest available year for a pollutant/source
#'
#' @param pollutant Character.
#' @param source Character. If NULL, uses default source.
#' @param year Numeric or character. Target year. Supports "midYYYYmidYYYY" format.
#' @param version Character. If NULL, uses default version.
#'
#' @return Integer. The closest available year.
#' @export
get_concentration_closest_year <- function(pollutant,
                                           source = NULL,
                                           year,
                                           version = NULL) {

  year <- .parse_year(year)
  available <- get_concentration_available_years(pollutant, source, version)

  if (length(available) == 0) {
    stop(glue::glue("No years available for {pollutant}/{source}"))
  }

  diffs <- abs(available - year)
  available[which.min(diffs)]
}


#' Get the default source for a pollutant
#'
#' @param pollutant Character.
#' @return Character. Default source name.
#' @export
get_concentration_default_source <- function(pollutant) {
  pollutant <- tolower(pollutant)
  src <- .default_sources[[pollutant]]
  if (is.null(src)) {
    stop(glue::glue("No default source defined for pollutant: {pollutant}"))
  }
  src
}


# --- Internal helpers ---------------------------------------------------------

.concentration_dir <- function(pollutant, source, version) {
  creahelpers::get_concentration_path(file.path(pollutant, source, version))
}

.concentration_path <- function(pollutant, source, version, filename) {
  creahelpers::get_concentration_path(file.path(pollutant, source, version, filename))
}

.get_source_config <- function(pollutant, source) {
  poll_sources <- .concentration_sources[[pollutant]]
  if (is.null(poll_sources)) {
    stop(glue::glue("Unknown pollutant: {pollutant}"))
  }
  config <- poll_sources[[source]]
  if (is.null(config)) {
    stop(glue::glue("Unknown source '{source}' for pollutant '{pollutant}'. ",
                     "Available: {paste(names(poll_sources), collapse=', ')}"))
  }
  config
}

.get_version_config <- function(pollutant, source, version) {
  source_config <- .get_source_config(pollutant, source)
  vc <- source_config$versions[[version]]
  if (is.null(vc)) {
    stop(glue::glue("Unknown version '{version}' for {pollutant}/{source}. ",
                     "Available: {paste(names(source_config$versions), collapse=', ')}"))
  }
  vc
}

.parse_year <- function(year) {
  if (is.null(year)) return(NULL)
  year <- as.character(year)
  if (stringr::str_detect(year, "^mid\\d{4}mid\\d{4}$")) {
    as.numeric(stringr::str_match(year, "mid(\\d{4})mid")[2])
  } else {
    as.numeric(year)
  }
}

.apply_temporal_scaling <- function(r, pollutant, source, base_year, target_year, grid_raster) {
  if (pollutant == "no2" && source == "larkin") {
    # Multiplicative OMI ratio for temporal adjustment
    omi_base <- get_concentration("no2", source = "omi", year = base_year)
    omi_target <- get_concentration("no2", source = "omi", year = target_year)

    if (!is.null(grid_raster)) {
      omi_base <- omi_base %>% terra::resample(terra::rast(grid_raster), method = "bilinear")
      omi_target <- omi_target %>% terra::resample(terra::rast(grid_raster), method = "bilinear")
    }

    # Smooth to reduce noise
    focal_w <- terra::focalMat(omi_base, d = 1, type = "circle")
    focal_w[focal_w > 0] <- 1
    omi_base_smooth <- terra::focal(omi_base, w = focal_w, fun = mean, na.rm = TRUE, pad = TRUE)
    omi_target_smooth <- terra::focal(omi_target, w = focal_w, fun = mean, na.rm = TRUE, pad = TRUE)

    ratio <- omi_target_smooth / omi_base_smooth
    r <- r * ratio
  }
  return(r)
}
