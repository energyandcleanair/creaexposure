library(testthat)

# ==============================================================================
# Tests for temporal scaling (.apply_temporal_scaling in concentration.R)
# ==============================================================================

# Helper: create a raster with resolution small enough for terra::focalMat(d=1)
make_raster <- function(vals, nrows = 20, ncols = 20,
                        xmin = 0, xmax = 10, ymin = 0, ymax = 10) {
  terra::rast(
    nrows = nrows, ncols = ncols,
    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
    vals = vals
  )
}


test_that("scaling multiplies NO2 larkin by OMI ratio", {
  base_rast  <- make_raster(10)
  omi_base   <- make_raster(100)
  omi_target <- make_raster(200)

  local_mocked_bindings(
    get_concentration = function(pollutant, source, year, ...) {
      if (source == "omi" && year == 2011) return(omi_base)
      if (source == "omi" && year == 2023) return(omi_target)
      stop("unexpected call")
    }
  )

  result <- creaexposure:::.apply_temporal_scaling(
    base_rast, "no2", "larkin", base_year = 2011, target_year = 2023, grid_raster = NULL
  )

  # Uniform OMI ratio of 200/100 = 2, so values should be ~20
  vals <- terra::values(result)
  expect_true(all(!is.na(vals)))
  expect_equal(mean(vals), 20, tolerance = 1)
})

test_that("scaling with uniform OMI ratio of 1 is identity", {
  base_rast <- make_raster(42)
  omi_same  <- make_raster(100)

  local_mocked_bindings(
    get_concentration = function(pollutant, source, year, ...) omi_same
  )

  result <- creaexposure:::.apply_temporal_scaling(
    base_rast, "no2", "larkin", base_year = 2011, target_year = 2023, grid_raster = NULL
  )

  vals <- terra::values(result)
  expect_equal(mean(vals), 42, tolerance = 0.5)
})

test_that("scaling is a no-op for non-larkin sources", {
  r <- make_raster(42)

  result <- creaexposure:::.apply_temporal_scaling(
    r, "no2", "omi", base_year = 2011, target_year = 2023, grid_raster = NULL
  )

  expect_equal(terra::values(result), terra::values(r))
})

test_that("scaling is a no-op for PM2.5", {
  r <- make_raster(15)

  result <- creaexposure:::.apply_temporal_scaling(
    r, "pm25", "vandonkelaar", base_year = 2020, target_year = 2023, grid_raster = NULL
  )

  expect_equal(terra::values(result), terra::values(r))
})

test_that("scaling preserves spatial properties", {
  base_rast  <- make_raster(5)
  omi_base   <- make_raster(50)
  omi_target <- make_raster(75)

  local_mocked_bindings(
    get_concentration = function(pollutant, source, year, ...) {
      if (source == "omi" && year == 2011) return(omi_base)
      if (source == "omi" && year == 2023) return(omi_target)
      stop("unexpected call")
    }
  )

  result <- creaexposure:::.apply_temporal_scaling(
    base_rast, "no2", "larkin", base_year = 2011, target_year = 2023, grid_raster = NULL
  )

  expect_equal(terra::nrow(result), terra::nrow(base_rast))
  expect_equal(terra::ncol(result), terra::ncol(base_rast))
  expect_equal(terra::ext(result), terra::ext(base_rast))
})

test_that("scaling halves values when OMI target is half of base", {
  base_rast  <- make_raster(100)
  omi_base   <- make_raster(200)
  omi_target <- make_raster(100)

  local_mocked_bindings(
    get_concentration = function(pollutant, source, year, ...) {
      if (source == "omi" && year == 2015) return(omi_base)
      if (source == "omi" && year == 2023) return(omi_target)
      stop("unexpected call")
    }
  )

  result <- creaexposure:::.apply_temporal_scaling(
    base_rast, "no2", "larkin", base_year = 2015, target_year = 2023, grid_raster = NULL
  )

  vals <- terra::values(result)
  expect_equal(mean(vals), 50, tolerance = 1)
})
