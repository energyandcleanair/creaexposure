library(testthat)

# ==============================================================================
# Tests for get_concentration() flow (concentration.R)
# ==============================================================================

# Helper: create a temporary .tif file and return its path
make_temp_tif <- function(vals = 10, nrows = 4, ncols = 4) {
  r <- terra::rast(nrows = nrows, ncols = ncols, vals = vals)
  path <- tempfile(fileext = ".tif")
  terra::writeRaster(r, path)
  path
}


# --- File not found -----------------------------------------------------------

test_that("get_concentration errors when file does not exist", {
  local_mocked_bindings(
    .concentration_path = function(...) "/nonexistent/path/fake.tif"
  )

  expect_error(
    get_concentration("no2", source = "larkin", year = 2011),
    "Concentration file not found"
  )
})

test_that("error message mentions process_* hint", {
  local_mocked_bindings(
    .concentration_path = function(...) "/nonexistent/path/fake.tif"
  )

  expect_error(
    get_concentration("no2", source = "larkin", year = 2011),
    "process_\\*\\(\\)"
  )
})


# --- Temporal scaling flow ----------------------------------------------------

test_that("scale_year == year skips scaling", {
  tif <- make_temp_tif(vals = 10)
  on.exit(unlink(tif))

  scaling_called <- FALSE
  local_mocked_bindings(
    .concentration_path = function(...) tif,
    .apply_temporal_scaling = function(...) {
      scaling_called <<- TRUE
      stop("should not be called")
    }
  )

  result <- get_concentration("no2", source = "larkin", year = 2011, scale_year = 2011)
  expect_false(scaling_called)
})

test_that("scale_year = NULL skips scaling", {
  tif <- make_temp_tif(vals = 10)
  on.exit(unlink(tif))

  scaling_called <- FALSE
  local_mocked_bindings(
    .concentration_path = function(...) tif,
    .apply_temporal_scaling = function(...) {
      scaling_called <<- TRUE
      stop("should not be called")
    }
  )

  result <- get_concentration("no2", source = "larkin", year = 2011, scale_year = NULL)
  expect_false(scaling_called)
})

test_that("scale_year != year triggers scaling", {
  tif <- make_temp_tif(vals = 10)
  on.exit(unlink(tif))

  scaling_called <- FALSE
  local_mocked_bindings(
    .concentration_path = function(...) tif,
    .apply_temporal_scaling = function(r, ...) {
      scaling_called <<- TRUE
      r
    }
  )

  result <- get_concentration("no2", source = "larkin", year = 2011, scale_year = 2023)
  expect_true(scaling_called)
})


# --- Unit assignment ----------------------------------------------------------

test_that("returned raster has correct unit", {
  tif <- make_temp_tif(vals = 10)
  on.exit(unlink(tif))

  local_mocked_bindings(
    .concentration_path = function(...) tif
  )

  r <- get_concentration("no2", source = "larkin", year = 2011)
  expect_equal(terra::units(r), enc2native("µg/m3"))
})

test_that("returned raster is named after the pollutant", {
  tif <- make_temp_tif(vals = 10)
  on.exit(unlink(tif))

  local_mocked_bindings(
    .concentration_path = function(...) tif
  )

  r <- get_concentration("no2", source = "larkin", year = 2011)
  expect_equal(names(r), "no2")
})


# --- Default resolution -------------------------------------------------------

test_that("source and version default correctly when omitted", {
  tif <- make_temp_tif(vals = 5)
  on.exit(unlink(tif))

  captured_args <- list()
  local_mocked_bindings(
    .concentration_path = function(pollutant, source, version, filename) {
      captured_args$source <<- source
      captured_args$version <<- version
      tif
    }
  )

  get_concentration("pm25", year = 2022)
  expect_equal(captured_args$source, "vandonkelaar")
  expect_equal(captured_args$version, "v5")
})


# --- Case insensitivity -------------------------------------------------------

test_that("pollutant is case-insensitive", {
  tif <- make_temp_tif(vals = 5)
  on.exit(unlink(tif))

  local_mocked_bindings(
    .concentration_path = function(...) tif
  )

  r <- get_concentration("PM25", source = "vandonkelaar", version = "v5", year = 2022)
  expect_equal(names(r), "pm25")
})
