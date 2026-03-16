library(testthat)

# ==============================================================================
# Tests for the concentration API (concentration.R)
# ==============================================================================

# --- .parse_year --------------------------------------------------------------

test_that(".parse_year handles standard years", {
  expect_equal(creaexposure:::.parse_year(2023), 2023)
  expect_equal(creaexposure:::.parse_year("2023"), 2023)
})

test_that(".parse_year handles mid-year format", {
  expect_equal(creaexposure:::.parse_year("mid2023mid2024"), 2023)
  expect_equal(creaexposure:::.parse_year("mid2019mid2020"), 2019)
})

test_that(".parse_year handles NULL", {
  expect_null(creaexposure:::.parse_year(NULL))
})


# --- Source registry ----------------------------------------------------------

test_that("default sources are defined for all pollutants", {
  expect_equal(get_concentration_default_source("pm25"), "vandonkelaar")
  expect_equal(get_concentration_default_source("no2"), "larkin")
  expect_equal(get_concentration_default_source("o3"), "geoschem")
})

test_that("unknown pollutant errors", {
  expect_error(get_concentration_default_source("co2"), "No default source")
})

test_that("unknown source errors", {
  expect_error(
    creaexposure:::.get_source_config("pm25", "nonexistent"),
    "Unknown source"
  )
})

test_that("unknown version errors", {
  expect_error(
    creaexposure:::.get_version_config("pm25", "vandonkelaar", "v99"),
    "Unknown version"
  )
})


# --- Path resolution ----------------------------------------------------------

test_that(".concentration_dir builds correct path", {
  dir <- creaexposure:::.concentration_dir("pm25", "vandonkelaar", "v5")
  expect_true(grepl("pm25/vandonkelaar/v5$", dir))
})

test_that(".concentration_path builds correct path", {
  path <- creaexposure:::.concentration_path("pm25", "vandonkelaar", "v5", "test.tif")
  expect_true(grepl("pm25/vandonkelaar/v5/test.tif$", path))
})


# --- File template ------------------------------------------------------------

test_that("vandonkelaar v5 file template works", {
  config <- creaexposure:::.get_version_config("pm25", "vandonkelaar", "v5")
  expect_equal(
    config$file_template(2023),
    "V5GL0502.HybridPM25.Global.202301-202312.tif"
  )
})

test_that("vandonkelaar v6 file template works", {
  config <- creaexposure:::.get_version_config("pm25", "vandonkelaar", "v6")
  expect_equal(
    config$file_template(2022),
    "V6GL02.02.CNNPM25.Global.202201-202212.tif"
  )
})

test_that("vandonkelaar v4 file template works with variants", {
  config <- creaexposure:::.get_version_config("pm25", "vandonkelaar", "v4")
  expect_equal(
    config$file_template(2019),
    "ACAG_PM25_GWR_V4GL03_201901_201912_0p01.tif"
  )
  expect_equal(
    config$file_template(2019, variant = "no_ssdust"),
    "ACAG_PM25_noDUSTnoSEASALT_GWR_V4GL03_201901_201912_0p01.tif"
  )
})

test_that("larkin file template is static", {
  config <- creaexposure:::.get_version_config("no2", "larkin", "default")
  expect_equal(config$file_template(2011), "no2_agg10_ugm3_wsg84.tif")
  expect_equal(config$file_template(2023), "no2_agg10_ugm3_wsg84.tif")
})

test_that("geoschem O3 file template uses variant", {
  config <- creaexposure:::.get_version_config("o3", "geoschem", "default")
  expect_equal(config$file_template(variant = "m3m"), "o3_m3m.tif")
  expect_equal(config$file_template(variant = "sm8h"), "o3_sm8h.tif")
  expect_error(config$file_template(variant = "bad"), "Unknown O3 variant")
})

test_that("geoschem O3 default variant is m3m", {
  config <- creaexposure:::.get_version_config("o3", "geoschem", "default")
  expect_equal(config$default_variant, "m3m")
})

test_that("merra2 file template works", {
  config <- creaexposure:::.get_version_config("pm25", "merra2", "default")
  expect_equal(config$file_template(2020), "pm25_merra2_2020.tif")
})

test_that("tap file template works", {
  config <- creaexposure:::.get_version_config("pm25", "tap", "china")
  expect_equal(config$file_template(2023), "2023.tif")
})


# --- Year regex matching ------------------------------------------------------

test_that("vandonkelaar v5 year regex extracts year", {
  config <- creaexposure:::.get_version_config("pm25", "vandonkelaar", "v5")
  m <- stringr::str_match("V5GL0502.HybridPM25.Global.202301-202312.tif", config$year_regex)
  expect_equal(m[, 2], "2023")
})

test_that("vandonkelaar v6 year regex extracts year", {
  config <- creaexposure:::.get_version_config("pm25", "vandonkelaar", "v6")
  m <- stringr::str_match("V6GL02.02.CNNPM25.Global.202201-202212.tif", config$year_regex)
  expect_equal(m[, 2], "2022")
})

test_that("vandonkelaar v4 year regex extracts year", {
  config <- creaexposure:::.get_version_config("pm25", "vandonkelaar", "v4")
  m <- stringr::str_match("ACAG_PM25_GWR_V4GL03_201901_201912_0p01.tif", config$year_regex)
  expect_equal(m[, 2], "2019")
})

test_that("merra2 year regex extracts year", {
  config <- creaexposure:::.get_version_config("pm25", "merra2", "default")
  m <- stringr::str_match("pm25_merra2_2020.tif", config$year_regex)
  expect_equal(m[, 2], "2020")
})

test_that("omi year regex extracts year", {
  config <- creaexposure:::.get_version_config("no2", "omi", "default")
  m <- stringr::str_match("no2_omi_2019.tif", config$year_regex)
  expect_equal(m[, 2], "2019")
})


# --- Fixed years --------------------------------------------------------------

test_that("larkin has fixed years", {
  config <- creaexposure:::.get_version_config("no2", "larkin", "default")
  expect_equal(config$fixed_years, 2011L)
})

test_that("get_concentration_available_years returns fixed years for larkin", {
  years <- get_concentration_available_years("no2", source = "larkin")
  expect_equal(years, 2011L)
})


# --- Closest year -------------------------------------------------------------

test_that("closest year returns exact match when available", {
  # larkin only has 2011
  expect_equal(
    get_concentration_closest_year("no2", source = "larkin", year = 2011),
    2011L
  )
})

test_that("closest year returns nearest when no exact match", {
  # larkin only has 2011 — any year should return 2011
  expect_equal(
    get_concentration_closest_year("no2", source = "larkin", year = 2023),
    2011L
  )
})

test_that("closest year handles mid-year format", {
  expect_equal(
    get_concentration_closest_year("no2", source = "larkin", year = "mid2020mid2021"),
    2011L
  )
})
