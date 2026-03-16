# testthat::source_test_helpers("tests", env = globalenv())
# testthat::source_test_helpers("../", env = globalenv())

library(testthat)

test_that("Mid-year works", {
  year <- "mid2023mid2024"
  available <- get_concentration_available_years("pm25", source = "vandonkelaar")
  expect_true(length(available) > 0, "No vandonkelaar PM2.5 files available on disk")
  expect_equal(
    get_concentration_closest_year("pm25", source = "vandonkelaar", year = year),
    get_concentration_closest_year("pm25", source = "vandonkelaar", year = 2023)
  )
})


test_that("Predictors exist mid-year", {
  year <- "mid2023mid2024"
  path <- creahelpers::get_concentration_path(
    file.path("pm25", "merra2", "default", sprintf("pm25_merra2_%s.tif", year))
  )
  expect_true(file.exists(path), "MERRA2 mid-year file not available on disk")
})
