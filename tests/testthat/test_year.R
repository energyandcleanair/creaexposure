# testthat::source_test_helpers("tests", env = globalenv())
# testthat::source_test_helpers("../", env = globalenv())

library(testthat)

test_that("Mid-year works", {
  year <- "mid2023mid2024"
  expect_equal(data.basemap_pm25_year(year), data.basemap_pm25_year(2023))
})


test_that("Predictors exist mid-year", {
  year <- "mid2023mid2024"
  expect_true(file.exists(creahelpers::get_concentration_path(sprintf("pm25_merra2_%s.tif", year))))
})
