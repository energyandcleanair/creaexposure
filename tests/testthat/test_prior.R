
library(testthat)

test_that("model prior returns prior", {
  year <- 2022
  res <- creaexposure::RES_2PT5_MIN
  bbox <- c(68.1, 6.7, 97.4, 35.5)
  bbox <- terra::ext(bbox, xy=T)
  pop <- creaexposure::data.pop(res=res, bbox=bbox)

  prior_direct <- data.basemap_pm25(pop=pop, res=res, year=year, use_cache=FALSE)

  prior_predicted <- creaexposure::build_map(res=res,
                                             pop=pop,
                                             polls = "pm25",
                                             model=MODEL_PRIOR,
                                             obs=tibble::tibble(),
                                             year=year,
                                             bbox=bbox,
                                             use_cache_predictors=FALSE)


  expect_equal(data.basemap_pm25_year(year), data.basemap_pm25_year(2023))
})
