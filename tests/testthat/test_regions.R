library(testthat)

test_that("Basic regions work", {


  regions <- get_all_regions()
  expect_equal(get_region_iso2s("CN"), "CN")
  expect_equal(get_region_iso2s("XX"), "XX")
  expect_equal(get_region_iso2s(list("myregion"=c("CN", "IN"))), c("CN", "IN"))


})


test_that("Add region to measurements works", {

  meas <- tibble(country = c("IN", "CN", "DE", "XX"), value = c(1, 2, 3, 4))

  # User specified regions
  regions <- list("asia" = c("IN", "CN"), "europe" = c("DE", "FR"))
  meas_w_obs <- add_region_to_obs(meas, regions)
  expect_equal(meas_w_obs$region, c("asia", "asia", "europe", NA))

  # Default
  meas_w_obs_default <- add_region_to_obs(meas)
  expect_equal(meas_w_obs$region, c("IN", "CN", "EU", NA))
})


test_that("Get bbox works", {

  bbox_cn1 <- get_bbox("CN")
  bbox_cn2 <- get_bbox(list("myregion"=c("CN")))
  expect_equal(bbox_cn1, bbox_cn2)

  bbox_asia <- get_bbox(list("myregion"=c("RU", "IN", "JP")))
  # Test that bbox_cn1 is contained within bbox_cn_in but not the other way around
  expect_true(all(bbox_cn1[1:2] > bbox_asia[1:2]))
  expect_true(all(bbox_cn1[3:4] < bbox_asia[3:4]))

})
