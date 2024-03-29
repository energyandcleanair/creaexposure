

# rasterOptions(tmpdir="/mnt/data/tmp/raster") # To avoid filling sda1
# terraOptions(tempdir="/mnt/data/tmp/raster") # To avoid filling sda1

# message("2pt5_min")
# res <- "2pt5_min"
# pop <- data.pop(res=res)

# pbapply::pblapply(seq(2018, 2022), function(year){
#   print(year)
#   adjust_global(pop, res, poll=c("pm25", "no2"), year=year)
# })

map_in <- creaexposure::build_map(res=creaexposure::RES_2PT5_MIN,
                                   suffix="_india",
                                   selected_regions = "IN",
                                   year=2022,
                                   remove_seasalt_dust_contribution = F,
                                   limit_distance_urban = F,
                                   obs_level = "station",
                                   model=MODEL_RF)



map_in <- creaexposure::build_map(res=creaexposure::RES_30_SEC,
                                   suffix="_india",
                                   selected_regions = "IN",
                                   year=2022,
                                   remove_seasalt_dust_contribution = F,
                                   limit_distance_urban = F,
                                   obs_level = "station",
                                   model=MODEL_RF)

a <- get_exposure_map(year=2020,
                 res=creaexposure::RES_2PT5_MIN,
                 model="rf",
                 limit_distance_urban=T,
                 remove_seasalt_dust_contribution=T,
                 results_folder = glue('results/{model}/'),
                 suffix = "test")



map_cn <- creaexposure::build_map(res=creaexposure::RES_2PT5_MIN,
                         suffix="_china",
                         selected_regions = "CN",
                         year=2023,
                         remove_seasalt_dust_contribution = F,
                         limit_distance_urban = F,
                         obs_level = "station",
                         model=MODEL_RF)

# Build more accurate version for certain countries --------------------------------------

lapply(seq(2020, 2023), function(year){
  creaexposure::build_map(res=creaexposure::RES_30_SEC,
                               suffix="_china",
                               selected_regions = "CN",
                               year=year,
                               model=MODEL_GAM)
})

lapply(seq(2020, 2023), function(year){
  creaexposure::build_map(res=creaexposure::RES_2PT5_MIN,
                               suffix="_china",
                               selected_regions = "CN",
                               year=year,
                               model=MODEL_RF)
})



creaexposure::adjust_country("30_sec", "_india", "IND", "IN", 2018)
creaexposure::adjust_country("30_sec", "_india", "IND", "IN", 2022)

creaexposure::adjust_country("30_sec", "_philippines", "PHL", "PH", 2018)
creaexposure::adjust_country("30_sec", "_philippines", "PHL", "PH", 2022)

creaexposure::adjust_country("30_sec", "_china", "CHN", "CN", 2020)
