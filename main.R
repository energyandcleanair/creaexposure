

# rasterOptions(tmpdir="/mnt/data/tmp/raster") # To avoid filling sda1
# terraOptions(tempdir="/mnt/data/tmp/raster") # To avoid filling sda1

# message("2pt5_min")
# res <- "2pt5_min"
# pop <- data.pop(res=res)

# pbapply::pblapply(seq(2018, 2022), function(year){
#   print(year)
#   adjust_global(pop, res, poll=c("pm25", "no2"), year=year)
# })

# Build more accurate version for certain countries --------------------------------------

lapply(seq(2020, 2023), function(year){
  aqexposuremap::adjust_global(res=aqexposuremap::RES_2PT5_MIN,
                               suffix="_china",
                               selected_regions = "CN",
                               year=year,
                               model=MODEL_RF)
})


aqexposuremap::adjust_country("30_sec", "_india", "IND", "IN", 2018)
aqexposuremap::adjust_country("30_sec", "_india", "IND", "IN", 2022)

aqexposuremap::adjust_country("30_sec", "_philippines", "PHL", "PH", 2018)
aqexposuremap::adjust_country("30_sec", "_philippines", "PHL", "PH", 2022)

aqexposuremap::adjust_country("30_sec", "_china", "CHN", "CN", 2020)
