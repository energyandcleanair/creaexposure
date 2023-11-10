library(raster)
library(terra)
library(tidyverse)
library(rcrea)
library(creahelpers)
library(mgcv)
library(tictoc)
library(snow)
library(tictoc)
library(countrycode)
library(pbapply)
library(glue)

readRenviron(".Renviron")

dir.create("cache", showWarnings = F)
dir.create("results", showWarnings = F)

# rasterOptions(tmpdir="/mnt/data/tmp/raster") # To avoid filling sda1
# terraOptions(tempdir="/mnt/data/tmp/raster") # To avoid filling sda1

source('data.R')
source('utils.R')
source('adjust_global.R')

message("2pt5_min")
res <- "2pt5_min"
pop <- data.pop(res=res)

pbapply::pblapply(seq(2018, 2022), function(year){
  print(year)
  adjust_global(pop, res, poll=c("pm25", "no2"), year=year)
})



# Build more accurate version for certain countries --------------------------------------
res_in = "30_sec"
bb_in <- data.gadm0() %>% subset(GID_0=="IND") %>% sf::st_as_sf() %>% sf::st_bbox()
pop_in <- data.pop(res=res_in, bb=bb_in)
adjust_global(pop_in, res_in, suffix="_india", selected_regions="IN", use_cache=T, year=2018)
adjust_global(pop_in, res_in, suffix="_india", selected_regions="IN", use_cache=T, year=2022)


res_ph = "30_sec"
bb_ph <- data.gadm0() %>% subset(GID_0=="PHL") %>% sf::st_as_sf() %>% sf::st_bbox()
pop_ph <- data.pop(res=res_ph, bb=bb_ph)
adjust_global(pop_ph, res_ph, suffix="_philippines", selected_regions="PH", year=2018)
adjust_global(pop_ph, res_ph, suffix="_philippines", selected_regions="PH", year=2022)
