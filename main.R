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

adjust_global(pop, res, poll="no2")



# Build more accurate version for PH --------------------------------------
res_ph = "30_sec"
bb_ph <- data.gadm0() %>% subset(GID_0=="PHL") %>% bbox()
pop_ph <- data.pop(res=res_ph, bb=bb_ph)
adjust_global(pop_ph, res_ph, suffix="_philippines", selected_regions="PH")
