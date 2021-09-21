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

rasterOptions(tmpdir="/mnt/data/tmp/raster") # To avoid filling sda1
terraOptions(tempdir="/mnt/data/tmp/raster") # To avoid filling sda1

source('data.R')
source('utils.R')
source('adjust_global.R')

message("2pt5_min")
res <- "2pt5_min"
pop <- data.pop(res=res)
adjust_global(pop, res, poll="no2")
