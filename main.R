library(raster)
library(terra)
library(tidyverse)
library(rcrea)
library(creahelpers)
readRenviron(".Renviron")

source('data.R')

obs.2018 <- data.get_obs(year=2019)
obs.2019 <- data.get_obs(year=2019)
obs.2020 <- data.get_obs(year=2020)

res <- "2pt5_min"
pop <- data.pop(res=res)


