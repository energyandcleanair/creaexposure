#' Querying ground observations from CREADB
#'
#' @return
#' @export
#'
#' @examples
data.get_obs <- function(polls=c("pm25","no2"), year=2020, use_cache=T){

  f <- sprintf("cache/pm_obs_%s_%s.RDS", year, paste(polls, collapse="_"))
  if(file.exists(f) && use_cache){
    return(readRDS(f))
  }else{
    min_measurements <- 365 * 0.8
    obs <- rcrea::measurements(poll=c("pm10","pm25","no2"),
                               with_geometry = T,
                               with_metadata = T,
                               process_id="city_day_mad",
                               collect=F,
                               date_from=paste0(year, "-01-01"),
                               date_to=paste0(year, "-12-31")) %>%
      dplyr::group_by(location_id, country, poll, unit, source, geometry) %>%
      dplyr::summarise(value=mean(value, na.rm=T),
                       count=n()) %>%
      dplyr::collect() %>%
      # Filter
      filter(count >= min_measurements) %>%
      # Keep the source with the most measurements
      group_by(location_id, poll) %>%
      arrange(desc(count)) %>%
      slice(1) %>%
      ungroup() %>%
      dplyr::mutate(geometry = sf::st_as_sfc(geometry))

    saveRDS(obs, f)
    return(obs)
  }
}



data.basemap_pm25 <- function(region){

  region_fs <- list(
    "NorthAmerica" = "V4NA03_PM25_NA_201801_201812-RH35.nc",
    "Europe" = "V4EU03_PM25_EU_201801_201812-RH35.nc",
    "China" = "V4CH03_PM25_CHi_201801_201812-RH35.nc",
    "Global" = "ACAG_PM25_GWR_V4GL03_201901_201912_0p01.nc" #TODO use 2018
  )

  if(!region %in% names(region_fs)){
    stop("Region should be in ", paste(fs, collapse=", "))
  }

  f_pm25_nc <- region_fs[[region]]
  f_pm25_tif <- gsub("\\.nc","\\.tif", f_pm25_nc)

  nc_to_tif <- function(f_pm25_nc, f_pm25_tif){
    pm25.base <- raster::raster(creahelpers::get_concentration_path(f_pm25_nc))

    if(region=="Global"){
      # For some reason the Global dataset is flipped
      pm25.base <- pm25.base %>%
        raster::t() %>% # nc dimensions aren't in raster expected order I suppose
        raster::flip(direction="y") %>%
        raster::flip(direction="x")
    }

    raster::writeRaster(pm25.base,
                        creahelpers::get_concentration_path(f_pm25_tif),
                        overwrite=T)
  }

  if(!file.exists(creahelpers::get_concentration_path(f_pm25_tif))) nc_to_tif(f_pm25_nc, f_pm25_tif)

  message("Using PM2.5 from ", f_pm25_tif)
  pm25_map <- raster(creahelpers::get_concentration_path(f_pm25_tif))
  crs(pm25_map) <- sp::CRS('+init=epsg:4326')
  return(pm25_map)
}


data.basemap_no2 <- function(){
  f_no2 <- "no2_agg8.grd"
  f_no2_wsg84 <- "no2_agg8_wsg84.tif"

  if(!file.exists(creahelpers::get_concentration_path(f_no2_wsg84))){
    no2_wsg84 <- terra::rast(creahelpers::get_concentration_path(f_no2)) %>%
      terra::project("EPSG:4326") # Required for hia
    terra::writeRaster(no2_wsg84, creahelpers::get_concentration_path(f_no2_wsg84),
                       overwrite=T)
  }
  raster(creahelpers::get_concentration_path(f_no2_wsg84))
}


data.gadm0 <- function(){
  creahelpers::get_adm(level=0, res="coarse")
}

data.mask <- function(){
  m <- creahelpers::get_adm(level=0, res="full")
  m[m$GID_0=="PHL",]
}

data.bb <- function(){
  data.gadm0() %>%
    subset(GID_0 %in% c('PHL')) %>% extent %>% multiply_by(1.1)
}

data.pop <- function(res="30_sec", bb=NULL, mask=NULL){
  r <- terra::rast(creahelpers::get_population_path(
    paste0('gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2020_',res,'.tif')))

  if(!is.null(bb)){
    r <- r %>% terra::crop(bb)
  }

  if(!is.null(mask)){
    r <- terra::mask(r, terra::vect(mask))
  }
  names(r) <- "gpw"
  r
}

data.grump <- function(pop, res, use_cache=T){

  f <- file.path("cache", paste0("grump_",res,".tif"))

  if(!use_cache | !file.exists(f)){
    # https://sedac.ciesin.columbia.edu/data/collection/grump-v1
    grump <- terra::rast(creahelpers::get_population_path('glurextents.bil')) %>%
      terra::resample(pop)
    terra::writeRaster(grump, f, overwrite=T)
    return(grump)
  }else{
    terra::rast(f)
  }
}


data.srtm <- function(pop, res, use_cache=T){

  f <- file.path("cache", paste0("srtm_",res,".tif"))

  if(!use_cache | !file.exists(f)){

    srtm <- terra::rast(creahelpers::get_population_path('srtm.....')) %>%
      terra::resample(pop)

    names(srtm) <- "srtm"

    # Some coastal pixels aren't covered by srtm but are by pop
    # There's a risk of missing significant part of the population as
    # pixel with NA values are excluded from GWR
    values(srtm)[is.na(values(srtm)) & !is.na(values(pop))] <- 0
    terra::writeRaster(srtm, f, overwrite=T)
    return(srtm)
  }else{
    terra::rast(f)
  }
}

#' Difference between local elevation and background (say 1deg around) elevation
#' This has been used by Donkelaar https://pubs.acs.org/doi/suppl/10.1021/acs.est.0c01764/suppl_file/es0c01764_si_001.pdf
#'
#' @param pop
#' @param res
#' @param use_cache
#'
#' @return
#' @export
#'
#' @examples
data.srtm_mean <- function(pop, res, d_deg=1, use_cache=T){

  name <- sprintf("srtm_mean_%sdeg", d_deg)
  f <- file.path("cache", sprintf("srtm_mean_%sdeg_%s.tif", d_deg, res))

  if(!use_cache | !file.exists(f)){

    srtm <- data.srtm(pop, res)
    w <- raster::focalWeight(raster(srtm), d_deg, "circle")
    w[w>0] <- 1
    srtm_background <- raster::focal(raster(srtm),
                                     w,
                                     na.rm=T,
                                     pad=T,
                                     fun=mean) %>%
      terra::rast()

    srtm_background[is.na(pop)] <- NA
    names(srtm_background) <- name
    terra::writeRaster(srtm_background, f, overwrite=T)
    return(srtm_background)
  }else{
    terra::rast(f)
  }
}





#' Download roads from OpenStreetMap
#'
#' @param use_cache
#'
#' @return
#' @export
#'
#' @examples
data.road_density <- function(res, pop, use_cache=T){

  f <- file.path("cache", paste0("road_density_",res,".tif"))

  if(!use_cache | !file.exists(f)){

    urls <- c("data/roads/GRIP4_density_tp1.zip"="https://dataportaal.pbl.nl/downloads/GRIP4/GRIP4_density_tp1.zip",
              "data/roads/GRIP4_density_tp2.zip"="https://dataportaal.pbl.nl/downloads/GRIP4/GRIP4_density_tp2.zip")

    lapply(names(urls), function(f){
      if(!file.exists(f)){
        download.file(url=urls[[f]], destfile=f)
        unzip(f, exdir = "data/roads")
      }
    })

    r <- raster("data/roads/grip4_tp1_dens_m_km2.asc") +
      raster("data/roads/grip4_tp2_dens_m_km2.asc")

    writeRaster(r, f)
    return(r)
  }else{
    raster(f)
  }
}


data.distance_urban <- function(pop, res, saturation_km=100, use_cache=T){

  f <- file.path("cache", paste0("distance_urban_",res,".tif"))

  if(!use_cache | !file.exists(f)){

    urban <- terra::rast(creahelpers::get_landcover_path("C3S-LC-L4-LCCS-Map-300m-P1Y-2019-v2.1.1.nc"))$lccs_class
    urban[is.na(urban)] <- 0
    urban_level <- 190
    urban[urban!=urban_level] <- NA
    dist <- terra::distance(urban)
    dist <- dist %>% terra::resample(pop)
    dist[is.na(pop)] <- NA
    raster::writeRaster(dist, f, overwrite=T)
  }else{
    dist <- terra::rast(f)
  }

  dist[dist>saturation_km*1000] <- saturation_km*1000
  return(dist)
}



data.landuse <- function(pop, res, use_cache=T){

  f <- file.path("cache", paste0("landuse_",res,".tif"))

  if(!use_cache | !file.exists(f)){
    landuse <- terra::rast(creahelpers::get_landcover_path("C3S-LC-L4-LCCS-Map-300m-P1Y-2019-v2.1.1.nc"))$lccs_class %>%
      terra::crop(extent(raster(pop))*1.1)


    # the original 22 classes in the product were integrated into eight
    # classes including farmland, woodland, grassland, sparse vegetation, bare land, urban area, water body, and snow and ice
    names(landuse) <- "landuse"

    # Reclassifying
    l <- data.landuse_factors()

    reclas_mat <- tibble(left=as.numeric(names(l)),
                         right=as.numeric(factor(unlist(l)))) %>%
      as.matrix()

    landuse <- landuse %>% terra::classify(reclas_mat)
    n_aggregate <- res(pop) / res(landuse)
    landuse <- landuse %>% terra::aggregate(fact=n_aggregate, fun=modal)
    landuse2 <- landuse %>% terra::resample(pop, method="near")

    values(landuse)[is.na(values(raster(pop)))] <- NA

    # Extend a bit within sea to have all stations / population covered
    # Get mode
    fill.na <- function(x, i=5, ...) {
      v <- x[!is.na(x)]
      if(length(v)>0){
        getmode <- function(v) {
          uniqv <- unique(v)
          uniqv[which.max(tabulate(match(v, uniqv)))]
        }
        return(getmode(v))
      }else{
        return(NA)
      }
    }

    w <- raster::focalWeight(raster(landuse), 0.3, "rectangle")
    w[] <- 1
    landuse <- raster::focal(raster(landuse),
                             w,
                             NAonly=T,
                             fun=fill.na)

    values(landuse)[is.na(values(raster(pop)))] <- NA
    raster::writeRaster(landuse, f, overwrite=T)
    return(terra::rast(landuse))
  }else{
    terra::rast(f)
  }
}

data.landuse_factors <- function(as_factor=F){
  # Original classification
  # Only 1-16 in PH
  # 0 No Data
  # 10 Cropland, rainfed
  # 20 Cropland, irrigated or post-flooding
  # 30 Mosaic cropland (>50%) / natural vegetation (tree, shrub, herbaceous cover) (<50%)
  # 40 Mosaic natural vegetation (tree, shrub, herbaceous cover) (>50%) / cropland (<50%)
  # 50 Tree cover, broadleaved, evergreen, closed to open (>15%)
  # 60 Tree cover, broadleaved, deciduous, closed to open (>15%)
  # 70 Tree cover, needleleaved, evergreen, closed to open (>15%)
  # 80 Tree cover, needleleaved, deciduous, closed to open (>15%)
  # 90 Tree cover, mixed leaf type (broadleaved and needleleaved)
  # 100 Mosaic tree and shrub (>50%) / herbaceous cover (<50%)
  # 110 Mosaic herbaceous cover (>50%) / tree and shrub (<50%)
  # 120 Shrubland
  # 130 Grassland
  # 140 Lichens and mosses
  # 150 Sparse vegetation (tree, shrub, herbaceous cover) (<15%)
  # 160 Tree cover, flooded, fresh or brakish water
  # 170 Tree cover, flooded, saline water
  # 180 Shrub or herbaceous cover, flooded, fresh/saline/brakish water
  # 190 Urban areas
  # 200 Bare areas
  # 210 Water bodies
  # 220 Permanent snow and ice

  l <- list(
    "10"="Cropland",
    "11"="Cropland",
    "12"="Cropland",
    "20"="Cropland",
    "30"="Cropland",
    "40"="Woodland",
    "50"="Woodland",
    "70"="Woodland", #Not in PH, beyond extent
    "80"="Woodland",
    "100"="Woodland",
    "110"="Woodland",
    "120"="Grassland",
    "121"="Grassland",
    "130"="Grassland",
    "150"="Sparse",
    "160"="Woodland",
    "170"="Woodland",
    "180"="Grassland", #Not in PH, beyond extent
    "190"="Urban",
    "210"=NA#"Water"
  )

  if(as_factor){
    return(factor(unlist(l)))
  }else{
    return(l)
  }
}


data.distance_coast <- function(pop, res, use_cache=T){

  f <- file.path("cache", paste0("distance_coast_",res,".tif"))
  if(!use_cache | !file.exists(f)){
    # Using NASA dist from coast
    # https://oceancolor.gsfc.nasa.gov/docs/distfromcoast/
    expand <- 1
    nasadist <- terra::rast(creahelpers::get_boundaries_path("distfromcoast/GMT_intermediate_coast_distance_01d.tif")) %>%
      terra::crop(extent(projectExtent(raster::raster(pop),
                                       crs(., proj4=T)))+expand) %>%
      terra::project(pop)

    nasadist[nasadist>0]=0
    nasadist=-nasadist
    terra::writeRaster(nasadist, f, overwrite=T)
    return(nasadist)
  }else{
    terra::rast(f)
  }
}

