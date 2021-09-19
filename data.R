#' Querying ground observations from CREADB
#'
#' @return
#' @export
#'
#' @examples
data.get_obs <- function(polls=c("pm25","no2"), year=2020, use_cache=T){

  f <- sprintf("cache/obs_%s_%s.RDS", year, paste(polls, collapse="_"))
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

data.basemap_pm25 <- function(pop, res, use_cache=T){

  f <- sprintf("cache/pm25_%s.tif", res)
  if(file.exists(f) && use_cache){
    terra::rast(f)
  }else{
    pm25_china <- data.basemap_pm25_region("China") %>% terra::resample(pop) %>% raster()
    pm25_eur <- data.basemap_pm25_region("Europe") %>% terra::resample(pop) %>% raster()
    pm25_usa <- data.basemap_pm25_region("NorthAmerica") %>% terra::resample(pop) %>% raster()
    pm25_global <- data.basemap_pm25_region("Global") %>% raster() %>%
      raster::resample(raster(pop)) #memory issue with terra

    pm25_stack <- stack(pm25_china, pm25_eur, pm25_usa, pm25_global)

    pm25 <- pm25_stack %>% raster::calc(
      fun=function(x, ...){
        local <- sum(x[1], x[2], x[3], na.rm=T)
        if(local==0) x[4] else local
      }
    )

    raster::writeRaster(pm25, f, overwrite=T)
    return(terra::rast(pm25))
  }
}

data.basemap_pm25_region <- function(region){

  region_fs <- list(
    "NorthAmerica" = "V4NA03_PM25_NA_201801_201812-RH35.nc",
    "Europe" = "V4EU03_PM25_EU_201801_201812-RH35.nc",
    "China" = "V4CH03_PM25_CHi_201801_201812-RH35.nc",
    "Global" = "ACAG_PM25_GWR_V4GL03_201801_201812_0p01.tif"
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
  pm25_map <- terra::rast(creahelpers::get_concentration_path(f_pm25_tif))
  crs(pm25_map) <- sp::CRS('+init=epsg:4326')
  return(pm25_map)
}


data.basemap_no2 <- function(pop, res, use_cache=T){

  f <- sprintf("cache/no2_%s.tif", res)
  if(file.exists(f) && use_cache){
    terra::rast(f)
  }else{
    f_no2 <- "no2_agg8.grd"
    f_no2_wsg84 <- "no2_agg8_wsg84.tif"

    if(!file.exists(creahelpers::get_concentration_path(f_no2_wsg84))){
      no2_wsg84 <- terra::rast(creahelpers::get_concentration_path(f_no2)) %>%
        terra::project("EPSG:4326") # Required for hia
      terra::writeRaster(no2_wsg84, creahelpers::get_concentration_path(f_no2_wsg84),
                         overwrite=T)
    }
    no2 <- terra::rast(creahelpers::get_concentration_path(f_no2_wsg84)) %>%
      terra::resample(pop)
    terra::writeRaster(no2, f)
    return(no2)
  }
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


#' Build a raster stack of all potential predictors used for model fitting
#'
#' @return
#' @export
#'
#' @examples
data.predictors <- function(pop, res){

  pm25.base <- data.basemap_pm25(pop, res, use_cache=T)
  no2.base <- data.basemap_no2(pop, res, use_cache=T)

  distance_coast <- data.distance_coast(pop, res, use_cache=T)
  # distance_urban <- data.distance_urban(pop, res, saturation_km=100, use_cache=T)
  grump <- data.grump(pop, res, use_cache=T)
  gadm0 <- data.gadm_raster(pop, res, level=0)
  gadm1 <- data.gadm_raster(pop, res, level=1)

  srtm <- data.srtm(pop, res, use_cache=T)
  srtm_05deg <- utils.focal_mean(srtm, d_deg=0.5, pop=pop, res=res, use_cache=T)
  # srtm_1deg <- utils.focal_mean(srtm, d_deg=1, pop=pop, res=res, use_cache=T)

  srtm_diff05deg <- srtm -srtm_05deg
  # srtm_diff1deg <- srtm -srtm_1deg

  # srtm_1deg <- utils.focal_mean(srtm, d_deg=1, res=res, use_cache=T)
  # srtm_2deg <- utils.focal_mean(srtm, d_deg=2, res=res, use_cache=T)

  pop_05deg <- utils.focal_mean(pop, d_deg=0.5, pop=pop, res=res, use_cache=T)
  pm25_ss_dust_frac <- data.pm25_ss_dust_frac(pop, res)
  lon <- data.lon(pop, res)
  lat <- data.lat(pop, res)

  # Not all predictors will be used but putting them together nonetheless
  predictors <- list(
    pm25_prior=pm25.base,
    no2_prior=no2.base,
    distance_coast=distance_coast,
    grump=grump,
    pop=pop,
    lon=lon,
    lat=lat,
    gadm0=gadm0,
    gadm1=gadm1,
    srtm=srtm,
    srtm_05deg=srtm_05deg,
    # srtm_1deg=srtm_1deg,
    srtm_diff05deg=srtm_diff05deg,
    # srtm_diff1deg=srtm_diff1deg,
    pop_05deg=pop_05deg,
    pm25_ss_dust_frac=pm25_ss_dust_frac
  )

  #rs: raster stack
  predictors <- lapply(predictors, raster) %>% raster::stack()
  names(predictors) <- names(predictors)

  # Adding surrogate variables
  # rs_predictors$distance_urban_inv <- 1/rs_predictors$distance_urban
  # rs_predictors$distance_urban_inv[is.infinite(rs_predictors$distance_urban_inv)] <- 1/400
  # rs_predictors$distance_urban_inv[rs_predictors$distance_urban_inv>1/400] <- 1/400
  # rs_predictors$road_density_log <- log(rs_predictors$road_density + 1)
  return(predictors)
}


data.grump <- function(pop, res, use_cache=T){

  f <- file.path("cache", paste0("grump_",res,".tif"))

  if(!use_cache | !file.exists(f)){
    # https://sedac.ciesin.columbia.edu/data/collection/grump-v1
    grump <- terra::rast(creahelpers::get_population_path('glurextents.bil'))
    n_aggregate <- floor(res(pop) / res(grump))
    if(n_aggregate>1){
      # We aggregate first to avoid memory issues
      grump <- grump %>% terra::aggregate(fact=n_aggregate, fun="modal")
    }
    grump.res <- grump %>% terra::resample(pop, method="near")
    terra::writeRaster(grump.res, f, overwrite=T)
    return(grump)
  }else{
    terra::rast(f)
  }
}

data.lon <- function(pop, res, use_cache=T){
  f <- file.path("cache", paste0("lon_",res,".tif"))
  if(!use_cache | !file.exists(f)){
    xy <- coordinates(raster(pop))
    r_lon <- raster(pop)
    r_lon[] <- xy[,1]
    r_lon <- r_lon %>% mask(raster(pop))
    terra::writeRaster(r_lon, f, overwrite=T)
  }else{
    terra::rast(f)
  }
}

data.lat <- function(pop, res, use_cache=T){
  f <- file.path("cache", paste0("lat_",res,".tif"))
  if(!use_cache | !file.exists(f)){
    xy <- coordinates(raster(pop))
    r_lat <- raster(pop)
    r_lat[] <- xy[,2]
    r_lat <- r_lat %>% mask(raster(pop))
    terra::writeRaster(r_lat, f, overwrite=T)
  }else{
    terra::rast(f)
  }
}

data.srtm <- function(pop, res, use_cache=T){

  f <- file.path("cache", paste0("srtm_",res,".tif"))

  if(!use_cache | !file.exists(f)){

    srtm <- terra::rast(creahelpers::get_elevation_path('SRTM_1km_GRD/srtmv4_30s/w001001.adf')) %>%
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


#'
#' #' Download roads from SEDAC
#' #' https://sedac.ciesin.columbia.edu/data/set/groads-global-roads-open-access-v1/data-download
#' #'
#' #' @param use_cache
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' data.road_density_groads <- function(res, pop, use_cache=T){
#'
#'   f <- file.path("cache", paste0("road_density_groads_",res,".tif"))
#'
#'   if(!use_cache | !file.exists(f)){
#'
#'     res_rasterizing <- "2pt5_min"
#'     pop_rasterizing <- data.pop(res=res_rasterizing)
#'
#'
#'     roads <- sf::read_sf(file.path(creahelpers::get_gis_dir(), "roads",
#'                                    "gROADS_v1.gdb"))
#'
#'     # roads <- roads %>% as("Spatial")
#'
#'     # r <- creahelpers::rasterize_lines(lines=roads,
#'     #                                   grid=raster(pop_rasterizing))
#'
#'     lines <- roads
#'     grid <- pop
#'     raster::crs(lines) <- raster::crs(raster::raster(pop))
#'
#'     # Cut lines along grid cells
#'     print("Polygonizing...")
#'     rs <- grid
#'     rs[] <- 1:ncell(rs)
#'     names(rs) <- "i_cell"
#'     rsp <- terra::as.polygons(rs)
#'     rsp <- as(rsp, "Spatial")
#'     print("Done")
#'
#'
#'     # Add temporary feature id for grouping
#'     lines$feature_id_tmp <- 1:nrow(lines)
#'
#'     print("Cutting along grid...")
#'     # sf much less memory intensive than raster::intersect
#'     # and faster
#'
#'     # Chunking it to avoid rgeos_binpredfunc_prepared: maximum returned dense matrix size exceeded
#'     cutting_successful <- F
#'     chunk_size <- 1E10
#'     while(!cutting_successful){
#'       tryCatch({
#'         rsp$chunk <- rsp$i_cell %/% chunk_size
#'         emission.sf <- sf::st_as_sf(lines)
#'         rp <- pbapply::pblapply(split(sf::st_as_sf(rsp), rsp$chunk),
#'                                 function(rsp_chunk){
#'                                   sf::st_intersection(emission.sf,rsp_chunk)
#'                                 }) %>%
#'           do.call("bind_rows",.)
#'         cutting_successful <- T
#'       }, error=function(e){
#'         if("size exceeded" %in% as.character(e)){
#'           chunk_size <- chunk_size / 100
#'           warning("Cutting failed: ", e, "\n Trying with smaller chunk size", )
#'         }else{
#'           stop(e)
#'         }
#'       })
#'     }
#'
#'     print("Done")
#'
#'     print("Calculating length...")
#'     rp$length <- sf::st_length(rp, byid=TRUE)
#'     print("Done")
#'
#'     # # Weighting accordingly
#'     # print("Weighting by length...")
#'     # rp <- rp %>%
#'     #   group_by(feature_id_tmp) %>%
#'     #   do(mutate(., emission=.$emission * length / sum(.$length)))
#'     # print("Done")
#'
#'     print("Rasterizing...")
#'     rp.sum <- rp %>%
#'       group_by(i_cell=as.integer(i_cell)) %>%
#'       summarise(length=sum(length, na.rm=T))
#'
#'     # Print into raster directly!
#'     cells_x <- rep(0,ncell(rs))
#'     cells_x[rp.sum$i_cell] <- rp.sum$length
#'     grid_result <- grid
#'     grid_result[] <- cells_x
#'
#'
#'
#'
#'
#'
#'
#'
#'     writeRaster(r, f)
#'     return(r)
#'   }else{
#'     raster(f)
#'   }
#'
#'
#' }



#' Download roads from GRIP4
#'
#' @param use_cache
#'
#' @return
#' @export
#'
#' @examples
data.road_density_grip <- function(res, pop, use_cache=T){

  f <- file.path("cache", paste0("road_density_grip_",res,".tif"))

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
    grump <- data.grump(pop, res)
    sea_level <- 0
    rural_level <- 1
    urban_level <- 2

    # Avoid computing distance in sea
    grump[grump==sea_level] <- -1
    grump[grump==rural_level] <- NA
    dist <- terra::distance(grump)
    dist[is.na(pop)] <- NA
    dist[dist==-1] <- NA
    raster::writeRaster(dist, f, overwrite=T)
  }else{
    dist <- terra::rast(f)
  }
  dist[dist>saturation_km*1000] <- saturation_km*1000
  return(dist)
}


data.gadm_raster <- function(pop, res, level, use_cache=T){

  f <- file.path("cache", paste0(sprintf("gadm%d_%s.tif", level, res)))

  if(!use_cache | !file.exists(f)){
    g <- creahelpers::get_adm(level, "full")
    gid_level <- sprintf("GID_%d",level)

    # gid <- sf::st_as_sf(g[gid_level]) %>%
    #   rename_at(gid_level, function(x)"gid") %>%
    #   mutate(gid=as.numeric(factor(gid)))

    r <- terra::rasterize(terra::vect(g), pop,
                          field=gid_level,
                          fun=first)

    # Extending a bit into water to be sure
    # we're not missing coastal cities
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

    r <- raster(r)
    w <- raster::focalWeight(r, 0.3, "rectangle")
    w[] <- 1
    r_filled <- creahelpers::focal.loop(r, w, fill.na, NAonly=T)

    values(r_filled)[is.na(values(raster(pop)))] <- NA
    raster::writeRaster(r_filled, f, overwrite=T)
    return(terra::rast(r_filled))
  }else{
    terra::rast(f)
  }
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
      terra::resample(pop)

    nasadist[nasadist>0]=0
    nasadist=-nasadist
    terra::writeRaster(nasadist, f, overwrite=T)
    return(nasadist)
  }else{
    terra::rast(f)
  }
}



#' Seasalt and dust proportion
#'
#' @param pop
#' @param res
#' @param use_cache
#'
#' @return
#' @export
#'
#' @examples
data.pm25_ss_dust_frac <- function(pop, res, use_cache=T){

  f <- sprintf("cache/pm25_ss_dust_frac_%s.tif", res)
  if(file.exists(f) && use_cache){
    terra::rast(f)
  }else{
    pm25_w_ssd <- raster(creahelpers::get_concentration_path("ACAG_PM25_GWR_V4GL03_201801_201812_0p01.tif"))
    pm25_wo_ssd <- raster(creahelpers::get_concentration_path("ACAG_PM25_noDUSTnoSEASALT_GWR_V4GL03_201801_201812_0p01.tif"))
    ss_dust_frac <- (pm25_w_ssd - pm25_wo_ssd)/ pm25_w_ssd
    ss_dust_frac <- ss_dust_frac %>% raster::resample(raster(pop))
    raster::writeRaster(ss_dust_frac, f, overwrite=T)
    return(ss_dust_frac)
  }
}
