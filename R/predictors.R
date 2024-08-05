#' Querying ground observations from CREADB
#'
#' @return
#' @export
#'
#' @examples
data.get_obs <- function(level = "city", polls=c("pm25","no2"), year=2020, use_cache=T){

  f <- glue("cache/obs_{level}_{year}_{paste(polls, collapse='_')}.RDS")
  dir.create(dirname(f), F, T)
  if(file.exists(f) && use_cache){
    return(readRDS(f))
  }else{
    min_measurements <- 365 * 0.8
    process_id <- dplyr::recode(level, "city"="city_day_mad", "station"="station_day_mad")
    obs <- rcrea::measurements(poll=c("pm10","pm25","no2"),
                               with_geometry = T,
                               with_metadata = T,
                               process_id = process_id,
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

    # We add Philipines data assuming year doesn't matter
    obs_ph <- read_csv(utils.get_data_file("ph_emb_derived_pm25.csv")) %>%
      mutate(location_id=MonitorID,
             country="PH",
             poll="pm25",
             unit="µg/m3",
             source="em",
             type=PM2.5_inferred_basis,
             value=pm25) %>%
      sf::st_as_sf(coords=c("Lon","Lat"), crs=4326) %>%
      select(location_id, country, poll, unit, source, type, value, geometry) %>%
      as.data.frame()

    obs <- bind_rows(tibble(obs), tibble(obs_ph))
    obs$type <- factor(obs$type)


    # NO2: 1 ppb = 1.88 µg/m3
    idx_ppm_to_ugm3 <- obs$unit=="ppm" & obs$poll=="no2"
    obs[idx_ppm_to_ugm3, "value"] <- obs[idx_ppm_to_ugm3, "value"] * 1880
    obs[idx_ppm_to_ugm3, "unit"] <- "µg/m3"

    idx_ppb_to_ugm3 <- obs$unit=="ppb" & obs$poll=="no2"
    obs[idx_ppb_to_ugm3, "value"] <- obs[idx_ppb_to_ugm3, "value"] * 1.88
    obs[idx_ppb_to_ugm3, "unit"] <- "µg/m3"

    # A weird ZA measurement
    obs=obs[obs$unit!="m/s",]

    saveRDS(obs, f)
    return(obs)
  }
}

#' Build a raster stack of all potential predictors used for model fitting
#'
#' @return
#' @export
#'
#' @examples
data.predictors <- function(pop, res, year, use_cache=T, suffix=""){

  pm25.base <- data.basemap_pm25(pop, res, year=year, use_cache=use_cache, suffix=suffix)
  no2.base <- data.basemap_no2(pop, res, use_cache=use_cache, suffix=suffix)

  distance_coast <- data.distance_coast(pop, res, use_cache=use_cache, suffix=suffix)
  distance_urban <- data.distance_urban(pop, res, use_cache=use_cache, suffix=suffix)
  grump <- data.grump(pop, res, use_cache=use_cache, suffix=suffix)
  gadm0 <- data.gadm_raster(pop, res, level=0, use_cache=use_cache, suffix=suffix)
  gadm1 <- data.gadm_raster(pop, res, level=1, use_cache=use_cache, suffix=suffix)

  pm25_merra2_diff <- data.pm25_merra2_diff(pop, res,
                                       year_i=data.basemap_pm25_year(year),
                                       year_f=year,
                                       use_cache=use_cache, suffix=suffix)
  no2_omi_diff <- data.no2_omi_diff(pop, res, year_f=year, use_cache=use_cache, suffix=suffix)
  pop_ratio_log <- data.pop_ratio_log(pop, res, use_cache=use_cache, suffix=suffix)
  srtm <- data.srtm(pop, res, use_cache=use_cache, suffix=suffix)
  srtm_05deg <- utils.focal_mean(srtm, d_deg=0.5, pop=pop, res=res, use_cache=use_cache, suffix=suffix)
  # srtm_1deg <- utils.focal_mean(srtm, d_deg=1, pop=pop, res=res, use_cache=use_cache)
  # type <- data.type(pop=pop)

  srtm_diff05deg <- srtm -srtm_05deg
  # srtm_diff1deg <- srtm -srtm_1deg

  # srtm_1deg <- utils.focal_mean(srtm, d_deg=1, res=res, use_cache=use_cache)
  # srtm_2deg <- utils.focal_mean(srtm, d_deg=2, res=res, use_cache=use_cache)

  pop_05deg <- utils.focal_mean(pop, d_deg=0.5, pop=pop, res=res, use_cache=use_cache, suffix=suffix)
  pm25_ss_dust_frac <- data.pm25_ss_dust_frac(pop, res, use_cache=use_cache, suffix=suffix)
  pm25_ss_dust <- data.pm25_ss_dust(pop, res, use_cache=use_cache, suffix=suffix)
  lon <- data.lon(pop, res, use_cache=use_cache, suffix=suffix)
  lat <- data.lat(pop, res, use_cache=use_cache, suffix=suffix)

  # Not all predictors will be used but putting them together nonetheless
  predictors <- list(
    pm25_prior=pm25.base,
    no2_prior=no2.base,
    distance_coast=distance_coast,
    distance_urban=distance_urban,
    grump=grump,
    pop=pop,
    lon=lon,
    lat=lat,
    gadm0=gadm0,
    gadm1=gadm1,
    no2_omi_diff=no2_omi_diff,
    pm25_merra2_diff=pm25_merra2_diff,
    pop_ratio_log=pop_ratio_log,
    srtm=srtm,
    srtm_05deg=srtm_05deg,
    # srtm_1deg=srtm_1deg,
    srtm_diff05deg=srtm_diff05deg,
    # srtm_diff1deg=srtm_diff1deg,
    pop_05deg=pop_05deg,
    pm25_ss_dust_frac=pm25_ss_dust_frac,
    pm25_ss_dust=pm25_ss_dust
    # type=type
  )

  #rs: raster stack
  predictors <- creahelpers::to_raster(predictors) %>% raster::stack()

  # Adding surrogate variables
  # rs_predictors$distance_urban_inv <- 1/rs_predictors$distance_urban
  # rs_predictors$distance_urban_inv[is.infinite(rs_predictors$distance_urban_inv)] <- 1/400
  # rs_predictors$distance_urban_inv[rs_predictors$distance_urban_inv>1/400] <- 1/400
  # rs_predictors$road_density_log <- log(rs_predictors$road_density + 1)
  return(predictors)
}


#' Find the closest year for pm2.5
#'
#' @param year either the year integer, or mid2023mid2024 if from July to June
#'
#' @return
#' @export
#'
#' @examples
data.basemap_pm25_year <- function(year){

  if(str_detect(year, "^mid\\d{4}mid\\d{4}$")){
    year <- as.numeric(str_match(year, "mid(\\d{4})mid")[2])
  }else{
    year <- as.numeric(year)
  }

  basemap_years <- seq(2018, 2022)
  basemap_year <- max(basemap_years[basemap_years<=year])
  return(basemap_year)
}


data.basemap_pm25 <- function(pop, res, year=2020, use_cache=T, suffix=""){

  basemap_year <- data.basemap_pm25_year(year)

  f <- sprintf("cache/pm25_%s_%s%s.tif", basemap_year, res, suffix)
  if(file.exists(f) && use_cache){
    terra::rast(f)
  }else{
    pm25 <- data.basemap_pm25_region("Global", year=basemap_year) %>% terra::resample(pop, method='bilinear')
    terra::writeRaster(pm25, filename = f, overwrite = T)
    return(pm25)
  }
}


data.basemap_pm25_region <- function(region, year=2020){

  if(region != "Global") stop("Now we only use Global from now on")
  f_pm25_nc <- glue("V5GL04.HybridPM25.Global.{year}01-{year}12.nc")
  f_pm25_tif <- gsub("\\.nc","\\.tif", f_pm25_nc)

  nc_to_tif <- function(f_pm25_nc, f_pm25_tif){
    pm25.base <- raster::raster(creahelpers::get_concentration_path(f_pm25_nc))

    if(region=="Global" & year==2018){
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


data.basemap_no2 <- function(pop, res, use_cache=T, suffix=""){

  f <- sprintf("cache/no2_ugm3_%s%s.tif", res, suffix)
  if(file.exists(f) && use_cache){
    terra::rast(f)
  }else{
    f_no2 <- "Global_LUR_NO2_2011_16b_2.tif"
    f_no2_agg10 <- "Global_LUR_NO2_2011_16b_2_agg10.tif"
    f_no2_wsg84 <- "no2_agg10_ugm3_wsg84.tif"

    if(!file.exists(creahelpers::get_concentration_path(f_no2_wsg84))){
      no2 <- terra::rast(creahelpers::get_concentration_path(f_no2))
      terra::NAflag(no2) <- 128
      no2_agg10 <- terra::aggregate(no2, fact=10, cores=parallel::detectCores()-1)
      terra::writeRaster(no2_agg10, creahelpers::get_concentration_path(f_no2_agg10),
                         overwrite=T)
      no2_wsg84 <- no2_agg10 %>% terra::project("epsg:4326") # Required for hia
      no2_wsg84 <- no2_wsg84 * 1.88 # ppb to µg/m3 !!!
      terra::writeRaster(no2_wsg84, creahelpers::get_concentration_path(f_no2_wsg84),
                         overwrite=T)
    }

    no2 <- terra::rast(creahelpers::get_concentration_path(f_no2_wsg84)) %>%
      terra::resample(pop) %>% terra::mask(pop)
    terra::writeRaster(no2, f, overwrite=T)
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

data.bbox <- function(){
  data.gadm0() %>%
    subset(GID_0 %in% c('PHL')) %>% extent %>% multiply_by(1.1)
}

data.pop <- function(res="30_sec", bbox=NULL, mask=NULL){
  r <- terra::rast(creahelpers::get_population_path(
    paste0('gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2020_',res,'.tif')))

  if(!is.null(bbox)){
    r <- r %>% terra::crop(bbox)
  }

  if(!is.null(mask)){
    r <- terra::mask(r, terra::vect(mask))
  }
  names(r) <- "gpw"
  r
}

data.pop_ratio_log <- function(pop, res, use_cache=T, suffix=""){

  f <- file.path("cache", paste0("pop_ratio_log_",res,suffix,".tif"))
  if(!use_cache | !file.exists(f)){

    r2010 <- terra::rast(creahelpers::get_population_path(
      sprintf('gpw_v4_population_density_rev11_2010_%s.tif',res))) %>%
      terra::resample(pop) %>% terra::mask(pop)

    r2020 <- terra::rast(creahelpers::get_population_path(
      sprintf('gpw_v4_population_density_rev11_2020_%s.tif',res))) %>%
      terra::resample(pop) %>% terra::mask(pop)

    r_ratio <- max(log(min(r2020/r2010,10)), -2)
    names(r_ratio) <- "pop_ratio"
    writeRaster(r_ratio, f, overwrite=T)
    return(r_ratio)
  }else{
    rast(f)
  }
}






data.grump <- function(pop, res, use_cache=T, suffix=""){

  f <- file.path("cache", paste0("grump_",res,suffix,".tif"))

  if(!use_cache | !file.exists(f)){
    # https://sedac.ciesin.columbia.edu/data/collection/grump-v1
    grump <- terra::rast(creahelpers::get_population_path('glurextents.bil'))
    n_aggregate <- min(floor(res(pop) / res(grump)))
    if(n_aggregate>1){
      # We aggregate first to avoid memory issues
      grump <- grump %>% terra::aggregate(fact=n_aggregate, fun="modal")
    }
    grump <- grump %>% terra::resample(pop, method="near")
    terra::writeRaster(grump, f, overwrite=T)
    return(grump)
  }else{
    terra::rast(f)
  }
}

data.lon <- function(pop, res, use_cache=T, suffix=""){
  f <- file.path("cache", paste0("lon_",res,suffix,".tif"))
  if(!use_cache | !file.exists(f)){
    xy <- coordinates(raster(pop))
    r_lon <- raster(pop)
    r_lon[] <- xy[,1]
    r_lon <- r_lon %>% raster::mask(raster(pop))
    terra::writeRaster(r_lon, f, overwrite=T)
    return(r_lon)
  }else{
    terra::rast(f)
  }
}

data.lat <- function(pop, res, use_cache=T, suffix=""){
  f <- file.path("cache", paste0("lat_",res,suffix,".tif"))
  if(!use_cache | !file.exists(f)){
    xy <- coordinates(raster(pop))
    r_lat <- raster(pop)
    r_lat[] <- xy[,2]
    r_lat <- r_lat %>% raster::mask(raster(pop))
    terra::writeRaster(r_lat, f, overwrite=T)
    return(r_lat)
  }else{
    terra::rast(f)
  }
}

data.pm25_merra2_diff <- function(pop, res, year_i=2020, year_f=2020, use_cache=T, suffix=""){

  f <- file.path("cache", sprintf("pm25_merra2_diff_%s_%s_%s%s.tif", year_i, year_f, res, suffix))

  if(!use_cache | !file.exists(f)){

    pm25_merra2_i <- terra::rast(creahelpers::get_concentration_path(sprintf("pm25_merra2_%s.tif", year_i))) %>%
      terra::resample(pop) %>% terra::mask(pop)

    pm25_merra2_f <- terra::rast(creahelpers::get_concentration_path(sprintf("pm25_merra2_%s.tif", year_f))) %>%
      terra::resample(pop) %>% terra::mask(pop)

    pm25_merra2_diff <- pm25_merra2_f - pm25_merra2_i
    pm25_merra2_diff <- pm25_merra2_diff * 1E10 # Scaling a bit...
    terra::writeRaster(pm25_merra2_diff, f, overwrite=T)
    return(pm25_merra2_diff)
  }else{
    terra::rast(f)
  }
}

data.no2_omi_diff <- function(pop, res, year_i=2011, year_f=2019, use_cache=T, suffix=""){

  f <- file.path("cache", sprintf("no2_omi_diff_%s_%s_%s%s.tif",year_i, year_f, res, suffix))

  if(!use_cache | !file.exists(f)){

    omi_i <- terra::rast(creahelpers::get_concentration_path(sprintf("no2_omi_%s.tif", year_i))) %>%
      terra::resample(pop) %>% terra::mask(pop)

    omi_f <- terra::rast(creahelpers::get_concentration_path(sprintf("no2_omi_%s.tif", year_f))) %>%
      terra::resample(pop) %>% terra::mask(pop)

    no2_omi_diff <- omi_f - omi_i
    no2_omi_diff <- no2_omi_diff / 1E15 # Scaling a bit...
    terra::writeRaster(no2_omi_diff, f, overwrite=T)
    return(no2_omi_diff)
  }else{
    terra::rast(f)
  }
}


data.srtm <- function(pop, res, use_cache=T, suffix=""){

  f <- file.path("cache", paste0("srtm_",res,suffix,".tif"))

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
data.srtm_mean <- function(pop, res, d_deg=1, use_cache=T, suffix=""){

  name <- sprintf("srtm_mean_%sdeg", d_deg)
  f <- file.path("cache", sprintf("srtm_mean_%sdeg_%s%s.tif", d_deg, res, suffix))

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



#' Download roads from SEDAC
#' https://sedac.ciesin.columbia.edu/data/set/groads-global-roads-open-access-v1/data-download
#'
#' @param use_cache
#'
#' @return
#' @export
#'
#' @examples
data.road_density_groads <- function(res, pop, use_cache=T, suffix=""){

  f <- file.path("cache", paste0("road_density_groads_",res,suffix,".tif"))

  if(!use_cache | !file.exists(f)){

    res_rasterizing <- "2pt5_min"
    pop_rasterizing <- data.pop(res=res_rasterizing)

    # Using this approach
    # https://gis.stackexchange.com/questions/119993/convert-line-shapefile-to-raster-value-total-length-of-lines-within-cell
    library(rgdal)
    library(raster)
    library(rgeos)
    library(tictoc)
    roads <- readOGR(file.path(creahelpers::get_gis_dir(), "roads",
                               "gROADS_v1.gdb"))

    # roads_utm <- spTransform(roads, CRS("+init=epsg:21037"))
    grid <- raster(raster(pop))
    icells <- which(!is.na(values(raster(pop))))

    tic()
    lengths <- sapply(icells[1:100], function(i) {

      if(is.na(pop[i])){
        return(NA)
      }

      tmp_rst <- grid
      tmp_rst[i] <- 1
      tmp_shp <- rasterToPolygons(tmp_rst)

      if (gIntersects(roads, tmp_shp)) {
        roads_crp <- crop(roads, tmp_shp)
        roads_crp_length <- gLength(roads_crp)
        return(roads_crp_length)
      } else {
        return(0)
      }
    })
    toc()

    r <- raster::rasteri

    # roads <- sf::read_sf(file.path(creahelpers::get_gis_dir(), "roads",
    # "gROADS_v1.gdb"))
    # roads <- roads %>% as("Spatial")

    # r <- creahelpers::rasterize_lines(lines=roads,
    #                                   grid=raster(pop_rasterizing))

    lines <- roads
    grid <- pop
    raster::crs(lines) <- raster::crs(raster::raster(pop))

    # Cut lines along grid cells
    print("Polygonizing...")
    rs <- grid
    rs[] <- 1:ncell(rs)
    names(rs) <- "i_cell"
    rsp <- terra::as.polygons(rs)
    rsp <- as(rsp, "Spatial")
    print("Done")


    # Add temporary feature id for grouping
    lines$feature_id_tmp <- 1:nrow(lines)

    print("Cutting along grid...")
    # sf much less memory intensive than raster::intersect
    # and faster

    # Chunking it to avoid rgeos_binpredfunc_prepared: maximum returned dense matrix size exceeded
    cutting_successful <- F
    chunk_size <- 1E10
    while(!cutting_successful){
      tryCatch({
        rsp$chunk <- rsp$i_cell %/% chunk_size
        emission.sf <- sf::st_as_sf(lines)
        rp <- pbapply::pblapply(split(sf::st_as_sf(rsp), rsp$chunk),
                                function(rsp_chunk){
                                  sf::st_intersection(emission.sf,rsp_chunk)
                                }) %>%
          do.call("bind_rows",.)
        cutting_successful <- T
      }, error=function(e){
        if("size exceeded" %in% as.character(e)){
          chunk_size <- chunk_size / 100
          warning("Cutting failed: ", e, "\n Trying with smaller chunk size", )
        }else{
          stop(e)
        }
      })
    }

    print("Done")

    print("Calculating length...")
    rp$length <- sf::st_length(rp, byid=TRUE)
    print("Done")

    # # Weighting accordingly
    # print("Weighting by length...")
    # rp <- rp %>%
    #   group_by(feature_id_tmp) %>%
    #   do(mutate(., emission=.$emission * length / sum(.$length)))
    # print("Done")

    print("Rasterizing...")
    rp.sum <- rp %>%
      group_by(i_cell=as.integer(i_cell)) %>%
      summarise(length=sum(length, na.rm=T))

    # Print into raster directly!
    cells_x <- rep(0,ncell(rs))
    cells_x[rp.sum$i_cell] <- rp.sum$length
    grid_result <- grid
    grid_result[] <- cells_x
    writeRaster(r, f)
    return(r)
  }else{
    raster(f)
  }


}



#' Download roads from GRIP4
#'
#' @param use_cache
#'
#' @return
#' @export
#'
#' @examples
data.road_density_grip <- function(res, pop, use_cache=T, suffix=""){

  f <- file.path("cache", paste0("road_density_grip_",res,suffix,".tif"))

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


data.distance_urban <- function(pop, res, use_cache=T, suffix=""){

  f <- file.path("cache", paste0("distance_urban_",res, suffix,".tif"))

  if(!use_cache | !file.exists(f)){

    f_gis <- creahelpers::get_landcover_path(glue("grumpv1/distance_urban_{res}.tif"))
    if(!file.exists(f_gis)){
      # Generate a global one first
      # Very slow though
      warning(glue("Generating {basename(f_gis)}. It might take quite some time..."))

      pop_all <- data.pop(res=res)
      # distance doesn't work at 30sec res, we do as best as we can....
      pop_2pt5_min <- data.pop(res=RES_2PT5_MIN)
      grump <- terra::rast(creahelpers::get_population_path('glurextents.bil'))
      # Resample using max to keep urban areas
      grump <- grump %>% terra::resample(pop_2pt5_min, method="max")

      # Compute distance
      sea_level <- 0
      rural_level <- 1
      urban_level <- 2
      grump[is.na(pop_2pt5_min)] <- NA

      dir.create(dirname(f_gis), showWarnings = F, recursive = T)
      dist <- terra::gridDist(grump, target=urban_level)

      dist_res <- dist %>% terra::resample(pop_all, method='bilinear')
      terra::writeRaster(dist_res, f_gis, overwrite=T)
    }

    dist <- terra::rast(f_gis) %>%
      terra::resample(pop) %>%
      terra::mask(pop)

    raster::writeRaster(dist, f, overwrite=T)
  }else{
    dist <- terra::rast(f)
  }
  return(dist)
}


data.gadm_raster <- function(pop, res, level, use_cache=T, suffix="", as_factor=T){

  f <- file.path("cache", paste0(sprintf("gadm%d_%s%s.tif", level, res, suffix)))

  gadm_raster <- if(!use_cache | !file.exists(f)){
    g <- creahelpers::get_adm(level, "full")
    gid_level <- sprintf("GID_%d",level)

    # gid <- sf::st_as_sf(g[gid_level]) %>%
    #   rename_at(gid_level, function(x)"gid") %>%
    #   mutate(gid=as.numeric(factor(gid)))

    r <- terra::rasterize(terra::vect(g),
                          pop,
                          field=gid_level,
                          fun="min")

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

    r <- creahelpers::to_raster(r)
    # We want 0.3 deg min. But if res is coarser, let's increase the focal
    w <- raster::focalWeight(r, max(min(res(r)) * 1.1, 0.3), "rectangle")
    w[] <- 1
    r_filled <- creahelpers::focal.loop(r, w, fill.na, NAonly=T)
    # values(r_filled)[is.na(values(raster(pop)))] <- NA
    raster::writeRaster(r_filled, f, overwrite=T)
    return(r_filled)
  }else{
    terra::rast(f)
  }

  if(as_factor){
    gadm_raster <- gadm_raster %>% terra::as.factor()
  }

  return(gadm_raster)
}

data.landuse <- function(pop, res, use_cache=T, suffix=""){

  f <- file.path("cache", paste0("landuse_",res,suffix,".tif"))

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


data.distance_coast <- function(pop, res, use_cache=T, suffix=""){

  f <- file.path("cache", paste0("distance_coast_",res,suffix,".tif"))
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
data.pm25_ss_dust_frac <- function(pop, res, use_cache=T, suffix=""){

  f <- sprintf("cache/pm25_ss_dust_frac_%s%s.tif", res, suffix)
  if(file.exists(f) && use_cache){
    terra::rast(f)
  }else{
    pm25_w_ssd <- terra::rast(creahelpers::get_concentration_path("ACAG_PM25_GWR_V4GL03_201901_201912_0p01.tif"))
    pm25_wo_ssd <- terra::rast(creahelpers::get_concentration_path("ACAG_PM25_noDUSTnoSEASALT_GWR_V4GL03_201901_201912_0p01.tif"))
    ss_dust_frac <- (pm25_w_ssd - pm25_wo_ssd)/ pm25_w_ssd
    ss_dust_frac <- ss_dust_frac %>% terra::resample(terra::rast(pop))
    terra::writeRaster(ss_dust_frac, f, overwrite=T)
    return(ss_dust_frac)
  }
}

#' Seasalt and dust contribution
#'
#' @param pop
#' @param res
#' @param use_cache
#'
#' @return
#' @export
#'
#' @examples
data.pm25_ss_dust <- function(pop, res, use_cache=T, suffix=""){

  f <- sprintf("cache/pm25_ss_dust_%s%s.tif", res, suffix)
  if(file.exists(f) && use_cache){
    terra::rast(f)
  }else{
    pm25_w_ssd <- terra::rast(creahelpers::get_concentration_path("ACAG_PM25_GWR_V4GL03_201901_201912_0p01.tif"))
    pm25_wo_ssd <- terra::rast(creahelpers::get_concentration_path("ACAG_PM25_noDUSTnoSEASALT_GWR_V4GL03_201901_201912_0p01.tif"))
    ss_dust <- pm25_w_ssd - pm25_wo_ssd
    ss_dust <- ss_dust %>% terra::resample(terra::rast(pop))
    terra::writeRaster(ss_dust, f, overwrite=T)
    return(ss_dust)
  }
}


data.type <- function(pop){
  # Tweak for Philippines
  # We added EMB measurements, plus Lauri's estimations using various methods
  # See https://github.com/energyandcleanair/202108_hia_philippines
  # We add a "measured"

  # Note: deprecated, plus doesn't it mean that non "measured" measurements will get ignored?
  type_raster <- pop %>% `names<-`("type")
  type_raster[!is.na(type_raster)] <- 1
  type_raster <- terra::as.factor(type_raster)
  levels(type_raster) <- terra::levels(type_raster)[[1]] %>% mutate(type="measured")
  return(type_raster)
}
