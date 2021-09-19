utils.focal_mean <- function(r, d_deg, pop, res, use_cache=T){

  name <- sprintf("%s_%sdeg", names(r), d_deg)
  f <- file.path("cache", sprintf("%s_%s.tif", name, res))

  if(!use_cache | !file.exists(f)){

    w <- raster::focalWeight(raster(srtm), d_deg, "circle")
    w[w>0] <- 1
    r.focal <- raster::focal(raster(srtm),
                                     w,
                                     na.rm=T,
                                     pad=T,
                                     fun=mean) %>%
      terra::rast()

    r.focal[is.na(pop)] <- NA
    names(r.focal) <- name
    terra::writeRaster(r.focal, f, overwrite=T)
    return(r.focal)
  }else{
    terra::rast(f)
  }
}

utils.to_raster <- function(x){
  if(class(x)[1]!="RasterLayer"){
    raster(x)
  }else{
    x
  }
}
