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

utils.add_predictors <- function(obs, predictors){
  cbind(obs, raster::extract(predictors, sf::st_as_sf(obs)))
}




#' Rasterize lines in parallel
#'
#' @param lines
#' @param grid
#'
#' @return
#' @export
#'
#' @examples
utils.rasterize_lines <- function(lines, grid){

  library(pbapply)
  library(SpaDES)
  library(parallel)

  grid <- utils.to_raster(grid) %>% raster()
  lines <- sf::st_as_sf(lines)

  # Add temporary feature id for later grouping
  lines$feature_id_tmp <- 1:nrow(lines)


  n <- 10 # Might need to increased this one if proc is killed
  grids <- splitRaster(grid, nx=n, ny=n)
  cl <- makeForkCluster(parallel::detectCores()-1)

  ress <- pblapply(grids, function(g){

      g[] <- 1:ncell(g)
      names(g) <- "i_cell"

      message("Polygonizing...")
      gsp <- as(g,'SpatialPolygonsDataFrame')
      gsf <- sf::st_as_sf(gsp)

      message("Cutting along grid...")
      g_lines <- suppressMessages(sf::st_intersection(lines, gsf))

      message("Calculating length...")
      g_lines$length <- suppressMessages(sf::st_length(g_lines, byid=TRUE))
      print("Done")

      message("Rasterizing...")
      g_length <- g_lines %>%
        group_by(i_cell=as.integer(i_cell)) %>%
        summarise(length=sum(length, na.rm=T))

      # Print into raster directly!
      cells_x <- rep(0, ncell(g))
      cells_x[g_length$i_cell] <- g_length$length
      g[] <- cells_x
      return(g)
    },
    cl=cl)
  stopCluster(cl)
  res <- do.call(merge, ress)
  return(res)
}
