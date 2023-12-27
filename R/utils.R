utils.focal_mean <- function(r, d_deg, pop, res, use_cache=T, suffix=""){

  name <- sprintf("%s_%sdeg", names(r), d_deg)
  f <- file.path("cache", sprintf("%s_%s%s.tif", name, res, suffix))

  if(!use_cache | !file.exists(f)){

    w <- raster::focalWeight(raster(r), d_deg, "circle")
    w[w>0] <- 1
    r.focal <- raster::focal(raster(r),
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

#' Remove constant variables from formula to avoid modeling errors
#'
#' @param formula
#' @param data
#'
#' @return
#' @export
#'
#' @examples
utils.remove_constant_variables_from_formula <- function(formula, data){

  # At the moment, only works for numeric values
  sd_safe <- function(x){
    if(!is.numeric(x)) return(NA)
    sd(x)
  }
  constant_variables <- names(which(sapply(data, sd_safe)==0))

  remove_vars_from_formula <- function(formula, vars){
    all_vars <- strsplit(as.character(formula[3]), " \\+ ")[[1]]
    for(x in vars){
      var_idx <- which(grepl(x, all_vars))
      if(length(var_idx)==0) next
      print(glue("Removing {x} from formula - constant"))
      all_vars <- all_vars[-var_idx]
    }

    formula <- as.formula(
      paste(formula[2], "~", paste(all_vars, collapse=" + "))
    )
    formula
  }

  formula <- remove_vars_from_formula(formula, vars=constant_variables)

  return(formula)
}
