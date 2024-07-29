get_grid <- function(res, use_terra=T, bb=NULL){
  pop <- data.pop(res=res, bb=bb)
  if(use_terra){
    return(terra::rast(creahelpers::to_rast(pop)))
  }else{
    return(raster::raster(creahelpers::to_raster(pop)))
  }
}
