get_grid <- function(res, use_terra=T, bbox=NULL){
  pop <- data.pop(res=res, bbox=bbox)
  if(use_terra){
    return(terra::rast(creahelpers::to_rast(pop)))
  }else{
    return(raster::raster(creahelpers::to_raster(pop)))
  }
}
