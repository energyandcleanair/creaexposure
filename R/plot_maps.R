#' Plotting maps in beautiful images
#'
#' @param maps list of maps, with names indicating poll
#' @param units named list with poll=unit
#'
#' @return list of ggplots
#' @export
plot_maps <- function(maps, units){

  plots <- list()
  for(poll in names(maps)){
    map <- maps[[poll]]
    unit <- units[[poll]]
    plots[[poll]] <- plot_map(map, poll, unit)
  }

  return(plots)
}


plot_map <- function(map,
                     poll,
                     unit,
                     tile_type = "cartolight",
                     scale_fill = scale_fill_viridis(option = "rocket", direction=-1)
                     ){

  library(raster)
  library(ggplot2)
  library(ggspatial)
  library(sf)
  library(prettymapr)


  # Make it a terra rast and convert to 3857
  map <- creahelpers::to_rast(map) %>%
    terra::project("EPSG:3857")

  # Get best zoom level from bbox
  zoom <- get_zoom_level(map)

  # Convert the raster to a data frame
  raster_df <- as.data.frame(map, xy = TRUE)
  names(raster_df)[3] <- "value"  # Rename the value column

  # Create the ggplot
  plt <- ggplot() +
    # Add the basemap
    annotation_map_tile(type = tile_type,
                        zoom = zoom) +
    # Add the raster layer
    geom_raster(data = raster_df, aes(x = x, y = y, fill = value), alpha=0.8) +
    # Customize the color scale
    scale_fill +

    # Cut exactly at raster bbox
    coord_sf(xlim = c(terra::ext(map)[1], terra::ext(map)[2]),
             ylim = c(terra::ext(map)[3], terra::ext(map)[4]),
             expand = FALSE) +
    # Add a north arrow and scale bar

    annotation_scale(location = "br", plot_unit="m") +

    # Customize the theme (hide latitude, longitude)
    theme_map() +

    # Add white background to the scale
    # theme(legend.key = element_rect(fill = "white")) +
    # Put legent underneath, hotizontally, using the whold width
    theme(legend.position = "bottom",
          legend.justification = "center",
          legend.direction = "horizontal",
          legend.key.width = unit(1, "cm"),
          legend.key.height = unit(0.5, "cm"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8)) +
    labs(fill = glue("{rcrea::poll_str(poll)} ({unit})")) +
    # make all background transparent
    theme(panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", colour = NA))

  return(plt)
}



#' Compute best zoom level based on the extent (in EPSG:3857)
#'
#' @param ext
#' @param factor
#'
#' @return
#' @export
#'
#' @examples
get_zoom_level <- function(map, delta=-1) {

  ext <- terra::ext(map)

  # Get the extent
  x_range <- ext[2] - ext[1]
  y_range <- ext[4] - ext[3]

  # Get the zoom level
  zoom <- round(21 - log2(max(x_range, y_range) / 256) + delta)
  return(zoom)
}
