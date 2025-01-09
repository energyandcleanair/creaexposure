#' Mask areas far from observation stations
#'
#' @param delta_raster Raster or SpatRaster to mask
#' @param obs sf object containing observation stations
#' @param outer_buffer_km Maximum distance from stations (areas beyond this will be masked out)
#' @param smooth Logical. If TRUE, creates a smooth transition zone
#' @param inner_buffer_km Distance from stations where values remain unchanged. Only used if smooth=TRUE.
#'                       Areas between inner and outer buffer will have gradual transition.
#' @param diagnostics_folder Path to save diagnostic plots. Set NULL to skip plotting
#' @param suffix String to append to diagnostic plot filename
#'
#' @return Masked Raster in original CRS
mask_far_from_stations <- function(delta_raster,
                                   obs,
                                   outer_buffer_km,
                                   smooth,
                                   inner_buffer_km,
                                   diagnostics_folder = "diagnostics",
                                   suffix = "") {
  r <- creahelpers::to_rast(delta_raster)
  # Convert buffers to meters
  outer_buffer_m <- outer_buffer_km * 1000
  inner_buffer_m <- if (!is.null(inner_buffer_km)) inner_buffer_km * 1000 else outer_buffer_m * 0.8

  # Create distance raster in Web Mercator for accurate distances
  obs_3857 <- terra::vect(obs %>% filter(!is.na(lat)), crs = "EPSG:4326") %>% terra::project("EPSG:3857")
  temp_r <- r
  temp_r_3857 <- terra::project(temp_r, "EPSG:3857")
  stations_raster <- terra::rasterize(obs_3857, temp_r_3857)
  dist_from_stations <- terra::distance(stations_raster)

  # Project distance raster back to original CRS
  dist_4326 <- terra::project(dist_from_stations, "EPSG:4326")

  if (!smooth) {
    # Simple binary mask
    w <- dist_4326 <= outer_buffer_m
    w[w == 0] <- 0
  } else {
    # Smooth transition
    w <- terra::rast(dist_4326, nlyrs = 1, names = "weight", vals = 0)

    # Assign weights based on distance thresholds
    w[dist_4326 <= inner_buffer_m] <- 1
    w[dist_4326 >= outer_buffer_m] <- 0
    transition_zone <- dist_4326 > inner_buffer_m & dist_4326 < outer_buffer_m
    w[transition_zone] <- (outer_buffer_m - dist_4326[transition_zone]) /
      (outer_buffer_m - inner_buffer_m)
  }

  # Diagnose mask if diagnostics folder is provided
  if (!is.null(diagnostics_folder)) {
    filepath <- file.path(
      diagnostics_folder,
      glue("mask_far_from_stations{suffix}.jpg")
    )

    # Set up multi-panel plot
    jpeg(filepath, width = 12, height = 6, units = "in", res = 300)
    par(mfrow = c(1, 2))

    # Plot 1: Distance from stations with station locations
    raster::plot(dist_4326 / 1000,
      main = "Distance from stations (km)",
      col = viridis::viridis(100)
    )
    plot(sf::st_geometry(sf::st_as_sf(obs)), add = TRUE, pch = 16, col = "red", cex = 0.4)

    # Plot 2: Resulting mask weights
    if (smooth) {
      raster::plot(w,
        main = glue("Mask weights\n(inner: {round(inner_buffer_km)}km, outer: {round(outer_buffer_km)}km)"),
        col = viridis::viridis(100)
      )
    } else {
      raster::plot(w,
        main = glue("Binary mask\n(threshold: {round(outer_buffer_km)}km)"),
        col = c("white", "grey60")
      )
    }
    plot(sf::st_geometry(sf::st_as_sf(obs)), add = TRUE, pch = 16, col = "red", cex = 0.4)

    dev.off()
  }

  return(creahelpers::to_raster(w))
}

#' Stations are predominantly in cities. We therefore do not update map in rural areas where
#' the model predictive power might be too weak
#'
#' @param delta_raster Raster to process
#' @param predictors Predictor variables
#' @param obs Observation data
#' @param quantile Quantile threshold for distance to urban areas
#'
#' @return Weight raster (1 = keep, 0 = mask)
#' @export
mask_far_from_urban <- function(delta_raster, predictors, obs, quantile = 0.95) {
  # Create weight raster (1 where we keep values, 0 where we mask)
  mask <- predictors$distance_urban <= quantile(obs$distance_urban, quantile, na.rm = TRUE)
  mask[mask == 0] <- 0
  return(mask)
}


#' Same as mask_far_from_urban but with a smooth transition zone
#'
#' @param delta_raster Raster to process
#' @param predictors Predictor variables
#' @param obs Observation data
#' @param quantile_inner Inner quantile threshold for smooth transition
#' @param quantile_outer Outer quantile threshold for smooth transition
#' @param diagnostics_folder Path to save diagnostic plots
#' @param suffix String to append to diagnostic plot filename
#'
#' @return Weight raster (1 = keep, 0 = mask, with smooth transition)
#' @export
mask_far_from_urban_smooth <- function(delta_raster,
                                       predictors,
                                       obs,
                                       quantile_inner = 0.90,
                                       quantile_outer = 0.95,
                                       diagnostics_folder = "diagnostics",
                                       suffix = "") {
  # Calculate inner and outer distance thresholds based on quantiles
  inner_distance <- quantile(obs$distance_urban, quantile_inner, na.rm = TRUE)
  outer_distance <- quantile(obs$distance_urban, quantile_outer, na.rm = TRUE)

  # Get the distance to urban areas raster
  dist <- predictors$distance_urban

  # Initialize the weighting raster
  w <- raster::raster(dist)

  # Assign weights based on distance thresholds
  w[dist <= inner_distance] <- 1
  w[dist >= outer_distance] <- 0
  transition_zone <- dist > inner_distance & dist < outer_distance
  w[transition_zone] <- (outer_distance - dist[transition_zone]) / (outer_distance - inner_distance)

  # Diagnose mask
  if (!is.null(diagnostics_folder)) {
    filepath <- file.path(diagnostics_folder, glue("mask_far_from_urban{suffix}.jpg"))
    jpeg(filepath, width = 8, height = 8, units = "in", res = 300)
    raster::plot(w)
    dev.off()
  }

  return(creahelpers::to_raster(w))
}
