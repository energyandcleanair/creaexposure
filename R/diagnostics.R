diagnose_obs_w_predictors <- function(obs_w_predictors, suffix="", folder="diagnostics"){

  dir.create(folder, showWarnings = FALSE, recursive = TRUE)

  d <- obs_w_predictors %>%
    select(location_id, country, poll, value, pm25_prior, no2_prior) %>%
    tidyr::pivot_longer(cols = c("pm25_prior", "no2_prior"),
                        names_transform = list(poll_prior = ~ gsub("_prior", "", .)),
                        names_to = "poll_prior", values_to = "prior") %>%
    filter(poll==poll_prior)

  ggplot(d, aes(x=prior, y=value)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    scale_x_continuous(limits = c(0, max(max(d$value), max(d$prior)))) +
    scale_y_continuous(limits = c(0, max(max(d$value), max(d$prior)))) +
    facet_wrap(poll~country) +
    ggrepel::geom_text_repel(aes(label=location_id), size=2)

  ggsave(paste0(folder, "/obs_w_predictors", suffix, ".png"))
}


diagnose_results <- function(results, polls, obs_w_predictors, suffix="", folder="diagnostics"){

  for(poll in polls){

    # Points points -----------------------------------------------------------
    map <- results[[poll]]
    data_sf <- sf::st_as_sf(obs_w_predictors) %>%
      filter(poll == !!poll) %>%
      select(location_id, country, poll, value_measured=value, pm25_prior, no2_prior) %>%
      tidyr::pivot_longer(cols = c("pm25_prior", "no2_prior"),
                          names_transform = list(poll_prior = ~ gsub("_prior", "", .)),
                          names_to = "poll_prior", values_to = "value_prior") %>%
      filter(poll == poll_prior)


    # Extract values at location, knowing that map is a raster
    data_sf$value_predicted <- raster::extract(map, data_sf)

    data_sf %>%
      select(country, starts_with("value")) %>%
      pivot_longer(cols = c(value_prior, value_predicted),
                   names_to = "variable",
                   names_prefix = "value_",
                   values_to = "value") %>%
      ggplot(aes(x=value_measured, y=value, color=variable)) +
      geom_point() +
      coord_equal() +
      geom_abline(intercept = 0, slope = 1) +
      facet_wrap(~country) -> plt

    ggsave(paste0(folder, "/results_points_", poll, suffix, ".png"))



    # Map view ----------------------------------------------------------------
    library(gstat)
    library(sp)
    library(tidyr)
    library(dplyr)
    library(rgdal)
    library(sf)
    library(raster)
    library(ggplot2)
    library(viridis)
    library(ggthemes)
    library(sf)

    # converting raster to df for ggplot
    map_spdf <- as(map, "SpatialPixelsDataFrame") %>% as.data.frame()
    colnames(map_spdf) <- c("value", "x", "y")

    obs_sf <- sf::st_as_sf(obs_w_predictors) %>%
      filter(poll == !!poll) %>%
      # Crop to the map
      sf::st_crop(map)


    ggplot() +
      geom_tile(data=map_spdf, aes(x=x, y=y, fill=value), alpha=0.8) +
      geom_sf(data=obs_sf, inherit.aes = F, aes(col=value), ) +
      scale_fill_viridis(option="D", limits=c(0, max(map_spdf$value))) +
      scale_color_viridis(option="D", limits=c(0, max(map_spdf$value))) +
      # ggrepel::geom_text_repel(data=obs_sf %>% mutate(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2]),
                               # aes(x=x, y=y, label=glue("{location_id}: {round(value)}"))) +
      # viridis::scale_fill_viridis(option="D") +
      # coord_equal() +
      ggthemes::theme_map() +
      theme(legend.position="bottom") +
      theme(legend.key.width=unit(5, "cm")) -> plt

    ggsave(paste0(folder, "/results_map_", poll, suffix, ".png"))
  }
}





