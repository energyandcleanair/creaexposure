
diagnose_results <- function(maps, polls, obs_w_predictors, suffix="", folder="diagnostics"){

  diagnostics <- list(
    plot_points = list(),
    map = list(),
    data = list()
  )


  for(poll in polls){

    # Points points -----------------------------------------------------------
    map <- maps[[poll]]
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
      select(location_id, country, starts_with("value")) %>%
      pivot_longer(cols = c(value_prior, value_predicted),
                   names_to = "variable",
                   names_prefix = "value_",
                   values_to = "value") %>%
      ggplot(aes(x=value_measured, y=value, color=variable)) +
      geom_point() +
      coord_equal() +
      geom_abline(intercept = 0, slope = 1) +
      ggrepel::geom_text_repel(aes(label=location_id), size=2) -> plt
      facet_wrap(~country) -> plt


    diagnostics[["plot_points"]][[poll]] <- plt
    diagnostics[["data"]][[poll]] <- data_sf

    ggsave(paste0(folder, "/results_points_", poll, suffix, ".png"))

    # Map view ----------------------------------------------------------------
    library(gstat)
    library(sp)
    library(tidyr)
    library(dplyr)
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

    diagnostics[["map"]][[poll]] <- plt
    ggsave(paste0(folder, "/results_map_w_points", poll, suffix, ".png"))
  }

  return(diagnostics)
}


diagnose_predictors <- function(predictors, formula=NULL, suffix="", folder="diagnostics"){


  if(!is.null(formula)){
    # Extract useful predictors from formula
    layers <- intersect(all.vars(formula), names(predictors))
  }else{
    layers <- names(predictors)
  }

  predictors <- predictors[[layers]]
  n_layers <- raster::nlayers(predictors)
  height <- ceiling(sqrt(n_layers)) * 4
  width <- height

  filepath <- file.path(folder, glue("predictors{suffix}.jpg"))
  jpeg(filepath, width=width, height=height, units="in", res=300)
  #TODO: Check why it only exports 16 layers when we have more
  raster::plot(predictors, col=viridis::viridis(n_layers),
                      layout=c(ceiling(sqrt(n_layers)),
                               ceiling(sqrt(n_layers))))
  dev.off()
}



