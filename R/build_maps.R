#' Build the ex
#'
#' @param res
#' @param pop
#' @param polls
#' @param suffix
#' @param regions
#' @param year
#' @param force_rebuild
#' @param use_cache_predictors
#' @param obs
#' @param obs_level
#' @param model
#' @param results_folder
#' @param diagnostics_folder
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
build_maps <- function(res,
                      pop=NULL,
                      polls = c("pm25", "no2"),
                      suffix = "",
                      regions = NULL,
                      year = 2022,
                      force_rebuild = T,
                      use_cache_predictors = T,
                      obs = NULL,
                      obs_level = "city",
                      model = MODEL_GAM,
                      results_folder = file.path("results", model),
                      diagnostics_folder = "diagnostics",
                      formula=NULL,

                      # Masking options
                      limit_distance_urban = F,
                      distance_urban_quantile_inner = 0.90,
                      distance_urban_quantile_outer = 0.95,

                      limit_distance_stations = F,
                      distance_stations_inner_km = 50,
                      distance_stations_outer_km = 100,
                      distance_stations_smooth = T,
                      ...) {

  logger::log_layout(creahelpers::log_layout_crea)

  ###################################################
  # Prepare data
  ###################################################
  # Get observations and predictors and these points
  obs <- default_if_null(obs, data.get_obs(level = obs_level, year = year, use_cache = T))
  bbox <- default_if_null(bbox, get_bbox(regions))
  pop <- default_if_null(pop, data.pop(res=res, bbox=bbox))

  # Given that predictors only rely on selected regions (for now), we use different suffix than for results
  suffix_predictors <- ifelse(is.null(regions), "", paste(c("_",tolower(get_region_names(regions))), collapse=""))
  predictors <- data.predictors(pop, res, year=year, suffix=suffix_predictors, use_cache=use_cache_predictors, model=model)
  obs_w_predictors <- utils.add_predictors(obs, predictors)
  obs_w_predictors <- add_region_to_obs(obs_w_predictors, regions)

  # Diagnose predictors (well mainly plot them)
  diagnose_predictors(predictors = predictors,
                      suffix = suffix,
                      formula = formula,
                      folder = diagnostics_folder)

  maps <- list()
  for(poll in polls){
    models.predict(model=model,
                   poll=poll,
                   regions=regions,
                   res=res,
                   year=year,
                   suffix=suffix,
                   predictors=predictors,
                   obs=obs_w_predictors,
                   results_folder=results_folder,
                   diagnostics_folder = diagnostics_folder,
                   force_rebuild=force_rebuild,
                   formula=formula,
                   # Masking
                   limit_distance_urban=limit_distance_urban,
                   distance_urban_quantile_inner=distance_urban_quantile_inner,
                   distance_urban_quantile_outer=distance_urban_quantile_outer,
                   limit_distance_stations=limit_distance_stations,
                   distance_stations_inner_km=distance_stations_inner_km,
                   distance_stations_outer_km=distance_stations_outer_km,
                   distance_stations_smooth=distance_stations_smooth,
                   ...) -> maps[[poll]]
  }

  # Run diagnostics on results
  diagnostics <- diagnose_results(
    maps = maps,
    polls = polls,
    obs_w_predictors = obs_w_predictors,
    suffix = suffix,
    folder = diagnostics_folder)


  # Postcompute e.g. plot maps
  units <- obs %>% distinct(poll, unit) %>% deframe()
  maps_plots <- plot_maps(maps=maps, units=units)

  return(list(
    maps = maps,
    diagnostics = diagnostics,
    maps_plots = maps_plots
  ))
}




build_map <- function(...) {
  # Use lifecycle's deprecate_warn to show a warning message
  deprecate_warn("0.4.0", "build_map()", "build_maps()")

  # Call the new function to maintain backward compatibility
  build_maps(...)
}
