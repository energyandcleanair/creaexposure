build_map <- function(res,
                      pop=NULL,
                      polls = c("pm25", "no2"),
                      suffix = "",
                      selected_regions = NULL,
                      year = 2022,
                      force_rebuild = T,
                      obs = NULL,
                      obs_level = "city",
                      model = MODEL_GAM,
                      results_folder = file.path("results", model),
                      ...) {

  logger::log_layout(creahelpers::log_layout_crea)

  ###################################################
  # Prepare data
  ###################################################
  # Get observations and predictors and these points
  obs <- default_if_null(obs, data.get_obs(level = obs_level, year = year, use_cache = T))
  pop <- default_if_null(pop, data.pop(res=res, bb=get_bb(selected_regions)))
  # Given that predictors only rely on selected regions (for now), we use different suffix than for results
  suffix_predictors <- ifelse(is.null(selected_regions), "", paste(c("_",tolower(selected_regions)), collapse=""))
  predictors <- data.predictors(pop, res, year = year, suffix = suffix_predictors, use_cache = T)
  obs <- utils.add_predictors(obs, predictors)

  regions <- get_regions(selected_regions = selected_regions)
  obs <- add_region_to_obs(obs, regions)

  # Visual check if some important countries haven't been covered
  # obs %>%  filter(poll=="pm25") %>%
  #   group_by(country, region) %>% summarise(count=n()) %>%
  #   mutate(name=countrycode(country, "iso2c", "country.name")) %>%
  #   arrange(desc(count)) %>%
  #   View()


  result <- list()
  for(poll in polls){
    models.predict(model=model,
                   poll=poll,
                   regions=regions,
                   res=res,
                   year=year,
                   suffix=suffix,
                   predictors=predictors,
                   obs=obs,
                   results_folder=results_folder,
                   force_rebuild=force_rebuild,
                   ...) -> result[[poll]]
  }

  return(result)
}
