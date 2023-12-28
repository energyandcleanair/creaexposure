adjust_global <- function(res,
                          pop=NULL,
                          poll = c("pm25", "no2"),
                          suffix = "",
                          selected_regions = NULL,
                          year = 2022,
                          use_cache = T,
                          obs = NULL,
                          models = MODEL_GAM) {
  ###################################################
  # Prepare data
  ###################################################
  # Get observations and predictors and these points
  obs <- default_if_null(obs, data.get_obs(year = year, use_cache = use_cache))
  pop <- default_if_null(pop, data.pop(res=res))
  predictors <- data.predictors(pop, res, year = year, suffix = suffix, use_cache = use_cache)
  obs <- utils.add_predictors(obs, predictors)

  regions <- get_regions(selected_regions = selected_regions)
  obs <- add_region_to_obs(obs, regions)

  # Visual check if some important countries haven't been covered
  # obs %>%  filter(poll=="pm25") %>%
  #   group_by(country, region) %>% summarise(count=n()) %>%
  #   mutate(name=countrycode(country, "iso2c", "country.name")) %>%
  #   arrange(desc(count)) %>%
  #   View()


  for(poll in polls){
    for(model in models){
      models.predict(model=models,
                     poll=poll,
                     regions=regions,
                     res=res,
                     year=year,
                     suffix=suffix,
                     predictors=predictors,
                     obs=obs)
    }
  }
}
