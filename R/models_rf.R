models.rf.predict <- function(obs,
                              poll,
                              predictors,
                              res,
                              year,
                              regions,
                              suffix,
                              force_rebuild = T,
                              results_folder = "results/rf",
                              diagnostics_folder = "diagnostics",
                              remove_seasalt_dust_contribution = F,
                              formula=NULL,

                              # Masking options
                              limit_distance_urban = F,
                              distance_urban_quantile = 0.95,
                              distance_urban_quantile_inner = 0.90,
                              limit_distance_stations = F,
                              distance_stations_inner_km = 50,
                              distance_stations_outer_km = 100,
                              distance_stations_smooth = T,
                              ...) {

  dir.create(results_folder, recursive = T, showWarnings = F)

  predict_fns <- list(
    "pm25" = models.rf.predict.pm25,
    "no2" = models.rf.predict.no2
  )

  predict_fns[[poll]](
    obs = obs,
    predictors = predictors,
    res = res,
    regions = regions,
    suffix = suffix,
    year = year,
    results_folder = results_folder,
    diagnostics_folder = diagnostics_folder,
    limit_distance_urban = limit_distance_urban,
    distance_urban_quantile = distance_urban_quantile,
    limit_distance_stations = limit_distance_stations,
    distance_stations_inner_km = distance_stations_inner_km,
    distance_stations_outer_km = distance_stations_outer_km,
    distance_stations_smooth = distance_stations_smooth,
    remove_seasalt_dust_contribution = remove_seasalt_dust_contribution,
    force_rebuild = force_rebuild,
    formula = formula
  )
}


models.rf.predict.generic <- function(formula,
                                      obs,
                                      region,
                                      predictors,
                                      prior,
                                      poll,
                                      res,
                                      year,
                                      suffix,
                                      limit_distance_stations,
                                      distance_stations_inner_km,
                                      distance_stations_outer_km,
                                      distance_stations_smooth,
                                      limit_distance_urban,
                                      distance_urban_quantile,
                                      remove_seasalt_dust_contribution,
                                      results_folder,
                                      diagnostics_folder = results_folder,
                                      ...) {
  # Prepare data
  data <- obs %>%
    filter(region == !!region) %>%
    filter(poll == !!poll)

  # Keep only variables in formula
  data <- data %>%
    dplyr::select_at(c("value", all.vars(formula))) %>%
    drop_na()

  # Define control parameters for training
  control <- caret::trainControl(method = "cv", number = 10) # Example: 10-fold cross-validation

  # Train Random Forest Model using caret

  model <- caret::train(formula, data = data, method = "rf", trControl = control)


  # Predict
  #TODO Understand why the raster predict != from vector predict (if they are)
  pred_raster <- raster::predict(
    predictors[[setdiff(all.vars(formula), "value")]],
    model
  )

  raster::plot(pred_raster,colNA="black")

  # Remove sea salt contribution
  if (remove_seasalt_dust_contribution & poll == "pm25") {
    pred_raster <- pred_raster * (1 - predictors$pm25_ss_dust_frac)
  }

  # Use delta predicted - prior so that we can smooth the transition
  delta_raster <- pred_raster - prior

  # Get weights from both masking functions
  w_urban <- 1
  w_stations <- 1

  if (limit_distance_urban) {
    w_urban <- mask_far_from_urban_smooth(
      delta_raster = delta_raster,
      predictors = predictors,
      obs = obs,
      quantile_inner = max(0, distance_urban_quantile - 0.05),
      quantile_outer = min(1, distance_urban_quantile + 0.05)
    )
  }

  if (limit_distance_stations) {
    w_stations <- mask_far_from_stations(
      delta_raster = delta_raster,
      obs = obs,
      inner_buffer_km = distance_stations_inner_km,
      outer_buffer_km = distance_stations_outer_km,
      smooth = distance_stations_smooth
    )
  }

  # Combine weights by taking the minimum
  w_combined <- min(w_urban, w_stations)

  # Apply combined weight to delta
  delta_raster <- delta_raster * w_combined

  # Take prior where pred is NA
  delta_raster[is.na(delta_raster)] <- 0

  # New predicted raster
  pred_raster <- prior + delta_raster


  # Diagnose post-transformations
  models.rf.diagnose(
    model = model,
    diagnostics_folder = diagnostics_folder,
    res = res,
    poll = poll,
    region = region,
    year = year,
    suffix = suffix,
    obs = obs,
    pred_raster = pred_raster
  )


  # Return predictions
  return(pred_raster)
}


models.rf.diagnose <- function(model, diagnostics_folder, res, poll, region, year, suffix, obs, pred_raster) {

  dir.create(diagnostics_folder, recursive = T, showWarnings = F)

  filepath <- models.rf.result_filepath(
    poll = poll,
    res = res,
    year = year,
    suffix = suffix,
    folder = diagnostics_folder,
    extension = "txt"
    )

  sink(file = filepath)
  print(model)
  print(summary(model))
  sink(file = NULL)

  model$finalModel$importance %>%
    as.data.frame() %>%
    arrange(desc(IncNodePurity)) %>%
    # head(20) %>%
    # add rownames as a column
    tibble::rownames_to_column("var") %>%
    mutate(var = fct_reorder(var, IncNodePurity)) %>%
    ggplot() +
    geom_col(aes(y = var, x = IncNodePurity)) +
    labs(title = glue("[DIAGNOSTIC] Predictor importance in Random Forest model for {rcrea::poll_str(poll)}"),
         x = "IncNodePurity",
         y = "Predictor") +
    rcrea::theme_crea_new() -> plt

  plt
  ggsave(filepath %>% gsub("\\.txt", ".png", .), plt, width = 10, height = 10)


  # Predicted vs observed (before last transformations)
  d_predicted <- predict(model, model$trainingData)
  d_observed <- model$trainingData$.outcome
  d_prior <- model$trainingData[,glue("{poll}_prior")]


  get_perf_str <- function(d_predicted, d_observed){
    n_pred <- sum(!is.na(d_predicted))
    get_rsq <- function (x, y) cor(x, y, use="complete.obs") ^ 2
    rsq <- get_rsq(d_observed, d_predicted)
    mae <- Metrics::mae(
      d_observed[!is.na(d_predicted)],
      d_predicted[!is.na(d_predicted)])
      glue("{n_pred} predicted values | R2 = {round(rsq,2)} | MAE = {round(mae,2)}\nCross-validation R2 = {round(model$results[model$results$mtry==model$bestTune$mtry,'Rsquared'],2)}")
  }



  bind_rows(
    tibble(x=d_observed, y=d_predicted, type='predicted'),
    tibble(x=d_observed, y=d_prior, type='prior'),
  ) %>%
  ggplot() +
    geom_point(aes(x=x, y=y, col=type)) +
    geom_abline(slope=1, intercept=0) +
    labs(x='Observed',
         y='Predicted or Prior',
         title='Before post-transformations',
         subtitle=get_perf_str(d_predicted, d_observed))

  ggsave(filepath %>% gsub("\\.txt", "_scatter_before.png", .), width = 10, height = 10)


  # Predicted vs observed (after last transformations)
  # obs_cn <- obs %>%
  #   filter(country=='CN', source=='mee')

  d_predicted <- obs %>%
    sf::st_as_sf() %>%
    raster::extract(pred_raster, .)
  d_observed <- obs$value
  d_prior <- obs[,glue("{poll}_prior")]



  bind_rows(
    tibble(x=d_observed, y=d_predicted, type='predicted'),
    tibble(x=d_observed, y=d_prior, type='prior'),
  ) %>%
    ggplot() +
    geom_point(aes(x=x, y=y, col=type)) +
    geom_abline(slope=1, intercept=0) +
    labs(x='Observed',
         y='Predicted | Prior',
         subtitle=get_perf_str(d_predicted, d_observed)) +
    facet_wrap(~type) +
    coord_equal()

  ggsave(filepath %>% gsub("\\.txt", "_scatter_after.png", .), width = 10, height = 10)

}


models.rf.result_filepath <- function(poll, res, year, suffix, folder, extension = "tif"){
  glue("{folder}/{poll}_{res}_{year}{suffix}.{extension}")
}


models.rf.predict.poll <- function(obs,
                                   predictors,
                                   regions,
                                   res,
                                   year,
                                   suffix,
                                   poll,
                                   formula,
                                   prior,
                                   results_folder,
                                   diagnostics_folder = results_folder,
                                   force_rebuild = T,
                                   limit_distance_urban = T,
                                   distance_urban_quantile = 0.95,
                                   limit_distance_stations = F,
                                   distance_stations_inner_km = 50,
                                   distance_stations_outer_km = 100,
                                   distance_stations_smooth = T,
                                   remove_seasalt_dust_contribution = T,
                                   ...) {


  # Check if file exists
  filepath <- models.rf.result_filepath(
    poll = poll,
    res = res,
    year = year,
    suffix = suffix,
    folder = results_folder
  )

  if(file.exists(filepath) & !force_rebuild){
    logger::log_info(glue("File {filepath} already exists. Returning it."))
    return(raster::raster(filepath))
  }

  pred <- pblapply(names(regions), function(region) {
    models.rf.predict.generic(
            formula = formula,
            obs = obs %>% filter(poll == !!poll),
            region = region,
            predictors = predictors,
            poll = poll,
            res = res,
            prior = prior,
            year = year,
            suffix=suffix,
            results_folder = results_folder,
            diagnostics_folder = diagnostics_folder,
            limit_distance_urban = limit_distance_urban,
            distance_urban_quantile = distance_urban_quantile,
            limit_distance_stations = limit_distance_stations,
            distance_stations_inner_km = distance_stations_inner_km,
            distance_stations_outer_km = distance_stations_outer_km,
            distance_stations_smooth = distance_stations_smooth,
            remove_seasalt_dust_contribution = remove_seasalt_dust_contribution,
          )
  }) %>%
    raster::stack() %>%
    terra::rast() %>%
    terra::app(mean, na.rm = T) %>%
    raster() %>%
    `names<-`("predicted")

  # Save
  logger::log_success(glue("Writing predicted map to {filepath}"))
  writeRaster(pred, filepath, overwrite = T)

  # Write parameters in csv
  params <- list(
    "poll" = poll,
    "res" = res,
    "year" = year,
    "suffix" = suffix,
    "limit_distance_urban" = limit_distance_urban,
    "distance_urban_quantile" = distance_urban_quantile,
    "remove_seasalt_dust_contribution" = remove_seasalt_dust_contribution
  )
  write.csv(params, file = filepath %>% gsub("\\.tif", ".csv", .))

  return(pred)
}


models.rf.predict.pm25 <- function(obs,
                                   predictors,
                                   regions, res, year, suffix,
                                   results_folder,
                                   diagnostics_folder,
                                   limit_distance_urban,
                                   distance_urban_quantile,
                                   limit_distance_stations,
                                   distance_stations_inner_km,
                                   distance_stations_outer_km,
                                   distance_stations_smooth,
                                   remove_seasalt_dust_contribution,
                                   force_rebuild,
                                   formula=NULL,
                                   ...) {


  pm25_formula <- if(is.null(formula)) {
    value ~ pm25_prior + pm25_merra2_diff + gadm1 + distance_urban +
      srtm + srtm_diff05deg + srtm_05deg + grump + distance_coast +
      pop_ratio_log + pop
  } else {
    formula
  }

  obs$gadm1 <- as.factor(obs$gadm1) # Should be done beforehand

  models.rf.predict.poll(
    formula = pm25_formula,
    obs = obs,
    predictors = predictors,
    regions = regions,
    res = res,
    year = year,
    suffix = suffix,
    poll = "pm25",
    prior = predictors$pm25_prior,
    results_folder = results_folder,
    diagnostics_folder = diagnostics_folder,
    force_rebuild=force_rebuild,
    limit_distance_urban = limit_distance_urban,
    distance_urban_quantile = distance_urban_quantile,
    remove_seasalt_dust_contribution = remove_seasalt_dust_contribution,
    limit_distance_stations = limit_distance_stations,
    distance_stations_inner_km = distance_stations_inner_km,
    distance_stations_outer_km = distance_stations_outer_km,
    distance_stations_smooth = distance_stations_smooth
  )
}


models.rf.predict.no2 <- function(obs, predictors, regions, res, year, suffix,
                                  results_folder,
                                  diagnostics_folder,
                                  limit_distance_urban,
                                  distance_urban_quantile,
                                  limit_distance_stations,
                                  distance_stations_inner_km,
                                  distance_stations_outer_km,
                                  distance_stations_smooth,
                                  remove_seasalt_dust_contribution,
                                  force_rebuild,
                                  formula=NULL) {


  no2_formula <- if(is.null(formula)){
    value ~ no2_prior + no2_omi_diff + gadm1 + srtm + grump + distance_coast + pop_ratio_log
  } else {
    formula
  }

  obs$gadm1 <- as.factor(obs$gadm1) # Should be done beforehand
  models.rf.predict.poll(
    formula = no2_formula,
    obs = obs,
    predictors = predictors,
    regions = regions,
    res = res,
    year = year,
    suffix = suffix,
    poll = "no2",
    prior = predictors$no2_prior,
    results_folder = results_folder,
    diagnostics_folder = diagnostics_folder,
    force_rebuild=force_rebuild,
    limit_distance_urban = limit_distance_urban,
    distance_urban_quantile = distance_urban_quantile,
    limit_distance_stations = limit_distance_stations,
    distance_stations_inner_km = distance_stations_inner_km,
    distance_stations_outer_km = distance_stations_outer_km,
    distance_stations_smooth = distance_stations_smooth,
    remove_seasalt_dust_contribution = remove_seasalt_dust_contribution
  )
}
