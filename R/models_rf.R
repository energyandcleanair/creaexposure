models.rf.predict <- function(obs,
                              poll,
                              predictors,
                              res,
                              year,
                              regions,
                              suffix,
                              force_rebuild = T,
                              results_folder = "results/rf",
                              remove_seasalt_dust_contribution = T,
                              limit_distance_urban = T,
                              distance_urban_quantile = 0.95,
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
    limit_distance_urban = limit_distance_urban,
    distance_urban_quantile = distance_urban_quantile,
    remove_seasalt_dust_contribution = remove_seasalt_dust_contribution,
    force_rebuild = force_rebuild
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

  # Export diagnostics
  models.rf.diagnose(
    model = model,
    diagnostics_folder = diagnostics_folder,
    res = res,
    poll = poll,
    region = region,
    year = year,
    suffix = suffix
  )

  # Predict
  pred <- raster::predict(predictors, model)

  # Remove far from urban
  if (limit_distance_urban) {
    pred <- utils.mask_far_from_urban(
      r = pred,
      predictors = predictors, obs = obs,
      quantile = distance_urban_quantile
    )
  }

  # Remove sea salt contribution
  if (remove_seasalt_dust_contribution & poll == "pm25") {
    pred <- pred * (1 - predictors$pm25_ss_dust_frac)
  }

  # Take prior where pred is NA
  pred[is.na(pred)] <- prior[is.na(pred)]

  # Return predictions
  return(pred)
}


models.rf.diagnose <- function(model, diagnostics_folder, res, poll, region, year, suffix) {
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
    head(20) %>%
    # add rownames as a column
    tibble::rownames_to_column("var") %>%
    mutate(var = fct_reorder(var, IncNodePurity)) %>%
    ggplot() +
    geom_col(aes(y = var, x = IncNodePurity)) -> plt


  ggsave(filepath %>% gsub("\\.txt", ".png", .), plt, width = 10, height = 10)
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
                                   force_rebuild = T,
                                   limit_distance_urban = T,
                                   distance_urban_quantile = 0.95,
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
            limit_distance_urban = limit_distance_urban,
            distance_urban_quantile = distance_urban_quantile,
            remove_seasalt_dust_contribution = remove_seasalt_dust_contribution,
          )
  }) %>%
    raster::stack() %>%
    terra::rast() %>%
    terra::app(mean, na.rm = T) %>%
    raster() %>%
    `names<-`("predicted")

  # Save
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
                                   regions, res, year, suffix, results_folder,
                                   limit_distance_urban,
                                   distance_urban_quantile,
                                   remove_seasalt_dust_contribution,
                                   force_rebuild,
                                   ...) {
  pm25_formula <- value ~ pm25_prior + pm25_merra2_diff +
    gadm1 + distance_urban +
    srtm + srtm_diff05deg + srtm_05deg +
    grump + distance_coast +
    pop_ratio_log + pop

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
    force_rebuild=force_rebuild,
    limit_distance_urban = limit_distance_urban,
    distance_urban_quantile = distance_urban_quantile,
    remove_seasalt_dust_contribution = remove_seasalt_dust_contribution
  )
}


models.rf.predict.no2 <- function(obs, predictors, regions, res, year, suffix, results_folder,
                                  limit_distance_urban,
                                  distance_urban_quantile,
                                  remove_seasalt_dust_contribution,
                                  force_rebuild) {
  no2_formula <- value ~ no2_prior + no2_omi_diff + gadm1 + srtm + grump + distance_coast + pop_ratio_log
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
    force_rebuild=force_rebuild,
    limit_distance_urban = limit_distance_urban,
    distance_urban_quantile = distance_urban_quantile,
    remove_seasalt_dust_contribution = remove_seasalt_dust_contribution
  )
}
