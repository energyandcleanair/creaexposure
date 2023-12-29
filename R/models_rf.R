models.rf.predict <- function(obs, poll, predictors, res, year, regions, suffix, results_folder="results/rf") {

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
    results_folder = results_folder
  )
}


models.rf.predict.generic <- function(formula, obs, region, predictors, prior, poll, res, year, diagnostics_folder="diagnostics/rf") {

  # Prepare data
  data <- obs %>%
    filter(region == !!region) %>%
    filter(poll == !!poll)

  # Keep only variables in formula
  data <- data %>%
    dplyr::select_at(c("value", all.vars(formula))) %>%
    drop_na()

  # Define control parameters for training
  control <- trainControl(method = "cv", number = 10) # Example: 10-fold cross-validation

  # Train Random Forest Model using caret
  model <- caret::train(formula, data = data, method = "rf", trControl = control)

  # Export diagnostics
  models.rf.diagnose(model=model, diagnostics_folder=diagnostics_folder, res=res, poll=poll, region=region,
                     year=year)

  # Predict
  pred <- raster::predict(predictors, model)

  # Remove far from urban
  pred <- utils.mask_far_from_urban(r=pred, predictors=predictors, obs=obs)

  # Take prior where pred is NA
  pred[is.na(pred)] <- prior[is.na(pred)]

  # Return predictions
  return(pred)
}

models.rf.diagnose <- function(model, diagnostics_folder, res, poll, region, year){
  dir.create(diagnostics_folder, recursive = T, showWarnings = F)
  sink(file = glue("{diagnostics_folder}/rf_{res}_{poll}_{region}_{year}.txt"))
  print(model)
  print(summary(model))
  sink(file = NULL)

  model$finalModel$importance %>%
    as.data.frame() %>%
    arrange(desc(IncNodePurity)) %>%
    head(20) %>%
    # add rownames as a column
    tibble::rownames_to_column("var") %>%
    mutate(var=fct_reorder(var, IncNodePurity)) %>%
    ggplot() +
    geom_col(aes(y=var, x=IncNodePurity)) -> plt


  ggsave(glue("{diagnostics_folder}/rf_{res}_{poll}_{region}_{year}.png"), plt, width=10, height=10)
}

#' Mosaic/Average predictions of different regions
#'
#' @param preds
#' @param error_relative_threshold
#'
#' @return
#' @export
#'
#' @examples
models.rf.combine_predictions <- function(preds, error_relative_threshold) {
  # We assume there is no overlap
  preds %>%
    raster::stack() %>%
    terra::rast() %>%
    terra::app(mean, na.rm = T) %>%
    raster()
}


models.rf.predict.poll <- function(obs, predictors, regions, res, year, suffix, poll, formula, prior, use_cache=F, results_folder) {

  preds <- pblapply(names(regions), function(region) {
    print(region)

    f <- glue("cache/{MODEL_RF}/{poll}_pred_{res}_{region}_{year}{suffix}.tif")
    dir.create(dirname(f), recursive = T, showWarnings = F)

    if (use_cache & file.exists(f)) {
      raster::stack(f) %>% `names<-`(c("predicted"))
    } else {
      tryCatch(
        {
          region_preds <- models.rf.predict.generic(
            formula=formula,
            obs = obs %>% filter(poll == !!poll),
            region = region,
            predictors = predictors,
            poll = poll,
            res = res,
            prior = prior,
            year=year
          )
          writeRaster(region_preds, f, overwrite = T)
          names(region_preds) <- c("predicted")
          region_preds
        },
        error = function(e) {
          warning("Failed to adjust for region ", region, ": ", e)
          return(NA)
        }
      )
    }
  })

  # Mosaic/combine predictions
  pred <- models.rf.combine_predictions(preds)

  mask <- predictors$distance_urban < quantile(obs$distance_urban, 0.95, na.rm = T) # 0.25 deg
  mask[mask == 0] <- NA
  pred <- pred %>%
    raster::mask(mask) %>%
    raster::mask(utils.to_raster(predictors$pop))

  writeRaster(pred, glue("{results_folder}/{poll}_{res}_{year}{suffix}.tif"),
    overwrite = T
  )

  return(pred)
}


models.rf.predict.pm25 <- function(obs, predictors, regions, res, year, suffix, results_folder) {

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
    results_folder = results_folder
  )
}


models.rf.predict.no2 <- function(obs, predictors, regions, res, year, suffix, results_folder) {
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
    results_folder = results_folder
  )
}
