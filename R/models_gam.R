models.gam.predict <- function(
    obs,
    poll,
    predictors,
    res,
    year,
    regions,
    suffix,
    force_rebuild = T,
    results_folder = "results/gam",
    remove_seasalt_dust_contribution = T,
    limit_distance_urban = T,
    distance_urban_quantile = 0.95,
    ...) {

  dir.create(results_folder, recursive = T, showWarnings = F)

  predict_fns <- list(
    "pm25" = models.gam.predict.pm25,
    "no2" = models.gam.predict.no2
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


models.gam.predict.generic <- function(obs_global,
                                       region,
                                       formula,
                                       predictors,
                                       poll,
                                       res,
                                       gadm0_levels,
                                       gadm1_levels,
                                       results_folder,
                                       ...) {
  #####################
  # Prepare data
  #####################
  data <- obs_global %>%
    filter(region == !!region) %>%
    cbind(sf::st_coordinates(.$geometry)) %>%
    dplyr::select_at(c("geometry", "gadm0", "value", all.vars(formula))) %>%
    drop_na()

  data <- data %>%
    mutate_at(intersect(names(.), c("gadm0", "gadm1")), as.factor)

  formula <- utils.remove_constant_variables_from_formula(data = data, formula = formula)

  #####################
  # Train
  #####################
  model <- gam(formula, data = data)
  summary(model)

  #####################
  # Quick diagnosis
  #####################
  sink(file = glue("{results_folder}/gam_{res}_{poll}_{region}.txt"))
  print(summary(model))
  sink(file = NULL)
  p <- predict(model, data, se.fit = T)
  data$predicted_model <- p$fit
  data$se <- p$se.fit

  ggplot(data) +
    geom_point(aes_string(as.character(formula)[2], "predicted_model", col = "se")) +
    geom_abline()

  #####################
  # Predict
  #####################
  mask_region <- raster::match(predictors$gadm0, unique(data$gadm0))
  mask_region[mask_region == 0] <- NA

  predfun <- function(model, data) {
    v <- mgcv::predict.gam(model, data, se.fit = TRUE)
    cbind(p = as.vector(v$fit), se = as.vector(v$se.fit))
  }

  predictors_slim <- raster::subset(predictors, intersect(names(predictors), all.vars(formula))) %>%
    raster::mask(mask_region)

  pred <- raster::predict(predictors_slim, model = model, fun = predfun, index = 1:2, na.rm = T)

  # Check predicted values by raster and model are equal (risks when factors involvedd)
  check <- data %>% utils.add_predictors(pred)
  if (!all(check$predicted_model == check$predicted)) {
    stop("Something wrong in prediction. Probably has to do with factor levels")
  }

  return(pred)
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
models.gam.combine_predictions <- function(preds, error_relative_threshold) {
  # We assume there is no overlap
  pblapply(preds, function(pred) {
    pred$valid <- abs(pred$predicted) > abs(pred$error) * error_relative_threshold
    pred$predicted[pred$valid == 0] <- NA
    pred$predicted
  }) %>%
    raster::stack() %>%
    terra::rast() %>%
    terra::app(mean, na.rm = T) %>%
    raster::raster()
}

models.gam.result_filepath <- function(poll, res, year, suffix, folder, extension = "tif"){
  glue("{folder}/{poll}_{res}_{year}{suffix}.{extension}")
}


models.gam.predict.pm25 <- function(obs, predictors, regions, res, year, suffix, results_folder, force_rebuild,
                                      limit_distance_urban, distance_urban_quantile, remove_seasalt_dust_contribution,
                                      ...

) {

  filepath <- models.rf.result_filepath(
    poll = "pm25",
    res = res,
    year = year,
    suffix = suffix,
    folder = results_folder
  )

  if(file.exists(filepath) & !force_rebuild){
    return(raster::raster(filepath))
  }

  obs_pm25 <- obs %>%
    filter(poll == "pm25") %>%
    mutate(diff_pm25 = value - pm25_prior)


  # IT HAS SIGNIFICANT IMPLICATIONS TO INCLUDE PRIOR
  # It means prior is systematically biased e.g. overestimating/underestimating high values
  # + risk of overfitting
  # We've manually checked smooth terms for s(pm25_prior) for every fitting and checked
  # it looked alright (e.g. low edf, consistent trend)
  pm25_formulas <- list(
    "CN" = diff_pm25 ~ s(pm25_prior, k = 3) + s(pm25_merra2_diff, k = 3) + s(distance_urban, k = 3) + s(lon, lat) + gadm1,
    "IN" = diff_pm25 ~ s(pm25_prior, k = 3) + s(pm25_ss_dust_frac) + s(lon, lat) + s(no2_prior),
    "NA" = diff_pm25 ~ s(srtm, k = 4) + s(pm25_prior, k = 3) + s(no2_prior) + s(lon, lat),
    "SEA" = diff_pm25 ~ s(pop_ratio_log, k = 3),
    "EU" = diff_pm25 ~ diff_pm25 ~ s(pop, k = 3) + s(distance_coast, k = 3) +
      s(distance_urban, k = 3) + s(pm25_prior, k = 3) + gadm0,
    "TR" = diff_pm25 ~ s(lat, lon),
    "JP" = diff_pm25 ~ s(pm25_prior) + s(lon, lat),
    "ZA" = diff_pm25 ~ s(distance_coast) + s(pm25_prior),
    "PH" = diff_pm25 ~ s(pm25_prior, k = 3) + grump + s(lon, lat),
    "TW" = diff_pm25 ~ s(lon, lat),
    "CL" = diff_pm25 ~ s(lon, lat) + grump,
    "default" = diff_pm25 ~ s(pm25_prior, k = 3) + s(pm25_ss_dust_frac) + s(lon, lat) + s(no2_prior)
  )

  pm25_preds <- pblapply(get_region_names(regions), function(region) {
      tryCatch(
        {
          formula <- default_if_null(pm25_formulas[[region]], pm25_formulas[["default"]])
          region_preds <- models.gam.predict.generic(
            obs_global = obs_pm25,
            region = region,
            formula = formula,
            predictors = predictors,
            poll = "pm25",
            res = res,
            results_folder = results_folder
          )
          names(region_preds) <- c("predicted", "error")
          region_preds
        },
        error = function(e) {
          warning("Failed to adjust for region ", region, ": ", e)
          return(NA)
        }
      )
    })

  pm25_diff <- models.gam.combine_predictions(pm25_preds, error_relative_threshold = 1.645)

  pm25_diff <- utils.mask_far_from_urban(
    r = pm25_diff,
    predictors = predictors,
    obs = obs,
    quantile = 0.95
  )

  pm25 <- raster::calc(raster::stack(list(pm25_diff, predictors$pm25_prior)), sum, na.rm = T) %>%
    raster::mask(creahelpers::to_raster(predictors$pop))

  # Remove far from urban
  if (limit_distance_urban) {
    pm25 <- utils.mask_far_from_urban(
      r = pm25,
      predictors = predictors,
      obs = obs,
      quantile = distance_urban_quantile
    )
  }

  # Remove sea salt contribution
  if (remove_seasalt_dust_contribution) {
    pm25 <- pm25 * (1 - predictors$pm25_ss_dust_frac)
  }

  writeRaster(pm25, filepath, overwrite=T)

  return(pm25)
}

models.gam.predict.no2 <- function(obs, predictors, regions, res, year, suffix, results_folder, force_rebuild,
                                      limit_distance_urban, distance_urban_quantile, remove_seasalt_dust_contribution,
                                      ...) {
  obs_no2 <- obs %>%
    filter(poll == "no2") %>%
    mutate(diff_no2 = value - no2_prior)

  # IT HAS SIGNIFICANT IMPLICATIONS TO INCLUDE PRIOR
  # It means prior is systematically biased e.g. overestimating/underestimating high values
  # + risk of overfitting
  # We've manually checked smooth terms for s(pm25_prior) for every fitting and checked
  # it looked alright (e.g. low edf, consistent trend)

  # We add a distance_urban to all of these so that
  # values far from urban center will be excluded (error too large)
  no2_formulas <- list(
    "CN" = diff_no2 ~ s(pop_ratio_log, k = 3) + s(lon, lat) + gadm1, # Shanghai is way lower than prior. GADM1 helps a lot
    "IN" = diff_no2 ~ s(lon, lat) + s(pop_ratio_log, k = 3),
    "NA" = diff_no2 ~ s(lon, lat) + s(no2_prior, k = 3),
    "SEA" = diff_no2 ~ 1, # Couldn't find any better
    "EU" = diff_no2 ~ s(no2_omi_diff, k = 3) + s(distance_urban, k = 3) +
      s(pop, k = 3) + s(distance_coast, k = 3) + gadm0 + s(lon, lat),
    "TR" = diff_no2 ~ s(no2_omi_diff, k = 3),
    "JP" = diff_no2 ~ s(no2_prior, k = 3) + s(pop_ratio_log, k = 3),
    "ZA" = diff_no2 ~ s(pop_ratio_log, k = 3),
    "TW" = diff_no2 ~ s(no2_prior, k = 3) + s(distance_urban, k = 3),
    "CL" = diff_no2 ~ s(distance_urban, k = 3),
    "default" = diff_no2 ~ s(lon, lat) + s(no2_prior, k = 3),
  )

  no2_preds <- pblapply(get_region_names(regions), function(region) {
    print(region)

    f <- glue("cache/{MODEL_GAM}/no2_pred_{res}_{region}_{year}{suffix}.tif")
    dir.create(dirname(f), recursive = T, showWarnings = F)

    if (file.exists(f)) {
      raster::stack(f) %>% `names<-`(c("predicted", "error"))
    } else {
      tryCatch(
        {
          formula <- default_if_null(no2_formulas[[region]], no2_formulas[["default"]])
          region_preds <- models.gam.predict.generic(
            obs_global = obs_no2,
            region = region,
            formula = formula,
            predictors = predictors,
            poll = "no2",
            res = res,
            results_folder = results_folder
          )
          writeRaster(region_preds, f)
          names(region_preds) <- c("predicted", "error")
          region_preds
        },
        error = function(e) {
          warning("Failed to adjust for region ", region, ": ", e)
          return(NA)
        }
      )
    }
  })

  no2_diff <- models.gam.combine_predictions(preds = no2_preds, error_relative_threshold = 1.645)

  no2_diff <- utils.mask_far_from_urban(
    r = no2_diff,
    predictors = predictors,
    obs = obs,
    quantile = 0.95
  )

  writeRaster(no2_diff, glue("{results_folder}/no2_adjustment_{res}_{year}{suffix}.tif"),
    overwrite = T
  )

  no2 <- raster::calc(raster::stack(c(predictors$no2_prior, no2_diff)), sum, na.rm = T) %>%
    raster::mask(creahelpers::to_raster(predictors$pop))
  writeRaster(no2, glue("{results_folder}/no2_adjusted_{res}_{year}{suffix}.tif"),
    overwrite = T
  )

  return(list(no2 = no2))
}
