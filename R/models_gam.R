models.gam.predict <- function(obs, poll, predictors, res, year, regions, suffix, results_folder="results/gam"){

  dir.create(results_folder, recursive = T, showWarnings = F)

  predict_fns <- list(
    "pm25"=models.gam.predict.pm25,
    "no2"=models.gam.predict.no2
  )

  predict_fns[[poll]](
    obs=obs,
    predictors=predictors,
    res=res,
    regions=regions,
    suffix=suffix,
    year=year,
    results_folder=results_folder
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
                                       results_folder) {
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

  formula <- utils.remove_constant_variables_from_formula(data=data, formula=formula)

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


models.gam.predict.pm25 <- function(obs, predictors, regions, res, year, suffix, results_folder){

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
    "CL" = diff_pm25 ~ s(lon, lat) + grump
  )

  pm25_preds <- pblapply(names(regions), function(region) {
    print(region)

    f <- glue("cache/{MODEL_GAM}/pm25_pred_{res}_{region}_{year}{suffix}.tif")
    dir.create(dirname(f), recursive = T, showWarnings = F)

    if (use_cache & file.exists(f)) {
      raster::stack(f) %>% `names<-`(c("predicted", "error"))
    } else {
      tryCatch(
        {
          region_preds <- models.gam.predict.generic(
            obs_global = obs_pm25,
            region = region,
            formula = pm25_formulas[[region]],
            predictors = predictors,
            poll = "pm25",
            res = res,
            results_folder = results_folder
          )
          writeRaster(region_preds, f, overwrite = T)
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


  pm25_diff <- models.gam.combine_predictions(pm25_preds, error_relative_threshold = 1.645)

  mask <- predictors$distance_urban < quantile(obs$distance_urban, 0.95, na.rm = T) # 0.25 deg
  mask[mask == 0] <- NA
  pm25_diff <- pm25_diff %>%
    raster::mask(mask) %>%
    raster::mask(creahelpers::to_raster(predictors$pop))

  writeRaster(pm25_diff, glue("{results_folder}/pm25_adjustment_{res}_{year}{suffix}.tif"),
              overwrite = T)

  pm25 <- raster::calc(raster::stack(list(pm25_diff, predictors$pm25_prior)), sum, na.rm = T) %>%
    raster::mask(creahelpers::to_raster(predictors$pop))

  writeRaster(pm25, glue("{results_folder}/pm25_{res}_{year}{suffix}.tif"),
              overwrite = T
  )

  pm25_noss <- pm25 * (1 - predictors$pm25_ss_dust_frac) # Remove sea salt & dust contribution
  writeRaster(pm25_noss, glue("{results_folder}/pm25_no_ss_dust_{res}_{year}{suffix}.tif"),
              overwrite = T
  )

  return(
    list(
      "pm25"=pm25,
      "pm25_noss"=pm25_noss
    )
  )
}

models.gam.predict.no2 <- function(obs, predictors, regions, res, year, suffix, results_folder){

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
    "CL" = diff_no2 ~ s(distance_urban, k = 3)
  )

  no2_preds <- pblapply(names(regions), function(region) {
    print(region)

    f <- glue("cache/{MODEL_GAM}/no2_pred_{res}_{region}_{year}{suffix}.tif")
    dir.create(dirname(f), recursive = T, showWarnings = F)

    if (file.exists(f)) {
      raster::stack(f) %>% `names<-`(c("predicted", "error"))
    } else {
      tryCatch(
        {
          region_preds <- models.gam.predict.generic(
            obs_global = obs_no2,
            region = region,
            formula = no2_formulas[[region]],
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

  no2_diff <- models.gam.combine_predictions(preds=no2_preds, error_relative_threshold = 1.645)
  mask <- predictors$distance_urban < quantile(obs$distance_urban, 0.95, na.rm = T) # 0.29 deg
  mask[mask == 0] <- NA
  no2_diff <- no2_diff %>%
    raster::mask(mask) %>%
    raster::mask(creahelpers::to_raster(predictors$pop))
  writeRaster(no2_diff, glue("{results_folder}/no2_adjustment_{res}_{year}{suffix}.tif"),
              overwrite = T
  )

  no2 <- raster::calc(raster::stack(c(predictors$no2_prior, no2_diff)), sum, na.rm = T) %>%
    raster::mask(creahelpers::to_raster(predictors$pop))
  writeRaster(no2, glue("{results_folder}/no2_adjusted_{res}_{year}{suffix}.tif"),
              overwrite = T
  )

  return(list(no2=no2))

}
