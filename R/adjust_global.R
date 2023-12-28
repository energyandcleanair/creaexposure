adjust_global <- function(pop,
                          res,
                          poll = c("pm25", "no2"),
                          suffix = "",
                          selected_regions = NULL,
                          year = 2022,
                          use_cache = T,
                          obs = NULL) {
  ###################################################
  # Prepare data
  ###################################################
  # Get observations and predictors and these points
  obs <- dplyr::coalesce(obs, data.get_obs(year = year, use_cache = use_cache))
  predictors <- data.predictors(pop, res, year = year, suffix = suffix, use_cache = use_cache)
  obs <- utils.add_predictors(obs, predictors)

  # obs$gadm0 <- as.factor(obs$gadm0)
  # obs$gadm1 <- as.factor(obs$gadm1)

  regions <- list(
    "NA" = c("US", "CA", "MX"),
    "IN" = "IN",
    "CN" = "CN",
    "EU" = c(
      "AD", "AT", "BA", "BE", "BG", "CH", "CY", "CZ",
      "DE", "DK", "EE", "ES", "FI", "FR", "GB", "GI", "GR",
      "HR", "HU", "IE", "IS", "IT", "LT", "LU", "LV", "ME",
      "MK", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI", "SK", "XK"
    ),
    "TR" = "TR",
    "SEA" = "TH", # No other countries with enough measurements
    "JP" = "JP",
    "ZA" = "ZA",
    "PH" = "PH",
    "TW" = "TW",
    "CL" = "CL",
    "AU" = "AU"
  )

  if (!is.null(selected_regions)) {
    regions <- regions[selected_regions]
  }

  obs$region <- NA
  for (r in names(regions)) {
    obs[obs$country %in% regions[[r]], "region"] <- r
  }

  # Visual check if some important countries haven't been covered
  # obs %>%  filter(poll=="pm25") %>%
  #   group_by(country, region) %>% summarise(count=n()) %>%
  #   mutate(name=countrycode(country, "iso2c", "country.name")) %>%
  #   arrange(desc(count)) %>%
  #   View()

  ###################################################
  # Adjust PM2.5
  ###################################################
  if ("pm25" %in% poll) {
    # Tweak for Philippines
    # We added EMB measurements, plus Lauri's estimations using various methods
    # See https://github.com/energyandcleanair/202108_hia_philippines
    obs$type <- factor(obs$type)
    type_raster <- pop
    type_raster[!is.na(type_raster)] <- which(levels(obs$type) == "measured")
    obs$type <- as.factor(as.numeric(obs$type))
    names(type_raster) <- "type"
    predictors <- predictors %>% addLayer(raster(type_raster))

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
      "PH" = diff_pm25 ~ s(pm25_prior, k = 3) + grump + s(lon, lat) + type,
      "TW" = diff_pm25 ~ s(lon, lat),
      "CL" = diff_pm25 ~ s(lon, lat) + grump
    )

    pm25_preds <- pblapply(names(regions), function(region) {
      print(region)

      f <- sprintf("cache/pm25_pred_%s_%s_%s%s.tif", res, region, year, suffix)
      if (use_cache & file.exists(f)) {
        raster::stack(f) %>% `names<-`(c("predicted", "error"))
      } else {
        tryCatch(
          {
            region_preds <- adjust_region(
              obs_global = obs_pm25,
              region = region,
              formula = pm25_formulas[[region]],
              predictors = predictors,
              poll = "pm25",
              res = res
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


    pm25_diff <- combine_preds(pm25_preds, error_relative_threshold = 1.645)

    mask <- predictors$distance_urban < quantile(obs$distance_urban, 0.95, na.rm = T) # 0.25 deg
    mask[mask == 0] <- NA
    pm25_diff <- pm25_diff %>%
      raster::mask(mask) %>%
      raster::mask(utils.to_raster(pop))

    writeRaster(pm25_diff, sprintf("results/pm25_adjustment_%s_%s%s.tif", res, year, suffix),
      overwrite = T
    )

    pm25 <- raster::calc(raster::stack(list(pm25_diff, predictors$pm25_prior)), sum, na.rm = T) %>%
      raster::mask(utils.to_raster(pop))

    writeRaster(pm25, sprintf("results/pm25_adjusted_%s_%s%s.tif", res, year, suffix),
      overwrite = T
    )

    pm25_noss <- pm25 * (1 - predictors$pm25_ss_dust_frac) # Remove sea salt & dust contribution
    writeRaster(pm25_noss, sprintf("results/pm25_adjusted_no_ss_dust_%s_%s%s.tif", res, year, suffix),
      overwrite = T
    )
  }


  ###################################################
  # Adjust NO2
  ###################################################
  if ("no2" %in% poll) {
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

      f <- sprintf("cache/no2_pred_%s_%s%s.tif", res, region, suffix)
      if (file.exists(f)) {
        raster::stack(f) %>% `names<-`(c("predicted", "error"))
      } else {
        tryCatch(
          {
            region_preds <- adjust_region(
              obs_global = obs_no2,
              region = region,
              formula = no2_formulas[[region]],
              predictors = predictors,
              poll = "no2",
              res = res
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

    no2_diff <- combine_preds(no2_preds, error_relative_threshold = 1.645)
    mask <- predictors$distance_urban < quantile(obs$distance_urban, 0.95, na.rm = T) # 0.29 deg
    mask[mask == 0] <- NA
    no2_diff <- no2_diff %>%
      raster::mask(mask) %>%
      raster::mask(utils.to_raster(pop))
    writeRaster(no2_diff, sprintf("results/no2_adjustment_%s_%s%s.tif", res, year, suffix),
      overwrite = T
    )

    no2 <- raster::calc(raster::stack(c(predictors$no2_prior, no2_diff)), sum, na.rm = T) %>%
      raster::mask(utils.to_raster(pop))
    writeRaster(no2, sprintf("results/no2_adjusted_%s_%s%s.tif", res, year, suffix),
      overwrite = T
    )
  }
}


adjust_region <- function(obs_global,
                          region,
                          formula,
                          predictors,
                          poll,
                          res,
                          gadm0_levels,
                          gadm1_levels) {
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
  sink(file = sprintf("results/gam_%s_%s_%s.txt", res, poll, region))
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

combine_preds <- function(preds, error_relative_threshold) {
  # We assume there is no overlap
  pblapply(preds, function(pred) {
    pred$valid <- abs(pred$predicted) > abs(pred$error) * error_relative_threshold
    pred$predicted[pred$valid == 0] <- NA
    pred$predicted
  }) %>%
    raster::stack() %>%
    terra::rast() %>%
    terra::app(mean, na.rm = T) %>%
    raster()
}
