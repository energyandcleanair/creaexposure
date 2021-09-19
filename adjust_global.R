adjust_global <- function(pop, res){

  # Get observations and predictors and these points
  obs <- data.get_obs(year=2019)
  predictors <- data.predictors(pop, res)
  obs <- cbind(obs,
               raster::extract(predictors, sf::st_as_sf(obs)))


  # Group regions
  obs %>% group_by(country) %>% summarise(count=n()) %>%
    mutate(name=countrycode(country, "iso2c", "country.name")) %>%
    arrange(desc(count)) %>%
    View()


  regions <- list(
    "NA"=c("US","CA","MX"),
    "IN"="IN",
    "CN"="CN",
    "EU"=c("AD","AT","BA","BE","BG","CH","CY","CZ",
                     "DE","DK","EE","ES","FI","FR","GB","GI","GR",
                     "HR","HU","IE","IS","IT","LT","LU","LV","ME",
                     "MK","MT","NL","NO","PL","PT","RO","RS","SE","SI","SK","TR","XK"),
    "SEA"="TH", # No other countries with enough measurements
    "JP"="JP",
    "ZA"="ZA"
  )

  obs$region <- NA
  for(r in names(regions)){
    obs[obs$country %in% regions[[r]], "region"] <- r
  }

  # Adjust PM2.5
  obs_pm25 <- obs %>% filter(poll=="pm25") %>% mutate(diff_pm25=value-pm25_prior)

  # IT HAS SIGNIFICANT IMPLICATIONS TO INCLUDE PRIOR
  # It means prior is systematically biased e.g. overestimating/underestimating high values
  # + risk of overfitting
  # We've manually checked smooth terms for s(pm25_prior) for every fitting and checked
  # it looked alright (e.g. low edf, consistent trend)

  pm25_formulas <- list(
    # "CN" = diff_pm25 ~  s(pm25_prior, k=3) + grump + s(lon,lat) + factor(gadm1),
    # "IN" = diff_pm25 ~  s(pm25_prior, k = 3) + s(pm25_ss_dust_frac) + s(lon,lat) + s(no2_prior),
    # "NA" = diff_pm25 ~  s(srtm, k = 4) + s(pm25_prior, k = 3) + s(no2_prior) +  s(lon, lat),
    # "SEA" = diff_pm25 ~  s(distance_coast) + s(pm25_prior, k=3),
    # "EU" = diff_pm25 ~  s(pop) + s(distance_coast) + grump + s(lon, lat) + s(pm25_prior),
    "JP" = diff_pm25 ~  s(pm25_prior) + s(lon,lat),
    "ZA" = diff_pm25 ~ s(distance_coast) + s(pm25_prior)
  )

  pm25_preds <- pblapply(names(pm25_formulas), function(region){
    print(region)
    adjust_region(
      obs_global = obs_pm25,
      region=region,
      formula=pm25_formulas[[region]],
      predictors=predictors,
      poll="pm25",
      res=res)
  })
}


adjust_region <- function(obs_global, region, formula, predictors, poll, res){

  #####################
  # Prepare data
  #####################
  data <- obs_global %>%
    filter(region==!!region) %>%
    cbind(sf::st_coordinates(.$geometry)) %>%
    dplyr::select_at(c("geometry","gadm0", all.vars(formula))) %>%
    drop_na()

  #####################
  # Train
  #####################
  model <- gam(formula, data=data)
  summary(model)

  #####################
  # Quick diagnosis
  #####################
  sink(file = sprintf("results/gam_%s_%s_%s.txt", res, poll, region))
  summary(model)
  sink(file=NULL)
  p <- predict(model, data, se.fit=T)
  data$predicted <- p$fit
  data$se <- p$se.fit

  ggplot(data) +
    geom_point(aes_string(as.character(formula)[2], "predicted", col="se")) +
    geom_abline()

  #####################
  # Predict
  #####################
  mask_region <- predictors$gadm0 %in% unique(data$gadm0)
  mask_region[mask_region==0] <- NA

  predfun <- function(model, data) {
    v <- predict.gam(model, data, se.fit=TRUE)
    cbind(p=as.vector(v$fit), se=as.vector(v$se.fit))
  }

  predictors_slim <- predictors[[all.vars(formula)]] %>%
    mask(mask_region)

  # You can use multiple cores to speed up the predict function
  # by calling it via the clusterR function (you may need to install the snow package)
  tic()
  beginCluster(parallel::detectCores()-1)
  pred <- clusterR(predictors_slim, predict, args=list(model=model, fun=predfun, index=1:2, na.rm=T, inf.rm=T))
  endCluster()
  toc()

  names(pred) <- c("predicted", "error")
  return(pred)
}
