models.prior.predict <- function(obs,
                                 poll,
                                 predictors,
                                 ...) {

  if (poll == "pm25") {
    prior <- predictors$pm25_prior
  } else if (poll == "no2") {
    prior <- predictors$no2_prior
  } else {
    stop("Unknown poll")
  }

  return(creahelpers::to_rast(prior))
}
