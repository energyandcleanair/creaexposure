models.predict <- function(model, ...){

  predict_fns <- list(
    MODEL_GAM = models.gam.predict
  )

  predict_fns[[model]](...)
}
