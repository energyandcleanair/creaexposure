models.predict <- function(model, ...){

  predict_fns <- list()
  predict_fns[[MODEL_GAM]] <- models.gam.predict
  predict_fns[[MODEL_RF]] <- models.rf.predict

  predict_fns[[model]](...)
}
