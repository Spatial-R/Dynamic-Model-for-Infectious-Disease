library(ABSEIR)
load("dataset.RData")
t2 <- system.time({
  weser.model2 <- SpatialSEIRModel(data_model = weser.data_model,
                                   exposure_model = weser.exposure_model,
                                   reinfection_model = weser.reinfection_model,
                                   distance_model = weser.CAR_model,
                                   transition_priors = weser.transition_priors,
                                   initial_value_container = weser.initial_values,
                                   sampling_control = weser.sampling_control.init,
                                   samples = 100, 
                                   verbose = FALSE)
  weser.model2 <- update(weser.model2, 
                         sampling_control = weser.sampling_control, 
                         verbose=FALSE)
  save(weser.model2, file = "weserModel2.rda")
}
)