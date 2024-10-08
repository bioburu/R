#----if dependency issues use 
install.packages('pak')
pak::pak(sprintf("rstudio/%s", c("reticulate", "tensorflow", "keras")))
if (is.null(reticulate::virtualenv_starter()))
  reticulate::install_python()
tensorflow::install_tensorflow()

library(tensorflow)
tf$config$list_physical_devices("GPU")
