if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsamtools")
install.packages('remotes')
if (!require(remotes)) { install.packages("remotes") }
remotes::install_github("fbreitwieser/pavian")

#---run kraken2 and view outputs
pavian::runApp(port=5000,maxUploadSize = (3000*1024^2))

