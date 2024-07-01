library(BiocParallel)

param<-MulticoreParam(workers=12,
                      progressbar=TRUE)
register(param)

#---maybe this will work?

