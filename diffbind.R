#----MUST BE ALIGNED TO MM10 OR HG38 
library(DiffBind)
library(GenomicRanges)
samples<-read.csv('/home/em_b/Desktop/diffbind_metadata.csv')
samples
dbObj<-dba(sampleSheet = samples)
dbObj
plot(dbObj)
Samples<-GRanges(samples)
test<-dba.analyze(samples,
                  method = DBA_ALL_METHODS,
                  bBlacklist = DBA_BLACKLIST_MM10,
                  bGreylist = DBA_BLACKLIST_MM10,
                  #bRetrieveAnalysis = TRUE,
                  bParallel = TRUE
                  )
cat('ERROR-----Unable to complete analysis.
Warning messages:
1: No contrasts added. There must be at least two sample groups with at least three replicates. 
2: No contrasts added. There must be at least two sample groups with at least three replicates. ')
