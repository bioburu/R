#---for blacklist regions ie mm39
library(excluderanges)
library(GenomicRanges)
library(AnnotationHub)
data(mm39.excluderanges)

ah <- AnnotationHub()
query_data <- subset(ah, preparerclass == "excluderanges")
query_data <- query(ah, c("excluderanges", "mm39"))
query_data
mm39<-query_data[['AH107321']]
mm39
#-----------------------------


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
