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

library(DiffBind)
db<-dba(sampleSheet='/home/em_b/Desktop/diffbind/diffbind_meta.csv')
plot(db)
counts<-dba.count(db,
                  summits=250,
                  bParallel =FALSE)
