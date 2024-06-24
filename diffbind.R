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
db<-dba(sampleSheet='/home/em_b/work_stuff/diffbind/diffbind_meta_h3k27me3.csv')
db
plot(db)
counts<-dba.count(db,
                  #peaks = NULL,
                  score = DBA_SCORE_READS_MINUS,
                  bParallel =FALSE)
df<-dba.peakset(counts,
                bRetrieve = TRUE)
df<-data.frame(df)
#write.csv(df,file='h3k27me3_matrix.csv')
