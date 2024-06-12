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
db
plot(db)
counts<-dba.count(db,
                  #filter = 0,
                  #summits =FALSE, 
                  bParallel =FALSE)
counts
plot(counts)
counts_norm<-dba.normalize(counts)
counts_norm$norm
counts_norm<-dba.contrast(counts_norm,
                          contrast = c('Treatment','none','treated'))
counts_norm$contrasts
counts_norm<-dba.analyze(counts_norm)
dba.show(counts_norm,bContrasts = TRUE)
plot(counts_norm,contrast=1)
dba.report(counts_norm)
counts_norm[["contrasts"]][[1]][["DESeq2"]][["de"]]
