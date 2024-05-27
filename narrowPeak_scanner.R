library(Gviz)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm39)
library(ChIPseeker)
library(BRGenomics)
txdb<-TxDb.Mmusculus.UCSC.mm39.knownGene
#---------------------
peak <- readPeakFile('/home/em_b/work_stuff/chipseq/thra_0hr_comb/peaks/thra_0hr_IgG_peaks.narrowPeak')
table(seqnames(peak))
peak
peak<-tidyChromosomes(peak,
                      keep.X=FALSE,
                      keep.Y=FALSE,
                      keep.M = FALSE,
                      keep.nonstandard = FALSE,
                      genome = 'mm39')
table(seqnames(peak))
peak<-data.frame(peak)
GR <- GRanges(seqnames=peak$seqnames,
              IRanges(peak$start,
                      peak$end))
GR
peakAnno <- annotatePeak(GR, tssRegion=c(-1000, 1000), 
                         TxDb=txdb, 
                         annoDb="org.Mm.eg.db")
View(data.frame(peakAnno))
df<-(data.frame(peakAnno))
df<-cbind(df,peak)

df<-df[order(df$V5, decreasing=TRUE),]
