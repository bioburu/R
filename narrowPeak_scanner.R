library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)
library(plyranges)
library(BRGenomics)
library(ggplot2)
library(UpSetR)
library(ChIPpeakAnno)
library(pathview)
txdb<-TxDb.Hsapiens.UCSC.hg38.knownGene
peak <- readPeakFile('/home/em_b/Desktop/gse129039/GSM3692173_ad11-thra1_peaks.narrowPeak.gz')
table(seqnames(peak))
peak
peak<-tidyChromosomes(peak,
                      keep.X=FALSE,
                      keep.Y=FALSE,
                      keep.M = FALSE,
                      keep.nonstandard = FALSE,
                      genome = 'hg38')
table(seqnames(peak))
peak<-data.frame(peak)
GR <- GRanges(seqnames=peak$seqnames,
              IRanges(peak$start,
                      peak$end))
peakAnno <- annotatePeak(GR, tssRegion=c(-1000, 1000), 
                         TxDb=txdb, 
                         annoDb="org.Hs.eg.db")
df<-(data.frame(peakAnno))
genes<-df[!duplicated(df$SYMBOL),]
cat(genes$SYMBOL)
View(df)
break 
