peak <- readPeakFile('/home/em_b/Downloads/GSM3692177_c13-thra1_peaks.narrowPeak.gz')
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
peakAnno <- annotatePeak(GR, tssRegion=c(-1000, 1000), 
                         TxDb=txdb, 
                         annoDb="org.Hs.eg.db")
View(data.frame(peakAnno))
