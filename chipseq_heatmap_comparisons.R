library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(ChIPseeker)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ReactomePA)
library(plyranges)
library(BRGenomics)
library(ggplot2)
library(UpSetR)
library(ChIPpeakAnno)
library(pathview)
#setwd('/home/em_b/work_stuff/chipseq2/thra_0hr_comb/peaks')
txdb<-TxDb.Mmusculus.UCSC.mm39.refGene
peak <- readPeakFile('/home/em_b/Desktop/chip_validations_final/in_house/peaks_q0.01/thra_0hr_comb_peaks.narrowPeak')
table(seqnames(peak))
peak
peak<-tidyChromosomes(peak,
                      keep.X=FALSE,
                      keep.Y=FALSE,
                      keep.M = FALSE,
                      keep.nonstandard = FALSE,
                      genome = 'mm39')
table(seqnames(peak))
View(data.frame(peak))
covplot(peak, weightCol="V5")
peak_Profile_Heatmap(peak,
                     upstream = 8000,
                     downstream = 8000,
                     weightCol = 'V5',
                     by = "gene",
                     type = "body",
                     TxDb = txdb,
                     nbin = 800,
                     palette = 'Greys',
                     ylab = 'untreated',
                     xlab = 'Peaks at TSS > TTS (+/-3kb)',
                     free_y = TRUE)
#-------------------------------------------------------------------------------
peak <- readPeakFile('/home/em_b/Desktop/chip_validations_final/in_house/peaks_q0.01/thra_T3_comb_peaks.narrowPeak')
table(seqnames(peak))
peak
peak<-tidyChromosomes(peak,
                      keep.X=FALSE,
                      keep.Y=FALSE,
                      keep.M = FALSE,
                      keep.nonstandard = FALSE,
                      genome = 'mm39')
table(seqnames(peak))
View(data.frame(peak))
covplot(peak, weightCol="V5")
peak_Profile_Heatmap(peak,
                     upstream = 8000,
                     downstream = 8000,
                     weightCol = 'V5',
                     by = "gene",
                     type = "body",
                     TxDb = txdb,
                     nbin = 800,
                     palette = 'Greys',
                     ylab = 'T3_treated',
                     xlab = 'Peaks at TSS > TTS (+/-3kb)',
                     free_y = TRUE)
#-------------------------------------------------------------------------------
files<-list(THRA_untreated='/home/em_b/Desktop/chip_validations_final/in_house/peaks_q0.01/thra_0hr_comb_peaks.narrowPeak',
            THRA_T3='/home/em_b/Desktop/chip_validations_final/in_house/peaks_q0.01/thra_T3_comb_peaks.narrowPeak')
files
peakAnnoList <- lapply(files, annotatePeak, TxDb=TxDb.Mmusculus.UCSC.mm39.knownGene,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
peakAnnoList
plotAnnoBar(peakAnnoList)
