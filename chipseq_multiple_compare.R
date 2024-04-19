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
setwd('/home/em_b/work_stuff/FCCC/chipseq/comb/macs/')
txdb<-TxDb.Mmusculus.UCSC.mm39.refGene
files<-list(B22.Me3_1='/home/em_b/work_stuff/FCCC/chipseq/comb/macs/E2/E2_A2_peaks.narrowPeak',
         B22.Me3_2='/home/em_b/work_stuff/FCCC/chipseq/comb/macs/G2/G2_B2_peaks.narrowPeak',
         B22.Me3_3='/home/em_b/work_stuff/FCCC/chipseq/comb/macs/I2/I2_C2_peaks.narrowPeak',
         T3.Me3_1='/home/em_b/work_stuff/FCCC/chipseq/comb/macs/E3/E3_A3_peaks.narrowPeak',
         T3.Me3_2='/home/em_b/work_stuff/FCCC/chipseq/comb/macs/G3/G3_B3_peaks.narrowPeak',
         T3.Me3_3='/home/em_b/work_stuff/FCCC/chipseq/comb/macs/I3/I3_C3_peaks.narrowPeak')
plotPeakProf2(files, upstream = 3000, downstream = 3000, conf = 0.95,
              by = "gene", type = "start_site", TxDb = txdb,
              facet = "row")
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
plotPeakProf2(files, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body",
              TxDb = txdb, facet = "row", nbin = 800)
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)
#genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
#vennplot(genes)
break 
peakHeatmap(peak = files,
            TxDb = txdb,
            upstream = rel(0.2),
            downstream = rel(0.2),
            by = "gene",
            type = "body",
            nbin = 800)+
  scale_fill_distiller(palette = 'Greys')
