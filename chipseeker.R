library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BiocManager)
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
setwd('/home/deviancedev01/chipseq/macs3')
txdb<-TxDb.Hsapiens.UCSC.hg38.refGene
peak <- readPeakFile('GSE215081_hap1_ifny_h3k4me3_summits.bed')
table(seqnames(peak))
#head(peak)
#peak<-tidyChromosomes(peak,
#                      keep.X=FALSE,
#                      keep.Y=FALSE,
#                      keep.M = FALSE,
#                      keep.nonstandard = FALSE,
#                      genome = 'mm39')
#table(seqnames(peak))
peaks.xls <- 'GSE215081_hap1_ifny_h3k4me3_peaks.xls'
peaks.xls <- read.delim(peaks.xls,comment.char="#")
list<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
        'chr20','chr21','chr22')
peaks.xls<-subset(peaks.xls, subset = chr %in% list)
table(peaks.xls$chr)
peaks_granges<-GRanges(peaks.xls)
#-----10*log(10)*q_value
covplot(peak,
        weightCol="V5",
        xlab = "Chromosome Size (bp)",
        ylab = "-10log(10)qvalue",
        title = "GSE215081 HAP1 cells + IFNy H3K4me3 coverage plot")
peak_Profile_Heatmap(peak,
                     upstream = 3000,
                     downstream = 3000,
                     weightCol = 'V5',
                     by = "gene",
                     type = "start_site",
                     TxDb = txdb,
                     nbin = 800,
                     palette = 'RdGy',
                     ylab = NULL,
                     xlab = 'Peaks at TSS (+/-3kb)',
                     free_y = TRUE)
peakHeatmap(peaks_granges,
            TxDb = txdb,
            nbin = 800,
            upstream=3000,
            downstream=3000,
            palette = 'RdGy',
            xlab = 'Peaks at TSS (+/-3kb)')
peak_Profile_Heatmap(peak,
                     upstream = 3000,
                     downstream = 3000,
                     weightCol = 'V5',
                     by = "gene",
                     type = "body",
                     TxDb = txdb,
                     nbin = 800,
                     palette = 'RdGy',
                     ylab = NULL,
                     xlab = 'Peaks at TSS > TTS(+/-3kb)',
                     free_y = TRUE)
peakHeatmap(peak = peaks_granges,
            TxDb = txdb,
            upstream = 3000,
            downstream = 3000,
            by = "gene",
            type = "body",
            nbin = 800)+
  scale_fill_distiller(palette = 'RdGy')
gene <- seq2gene(peaks_granges,
                 tssRegion = c(-1000, 1000),
                 flankDistance = 3000,
                 TxDb=txdb)
reactome <- enrichPathway(gene,
                          organism = 'human')
dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='Reactome pathways',
        label_format=50)
edox <- setReadable(reactome, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox,
         foldChange=NULL,
         node_label='all',
         colorEdge=TRUE,
         cex_label_gene=0.5,
         circular=FALSE,
         cex_gene=0.1)
cnetplot(edox,
         foldChange=NULL,
         node_label='all',
         colorEdge=TRUE,
         cex_label_gene=0.5,
         circular=TRUE,
         cex_gene=0.1)
peakAnno <- annotatePeak(peaks_granges, tssRegion=c(-1000, 1000), 
                         TxDb=txdb, 
                         annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
peakanno_granges <- as.GRanges(peakAnno)
promoters <- peakanno_granges[peakanno_granges$annotation == 'Promoter',]
promoters <- unique(promoters$geneId)
head(promoters)
promoter_go <- enrichGO(gene = promoters,
                            OrgDb = org.Hs.eg.db,
                            ont = "BP")
dotplot(promoter_go,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=10,
        title='GO:Biological Process',
        label_format=50)
