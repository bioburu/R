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
setwd('/home/em_b/work_stuff/chipseq2/thra_0hr_comb/peaks_p0.01')
txdb<-TxDb.Mmusculus.UCSC.mm39.refGene
peak <- readPeakFile('thra_0hr_comb_peaks.narrowPeak')
table(seqnames(peak))
peak
peak<-tidyChromosomes(peak,
                      keep.X=FALSE,
                      keep.Y=FALSE,
                      keep.M = FALSE,
                      keep.nonstandard = FALSE,
                      genome = 'mm39')
table(seqnames(peak))
covplot(peak, weightCol="V5")
summary(peak$V5)
#covplot(peak, weightCol="V5", chrs=c("chr1", "chr19"))
gc()
plotAvgProf2(peak, TxDb=txdb, upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency",
             conf = 0.95)
plotPeakProf2(peak = peak, upstream = 3000, downstream = 3000,
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = txdb, weightCol = "V5",ignore_strand = F)
gc()
macsPeaks <- 'thra_0hr_comb_peaks.xls'
#htmapdf <- read.delim(macsPeaks,comment.char="#")
#summary(htmapdf$X.log10.qvalue.)
#df<-subset(df,X.log10.qvalue.>84)
#htmapdf<-GRanges(htmapdf)
#htmapdf<-tidyChromosomes(htmapdf,
#                         keep.X=FALSE,
#                         keep.Y=FALSE,
#                         keep.M = FALSE,
#                         keep.nonstandard = FALSE)
#table(htmapdf@seqnames)
#peakHeatmap(htmapdf,
#            TxDb = TxDb.Mmusculus.UCSC.mm39.knownGene,
#            nbin = 800,
#            upstream=3000,
#            downstream=3000)
#peakHeatmap(peak = htmapdf,
#            TxDb = txdb,
#            upstream = rel(0.2),
#            downstream = rel(0.2),
#            by = "gene",
#            type = "body",
#            nbin = 800)+
#  scale_fill_distiller(palette = 'Greys')
#--------------enter in peaks.xls
peaks.xls <- 'thra_0hr_comb_peaks.xls'
#------------read in peaks file
peaks.xls <- read.delim(peaks.xls,comment.char="#")
table(peaks.xls$chr)
#----------remove chromosomes not needed
list<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19')
peaks.xls<-subset(peaks.xls, subset = chr %in% list)
table(peaks.xls$chr)
reactome<-GRanges(peaks.xls)
cat('Reactome is for all peaks general output')
gene <- seq2gene(reactome,
                 tssRegion = c(-1000, 1000),
                 flankDistance = 3000,
                 TxDb=txdb)
length(gene)
head(gene)
str(gene)
reactome <- enrichPathway(gene,
                          organism = 'mouse')
reactome
dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=15,
        title='q-values derived from modeling analysis of chip-seq 3
peak calls',
        label_format=50)
Reactome<-data.frame(reactome)
str(gene)
edox <- setReadable(reactome, 'org.Mm.eg.db', 'ENTREZID')
cnetplot(edox,
         foldChange=gene,
         node_label='all',
         colorEdge=TRUE)
cnetplot(edox,
         foldChange=gene,
         node_label='all',
         colorEdge=TRUE,
         circular=TRUE)
#-------------------------------------------------------------------------------
table(peaks.xls$chr)
GR <- GRanges(seqnames=peaks.xls[,"chr"],
              IRanges(peaks.xls[,"start"],
                      peaks.xls[,"end"]))
peakAnno <- annotatePeak(GR, tssRegion=c(-1000, 1000), 
                         TxDb=txdb, 
                         annoDb="org.Mm.eg.db")
go <- as.GRanges(peakAnno)
table(go$annotation)
plotAnnoPie(peakAnno)
#-------------------------------------------------------------------------------
promoters <- go[go$annotation == 'Promoter',]
promoters <- unique(promoters$geneId)
head(promoters)
promoter_result <- enrichGO(gene = promoters,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP")
GO_plot <- pairwise_termsim(promoter_result)
emapplot(GO_plot,
         color='p.adjust',
         layout.params = list(layout = 'dh'),
         edge.params = list(min = 0.5),
         cex.params = list(line = 0.1),
         pie='equal',
         shadowtext=TRUE,
         repel=TRUE,
         cex_label_category=1,
         showCategory = 50)
goplot(GO_plot,
       showCategory = 17,
       color = "p.adjust",
       layout = "sugiyama",
       geom = "text")
View(data.frame(promoter_result))
#-------------------------------------------------------------------------------
distal_intergenic <- go[go$annotation == "Distal Intergenic",]
distal_intergenic <- unique(distal_intergenic$geneId)
head(distal_intergenic)
DI_result <- enrichGO(gene = distal_intergenic,
                      #universe = allGeneGR,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP")
GO_plot <- pairwise_termsim(DI_result)
emapplot(GO_plot,
         color='p.adjust',
         layout.params = list(layout = 'circle'),
         edge.params = list(min = 0.5),
         cex.params = list(line = 0.1),
         pie='equal',
         shadowtext=TRUE,
         repel=TRUE,
         cex_label_category=0.7,
         showCategory = 50)
goplot(GO_plot,
  showCategory = 17,
  color = "p.adjust",
  layout = "sugiyama",
  geom = "text")
View(data.frame(DI_result))
break
