library(Gviz)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm39)
library(ChIPseeker)
library(BRGenomics)
library(GenomicRanges)
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
txdb<-TxDb.Mmusculus.UCSC.mm39.knownGene
setwd('/home/em_b/work_stuff/chipseq/final/peaks_q_0.05')
peak <- readPeakFile('D1_J1_peaks.narrowPeak')
table(seqnames(peak))
peak
peak<-tidyChromosomes(peak,
                      keep.X=FALSE,
                      keep.Y=FALSE,
                      keep.M = FALSE,
                      keep.nonstandard = FALSE,
                      genome = 'mm39')
table(seqnames(peak))
peak<-peak[order(peak$V5, decreasing=TRUE),]
peak<-data.frame(peak)
summary(peak$V5)
cat('q_value set to 0.01 on macs. cutoff for 10*-log10qvalue at V5 is 0.001==30')
peak<-subset(peak, V5 >30)
summary(peak$V5)
GR <- GRanges(seqnames=peak$seqnames,
              IRanges(peak$start,
                      peak$end))
peakAnno <- annotatePeak(GR, tssRegion=c(-1000, 1000), 
                         TxDb=txdb, 
                         annoDb="org.Mm.eg.db")
anno_df<-data.frame(peakAnno)
reactome <- enrichPathway(anno_df$geneId,
                          organism = 'mouse')
reactome
#--if subsetting use
#pathways<-c('R-MMU-6809371')
#reactome@result<-reactome@result[reactome@result$ID%in%pathways,]

dotplot(reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=15,
        label_format=50,
        title='')
#View(data.frame(reactome))
edox <- setReadable(reactome, 'org.Mm.eg.db', 'ENTREZID')
cnetplot(edox,
         #foldChange=gene)#,
         node_label='all',
         colorEdge=TRUE)
cnetplot(edox,
         #foldChange=gene,
         node_label='all',
         colorEdge=TRUE,
         circular=TRUE)
go <- enrichGO(gene = anno_df$geneId,
               OrgDb = org.Mm.eg.db,
               ont = "BP")
#View(data.frame(go))
#---subset if necessary
#pathways<-c('GO:0007416','GO:0048167','GO:0035249','GO:0021953','GO:1905606',
#            'GO:0099174','GO:0099643','GO:0097106','GO:0046928','GO:0021549',
#            'GO:0045666','GO:0034329') 
#go@result = go@result[go@result$ID %in% pathways,]
#go@result
GO_plot <- pairwise_termsim(go)
#View(data.frame(GO_plot))
emapplot(GO_plot,
         color='p.adjust',
         layout.params = list(layout = 'dh'),
         edge.params = list(min = 0.5),
         cex.params = list(line = 0.1),
         pie='circle',
         shadowtext=TRUE,
         repel=TRUE,
         cex_label_category=1,
         showCategory = 50)
goplot(GO_plot,
       showCategory = 17,
       color = "p.adjust",
       layout = "sugiyama",
       geom = "text")
summary(peak)
break 
