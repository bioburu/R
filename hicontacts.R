#https://jserizay.com/OHCA/docs/devel/pages/topological-features.html
library(ShortRead)
library(HiCool)
library(HiContacts)
library(GenomicRanges)
library(HiCExperiment)
library(BiocParallel)
library(HiCcompare)
library(HiCDOC)
library(WGCNA)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(rtracklayer)
library(org.Mm.eg.db)
library(EnsDb.Mmusculus.v79)
library(ChIPpeakAnno)
setwd('/home/em_b/Desktop/hicool/GSM7461656/SRR24843484')
cf<-CoolFile('/home/em_b/Desktop/hicool/GSM7461656/SRR24843484/HiCool_output/matrices/threefiles^mapped-mm10^QNB15W.mcool')
pairs_file<- PairsFile('/home/em_b/Desktop/hicool/GSM7461656/SRR24843484/HiCool_output/pairs/threefiles^mapped-mm10^QNB15W.pairs')
availableResolutions(cf)
hic<-import(cf,
            resolution=1000,
            pairs=pairs_file,
            format = 'mcool',
            focus=c('chr11:115496001-118000000'))
hic<-normalize(hic)
hic <- detrend(hic)
hic
hic<-HiCool::getLoops(hic,resolution=16000)
hic
metadata(hic)$chromosight_args
loops<-topologicalFeatures(hic,'loops')
loops
loops<-loops[loops$score>=0.4&loops$qvalue<=1e-6]
loops
#GenomicInteractions::export.bedpe(loops,'loops.bedpe')
#hic<-zoom(hic,1000)
hic
hic<-getDiamondInsulation(hic,BPPARAM = BiocParallel::bpparam())
hic
hic<-getBorders(hic,weak_threshold = 0.2,strong_threshold = 0.5)
hic
borders<-topologicalFeatures(hic,'borders')
borders
phasing_track<-BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
hic<-getCompartments(hic,
                     chromosomes = 'chr11',
                     genome = phasing_track,
                     resolution=16000)
hic
plotMatrix(hic,
           use.scores = 'balanced',
           scale = 'log10',
           #limits = c(-1, 1),
           #cmap = bwrColors(),
           caption = TRUE,
           loops=loops,
           borders=borders,
           symmetrical=TRUE,
           rasterize=TRUE)
detrend(hic)
loopsdf<-data.frame(topologicalFeatures(hic,'loops'))
loopsdf<-loopsdf[order(loopsdf$qvalue, decreasing=FALSE),]
#topologicalFeatures(hic,'loops')
#seqlevelsStyle(loops)
annoData<-toGRanges(EnsDb.Mmusculus.v79)
annoData
seqlevelsStyle(loops) <- seqlevelsStyle(annoData)
#colnames(loopsdf)
colnames(loopsdf)[1:3]<-c('seqnames','start','end')
Loops <- toGRanges(loopsdf)
anno<-annotatePeakInBatch(Loops,AnnotationData = annoData)
anno <- addGeneIDs(anno, orgAnn="org.Mm.eg.db", 
                   feature_id_type="ensembl_gene_id",
                   IDs2Add=c("symbol"))
anno
enrichme<-data.frame(anno$symbol)
enrichme<-na.omit(enrichme)
enrichme
#write.csv(enrichme,file = 'test.csv')
break 
cowplot::plot_grid(
  plotMatrix(hic, use.scores = 'count', caption = FALSE),
  plotMatrix(hic, use.scores = 'balanced', caption = FALSE),
  plotMatrix(hic, use.scores = 'ICE', caption = FALSE), 
  nrow = 1)
cowplot::plot_grid(
  plotMatrix(hic, use.scores = 'balanced', caption = FALSE),
  plotMatrix(hic, use.scores = 'ICE', caption = FALSE),
  plotMatrix(hic, use.scores = 'expected', caption = FALSE),
  plotMatrix(hic, use.scores = 'detrended', caption = FALSE), 
  nrow = 2)
cowplot::plot_grid(
  plotMatrix(hic, use.scores = 'balanced', scale = 'log10', limits = c(-3.5, -1.2), caption = FALSE),
  plotMatrix(hic, use.scores = 'ICE', scale = 'log10', limits = c(-3.5, -1.2), caption = FALSE),
  plotMatrix(hic, use.scores = 'detrended', scale = 'linear', limits = c(-1, 1), cmap = bwrColors(), caption = FALSE), 
  nrow = 1)
#----detrended linear plots may be best resolution for topological features
plotMatrix(hic, use.scores = 'detrended', scale = 'linear', limits = c(-1, 1), cmap = bwrColors(), caption = FALSE)

plotMatrix(hic,
           loop=loops,
           limits=c(-4,-1.2),
           caption=TRUE)
plotMatrix(hic,
           use.scores = 'detrended',
           scale = 'linear',
           limits = c(-1, 1),
           cmap = bwrColors(),
           caption = FALSE,
           loops=loops)
plotMatrix(hic,
           use.scores = 'detrended',
           scale = 'linear',
           limits = c(-1, 1),
           cmap = bwrColors(),
           caption = FALSE,
           loops=loops,
           borders=borders)
#hic<-autocorrelate(hic)
break 

plotSaddle(hic)
hic
topologicalFeatures(hic,'compartments')
metadata(hic)$eigens
coverage(metadata(hic)$eigens, weight = 'eigen') |> export('hic_eigen.bw')
topologicalFeatures(hic, "compartments") |> export('hic_compartments.gff3')
plotMatrix(hic)
detrend(hic)
hic<-autocorrelate(hic)
hic
p1 <- plotMatrix(hic,
                 use.scores = 'autocorrelated',
                 scale = 'linear',
                 limits = c(-1, 1),
                 caption = FALSE)
p1
eigen <- coverage(metadata(hic)$eigens, weight = 'eigen')[4]
eigen
x<-cumsum(runLength(eigen))
x
y<-runValue(eigen)
y
eigen_df <- data.frame(pos = x, eigen = y)
eigen_df<-eigen_df[,-c(1:2,4:5)]
head(eigen_df)
#-----------------------------------------------------------------------------
str(eigen_df)
barplot(eigen_df$eigen.value)
plot(eigen_df)#,aes=c(x=eigen_df$pos.value,y=eigen_df$eigen.value))
p1
break 
wrap_plots(p1, p2, ncol = 1, heights = c(10, 1))
break 
p2 <- ggplot(eigen_df, aes(x = pos, y = eigen)) + 
  geom_area() + 
  theme_void() + 
  coord_cartesian(expand = TRUE) + 
  labs(x = "Genomic position", y = "Eigenvector value")
p2
wrap_plots(p1, p2, ncol = 1, heights = c(10, 1))
break 
