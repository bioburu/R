#library(ShortRead)
library(HiCool)
library(HiContacts)
library(GenomicRanges)
library(HiCExperiment)
library(BiocParallel)
#library(HiCcompare)
#library(HiCDOC)
library(WGCNA)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(rtracklayer)
library(org.Mm.eg.db)
library(regioneR)
library(EnsDb.Mmusculus.v79)
library(ChIPpeakAnno)
library(tibble)
library(pathfindR)
setwd('/home/em_b/Desktop/hicool/SRR26855487_il4_macro_DpnII_HinfI/HiCool')
cf<-CoolFile('/home/em_b/Desktop/hicool/SRR26855487_il4_macro_DpnII_HinfI/HiCool/matrices/SRR26855487_filtered^mapped-mm10^B3IIMF.mcool')
pairs_file<- PairsFile('/home/em_b/Desktop/hicool/SRR26855487_il4_macro_DpnII_HinfI/HiCool/pairs/SRR26855487_filtered^mapped-mm10^B3IIMF.pairs')
availableResolutions(cf)
#-------------------------------------------------------------------------------
hic<-import(cf,
            resolution=1000,
            pairs=pairs_file,
            format = 'mcool',
            focus=c('chr11:115496001-118000000'))
hic<-normalize(hic)
hic <- detrend(hic)
hic<-HiCool::getLoops(hic,resolution=16000)
#metadata(hic)$chromosight_args
loops<-topologicalFeatures(hic,'loops')
#View(data.frame(loops))
loops<-loops[loops$score>=0.4&loops$qvalue<=1e-6]
#GenomicInteractions::export.bedpe(loops,'loops.bedpe')
#hic<-zoom(hic,1000)
hic<-getDiamondInsulation(hic,BPPARAM = BiocParallel::bpparam())
hic<-getBorders(hic,weak_threshold = 0.2,strong_threshold = 0.5)
borders<-topologicalFeatures(hic,'borders')
#View(data.frame(borders))
phasing_track<-BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
hic<-getCompartments(hic,
                     chromosomes = 'chr11',
                     genome = phasing_track,
                     resolution=16000)
hic
#-------------------------------------------------------------------------------
loopsdf<-data.frame(topologicalFeatures(hic,'loops'))
loopsdf<-loopsdf[order(loopsdf$qvalue, decreasing=FALSE),]
annoData<-toGRanges(EnsDb.Mmusculus.v79)
annoData
seqlevelsStyle(loops) <- seqlevelsStyle(annoData)
colnames(loopsdf)[1:3]<-c('seqnames','start','end')
Loops1 <- toGRanges(loopsdf)
anno<-annotatePeakInBatch(Loops1,AnnotationData = annoData)
anno <- addGeneIDs(anno, orgAnn="org.Mm.eg.db", 
                   feature_id_type="ensembl_gene_id",
                   IDs2Add=c("symbol"))
#View(data.frame(anno))
Loops1<-data.frame(Loops1)
Loops1<-add_column(Loops1, frag1_genes=anno$symbol, .after = 'seqnames')
#-------------------------------------------------------------------------------
gc()
hic<-autocorrelate(hic)
hic
summary(scores(hic, 'autocorrelated'))
ps <- distanceLaw(hic)
ps
plotPs(ps, ggplot2::aes(x = binned_distance, y = norm_p))
cat('plot shows contact frequencies are more abundant at shorter distances')
plotPsSlope(ps, ggplot2::aes(x = binned_distance, y = slope))
cat('this plot is useful if comparing two conditions and overlaying')
gc()
cowplot::plot_grid(
  #  plotMatrix(hic, 
  #             use.scores = 'count',
  #             caption = FALSE,
  #             loops=loops,
  #             borders=borders,
  #             symmetrical=TRUE,
  #             rasterize=TRUE),
  plotMatrix(hic,
             use.scores = 'balanced',
             caption = FALSE,
             loops=loops,
             borders=borders,
             symmetrical=TRUE,
             rasterize=TRUE),
  plotMatrix(hic,
             use.scores = 'ICE',
             caption = FALSE,
             loops=loops,
             borders=borders,
             symmetrical=TRUE,
             rasterize=TRUE),
  plotMatrix(hic,
             use.scores = 'expected',
             caption = FALSE,
             loops=loops,
             borders=borders,
             symmetrical=TRUE,
             rasterize=TRUE),
  plotMatrix(hic,
             use.scores = 'detrended',
             caption = FALSE,
             loops=loops,
             borders=borders,
             symmetrical=TRUE,
             rasterize=TRUE),
  nrow = 2)
cowplot::plot_grid(
  plotMatrix(hic,
             use.scores='balanced',
             limits=c(-4,-1),
             maxDistance=2000000,
             caption=FALSE),
  plotMatrix(hic,
             use.scores='ICE',
             limits=c(-4,-1),
             maxDistance=2000000,
             caption=FALSE),
  plotMatrix(hic,
             use.scores='expected',
             limits=c(-4,-1),
             maxDistance=2000000,
             caption=FALSE),
  plotMatrix(hic,
             use.scores='detrended',
             limits=c(-4,-1),
             maxDistance=2000000,
             caption=FALSE),
  nrow = 2)
plotMatrix(hic, use.scores = 'autocorrelated', limits = c(-1, 1), scale = 'linear')
#----check resolutions here 
hic_despeckled<-zoom(hic, 8000)
hic_despeckled <- despeckle(hic_despeckled)
plotMatrix(hic_despeckled, use.scores = 'balanced.despeckled', scale = 'log10', limits = c(-4, -1))
#--this is the best plot so far
plotMatrix(hic, use.scores = 'detrended', scale = 'linear', limits = c(-1, 1), dpi = 240)

gene_interactions<-data.frame(anno)
Gene_symbols<-gene_interactions$symbol
q_values<-gene_interactions$qvalue
gene_annotations<-data.frame(cbind(Gene_symbols,q_values))
gene_annotations<-na.omit(gene_annotations)
gene_annotations$q_values<-as.numeric(gene_annotations$q_values)
head(gene_annotations)
str(gene_annotations)
#write.csv(gene_annotations,file = 'gene_annotations.csv')
gene_list <- fetch_gene_set(
  gene_sets = "mmu_KEGG",
  min_gset_size = 10,
  max_gset_size = 300)
#------------------------------------------------------------------------
mmu_kegg_gsets <- gene_list[[1]]
head(mmu_kegg_gsets)
#------------------------------------------------------------
mmu_kegg_descriptions <- gene_list[[2]]
head(mmu_kegg_descriptions)
test<-run_pathfindR(
  input = gene_annotations,
  gene_sets = 'mmu_KEGG',
  convert2alias = FALSE,
  custom_genes = mmu_KEGG_gsets,
  custom_descriptions = mmu_KEGG_descriptions,
  pin_name_path = 'mmu_STRING'
)
data<-cluster_enriched_terms(test)
term_gene_graph(test)
term_gene_graph(test,num_terms = c(5),use_description = TRUE)
enrichment_chart(test, plot_by_cluster = TRUE,top_terms = 20,num_bubbles = 8)
data<-data[order(data$highest_p), ]
#write.csv(data,file = 'bulk.tiss_mmu_KEGG_results.csv')
plotSaddle(hic)
aggr_centros<-aggregate(hic,
                        targets=loops,
                        flankingBins=5)
plotMatrix(aggr_centros,
           #use.scores='detrended',
           #limits=c(-1,-1),
           #scale='linear',
           cmap=bgrColors())
break 

