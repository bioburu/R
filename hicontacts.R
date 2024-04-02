library(HiCool)
library(HiContacts)
library(GenomicRanges)
library(HiCExperiment)
library(BiocParallel)
library(HiCcompare)
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
setwd('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output')
#-------------------------------------------------------------------------------
BL6_cf<-CoolFile('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/matrices/SRR26855486_BL6.mcool')
BALBc_cf<-CoolFile('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/matrices/SRR26855487_BALBc.mcool')
availableResolutions(BL6_cf)
availableResolutions(BALBc_cf)
#-------------------------------------------------------------------------------
BL6_pairs<- PairsFile('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/pairs/SRR26855486_BL6.pairs')
BALBc_pairs<- PairsFile('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/pairs/SRR26855487_BALBc.pairs')
#-------------------------------------------------------------------------------
BL6_hic<-import(BL6_cf,
            resolution=1000,
            pairs=BL6_pairs,
            format = 'mcool',
            focus=c('chr12:24000000-28000000'))
BALBc_hic<-import(BALBc_cf,
                resolution=1000,
                pairs=BALBc_pairs,
                format = 'mcool',
                focus=c('chr12:24000000-28000000'))
gc()
#cowplot::plot_grid(
#  plotMatrix(BL6_hic,
#             use.scores = 'balanced',
#             #scale = 'log10',
##             limits = c(-4, -1),
#             loops=BL6_loops,
#             borders=BL6_borders),
#  plotMatrix(BALBc_hic,
#             use.scores = 'balanced',
#             #scale = 'log10',
#             limits = c(-4, -1),
#             loops=BALBc_loops,
#             borders=BALBc_borders),
#  nrow = 1)
#-------------------------------------------------------------------------------
BL6_hic<-normalize(BL6_hic)
BALBc_hic<-normalize(BALBc_hic)
#-------------------------------------------------------------------------------
BL6_hic <- detrend(BL6_hic)
BALBc_hic <- detrend(BALBc_hic)
#-------------------------------------------------------------------------------
BL6_hic<-autocorrelate(BL6_hic)
BALBc_hic<-autocorrelate(BALBc_hic)
summary(scores(BL6_hic, 'autocorrelated'))
summary(scores(BALBc_hic, 'autocorrelated'))
#-------------------------------------------------------------------------------
BL6_hic <- despeckle(BL6_hic)
BALBc_hic <- despeckle(BALBc_hic)
gc()
#-------------------------------------------------------------------------------
setwd('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/BL6_loops')
BL6_hic<-HiCool::getLoops(BL6_hic,resolution=16000)
setwd('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/BALBc_loops')
BALBc_hic<-HiCool::getLoops(BALBc_hic,resolution=16000)
#metadata(hic)$chromosight_args
BL6_loops<-topologicalFeatures(BL6_hic,'loops')
BALBc_loops<-topologicalFeatures(BALBc_hic,'loops')
#View(data.frame(loops))
BL6_loops<-BL6_loops[BL6_loops$score>=0.4&BL6_loops$qvalue<=1e-6]
BALBc_loops<-BALBc_loops[BALBc_loops$score>=0.4&BALBc_loops$qvalue<=1e-6]
setwd('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output')
#GenomicInteractions::export.bedpe(BL6_loops,'loops.bedpe')
#hic<-zoom(hic,1000)

#-------------------------------------------------------------------------------
BL6_hic<-getDiamondInsulation(BL6_hic,BPPARAM = BiocParallel::bpparam())
BALBc_hic<-getDiamondInsulation(BALBc_hic,BPPARAM = BiocParallel::bpparam())
BL6_hic<-getBorders(BL6_hic,weak_threshold = 0.2,strong_threshold = 0.5)
BALBc_hic<-getBorders(BALBc_hic,weak_threshold = 0.2,strong_threshold = 0.5)
BL6_borders<-topologicalFeatures(BL6_hic,'borders')
BALBc_borders<-topologicalFeatures(BALBc_hic,'borders')
#View(data.frame(borders))
#-------------------------------------------------------------------------------
phasing_track<-BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
BL6_hic<-getCompartments(BL6_hic,
                     chromosomes = 'chr12',
                     genome = phasing_track,
                     resolution=16000)
BALBc_hic<-getCompartments(BALBc_hic,
                         chromosomes = 'chr12',
                         genome = phasing_track,
                         resolution=16000)
BL6_hic
BALBc_hic
gc()
#-------------------------------------------------------------------------------
BL6_ps <- distanceLaw(BL6_pairs)
BL6_ps$sample <- 'BL6_il4_macrophages'
#bl6_ps<-mutate(BL6_ps,sample='BL6_il4_macrophages')
gc()
BALBc_ps <- distanceLaw(BALBc_pairs)
BALBc_ps$sample <- 'BALBc_il4_macrophages'
#balbc_ps<-mutate(BALBc_ps,sample='BALBc_il4_macrophages')
gc()
plot_ps<-rbind(BL6_ps,BALBc_ps)
plot_ps
plotPs(plot_ps, mapping=ggplot2::aes(x = binned_distance, y = norm_p,color=sample))
#plotPsSlope(plot_ps, mapping=ggplot2::aes(x = binned_distance, y = slope,color=sample))
gc()
#-------------------------------------------------------------------------------
cowplot::plot_grid(
  plotMatrix(BL6_hic,
             use.scores = 'balanced.despeckled',
             scale = 'log10',
             limits = c(-4, -1),
             loops=BL6_loops,
             borders=BL6_borders),
  plotMatrix(BALBc_hic,
             use.scores = 'balanced.despeckled',
             scale = 'log10',
             limits = c(-4, -1),
             loops=BALBc_loops,
             borders=BALBc_borders),
  nrow = 1)
break 
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
break 
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
#---try to get into form for use in multiHicCompare?????
library(purrr)
library(HiCExperiment)
library(multiHiCcompare)
library(tidyverse)
BL6_cf<-CoolFile('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/matrices/SRR26855486_BL6.mcool')
BL6_pairs<- PairsFile('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/pairs/SRR26855486_BL6.pairs')
BL6_hic<-import(BL6_cf,
                resolution=1000,
                pairs=BL6_pairs,
                format = 'mcool',
                focus=c('chr12:24000000-28000000'))
test<-data.frame(BL6_hic@scores@listData[["count"]])
BL6_hic
summary(test)


