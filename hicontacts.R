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
           caption = FALSE,
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
gene_interactions<-data.frame(anno)
Gene_symbols<-gene_interactions$symbol
head(Gene_symbols)
q_values<-gene_interactions$qvalue
head(q_values)
enrichme<-data.frame(cbind(Gene_symbols,q_values))
enrichme<-na.omit(enrichme)
enrichme$q_values<-as.numeric(enrichme$q_values)
head(enrichme)
str(enrichme)

#write.csv(enrichme,file = 'test.csv')
#------get gene sets
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
#matrix<-matrix[,-2]
test<-run_pathfindR(
  input = enrichme,
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
#----for individual cluster plots
#visualize_terms(result_df = data,
#                hsa_KEGG = FALSE,
#                pin_name_path = 'mmu_STRING')
#visualize_term_interactions(data,
#                            pin_name_path = 'mmu_STRING',
#                            show_legend = FALSE)
#----------------------
#fishers exact test can be performed to confirm cluster enrichements using a contingency table 
#https://carpentries-incubator.github.io/bioc-rnaseq/06-gene-set-analysis.html
#--------------------------------------------------------
break 
cowplot::plot_grid(
  plotMatrix(hic, use.scores = 'balanced', scale = 'log10', limits = c(-3.5, -1.2), caption = FALSE),
  plotMatrix(hic, use.scores = 'ICE', scale = 'log10', limits = c(-3.5, -1.2), caption = FALSE),
  plotMatrix(hic, use.scores = 'detrended', scale = 'linear', limits = c(-1, 1), cmap = bwrColors(), caption = FALSE), 
  nrow = 1)
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
plotSaddle(hic)
break 
#-------------------------------------------------------------------------
#-----------------------------------------------------------------
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

