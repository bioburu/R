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
mcool_path<-file.path('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/matrices/SRR26855486_BL6.mcool')
pairs_path<-file.path('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/pairs/SRR26855486_BL6.pairs')
log_path<-file.path('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/logs/SRR26855486_BL6.log')
bl6_cf<-CoolFile(mcool_path,
             pairs_path,
             metadata = list(log=log_path),
             resolution=16000)
bl6_cf
bl6_hic<-import(bl6_cf,
            focus=c('chr18:34368001-38736000'))
bl6_hic
gc()
#-------------------------------------------------------------------------------
mcool_path<-file.path('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/matrices/SRR26855487_BALBc.mcool')
pairs_path<-file.path('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/pairs/SRR26855487_BALBc.pairs')
log_path<-file.path('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/logs/SRR26855487_BALBc.log')
balbc_cf<-CoolFile(mcool_path,
                 pairs_path,
                 metadata = list(log=log_path),
                 resolution=16000)
balbc_cf
balbc_hic<-import(balbc_cf,
                focus=c('chr18:34368001-38736000'))
balbc_hic
gc()
p1<-plotMatrix(bl6_hic)
p2<-plotMatrix(balbc_hic)
p1+p2
#-------------------------------------------------------------------------------
setwd('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/BL6_loops')
bl6_hic<-HiCool::getLoops(bl6_hic,
                          resolution=16000,
                          ncores = 15)
bl6_hic
setwd('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/BALBc_loops')
balbc_hic<-HiCool::getLoops(balbc_hic,
                            resolution=16000,
                            ncores = 15)
balbc_hic
metadata(bl6_hic)$chromosight_args
bl6_loops<-topologicalFeatures(bl6_hic,'loops')
balbc_loops<-topologicalFeatures(balbc_hic,'loops')
bl6_loops<-bl6_loops[bl6_loops$score>=0.4&bl6_loops$qvalue<=1e-6]
balbc_loops<-balbc_loops[balbc_loops$score>=0.4&balbc_loops$qvalue<=1e-6]
setwd('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output')
#GenomicInteractions::export.bedpe(bl6_loops,'loops.bedpe')
bl6_hic<-getDiamondInsulation(bl6_hic,
                              BPPARAM = BiocParallel::bpparam(),
                              window_size = 16000)
balbc_hic<-getDiamondInsulation(balbc_hic,
                                BPPARAM = BiocParallel::bpparam(),
                                window_size = 16000)
bl6_hic<-getBorders(bl6_hic,weak_threshold = 0.2,strong_threshold = 0.5)
balbc_hic<-getBorders(balbc_hic,weak_threshold = 0.2,strong_threshold = 0.5)
bl6_borders<-topologicalFeatures(bl6_hic,'borders')
balbc_borders<-topologicalFeatures(balbc_hic,'borders')
phasing_track<-BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
bl6_hic<-getCompartments(bl6_hic,
                     chromosomes = 'chr18',
                     genome = phasing_track,
                     resolution=16000)
balbc_hic<-getCompartments(balbc_hic,
                         chromosomes = 'chr18',
                         genome = phasing_track,
                         resolution=16000)
bl6_hic
balbc_hic
gc()
#-------------------------------------------------------------------------------
#BL6_ps <- distanceLaw(PairsFile('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/pairs/SRR26855486_BL6.pairs'))
#BL6_ps$sample <- 'BL6_il4_macrophages'
#gc()
#BALBc_ps <- distanceLaw(PairsFile('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/pairs/SRR26855487_BALBc.pairs'))
#BALBc_ps$sample <- 'BALBc_il4_macrophages'
#gc()
#plot_ps<-rbind(BL6_ps,BALBc_ps)
#plot_ps
#plotPs(plot_ps, mapping=ggplot2::aes(x = binned_distance, y = norm_p,color=sample))
gc()
#-------------------------------------------------------------------------------
bl6_hic<-normalize(bl6_hic)
balbc_hic<-normalize(balbc_hic)
bl6_hic <- detrend(bl6_hic)
balbc_hic <- detrend(balbc_hic)
bl6_hic<-autocorrelate(bl6_hic)
balbc_hic<-autocorrelate(balbc_hic)
bl6_hic <- despeckle(bl6_hic)
balbc_hic <- despeckle(balbc_hic)
gc()
cowplot::plot_grid(
  plotMatrix(bl6_hic,
             use.scores = 'balanced.despeckled',
             scale = 'log10',
             limits = c(-4, -1),
             loops=bl6_loops,
             borders=bl6_borders),
  plotMatrix(balbc_hic,
             use.scores = 'balanced.despeckled',
             scale = 'log10',
             limits = c(-4, -1),
             loops=balbc_loops,
             borders=balbc_borders),
  nrow = 1)
annoData<-toGRanges(EnsDb.Mmusculus.v79)
annoData
#-------------------------------------------------------------------------------
bl6_loopsdf<-data.frame(topologicalFeatures(bl6_hic,'loops'))
bl6_loopsdf<-bl6_loopsdf[order(bl6_loopsdf$qvalue, decreasing=FALSE),]
colnames(bl6_loopsdf)[c(1:3,5)]<-c('seqnames','start','end','strand')
str(bl6_loopsdf)
summary(bl6_loopsdf$seqnames)
cat('remove chrX, chrY, and chrX_GL456233_random')
to_remove<-c('chrX','chrY','chrX_GL456233_random')
bl6_loopsdf<-bl6_loopsdf[!(bl6_loopsdf$seqnames %in% to_remove),]
summary(bl6_loopsdf$seqnames)
summary(bl6_loopsdf$seqnames2)
Loops1 <- toGRanges(bl6_loopsdf)
Loops1
anno<-annotatePeakInBatch(Loops1,AnnotationData = annoData)
anno <- addGeneIDs(anno, orgAnn="org.Mm.eg.db", 
                   feature_id_type="ensembl_gene_id",
                   IDs2Add=c("symbol"))
Loops1<-data.frame(Loops1)
Loops1<-add_column(Loops1, frag1_genes=anno$symbol, .after = 'end')
bl6_loopsdf<-Loops1
summary(bl6_loopsdf$seqnames1)
summary(bl6_loopsdf$seqnames2)
colnames(bl6_loopsdf)[c(1:3,5:6)]<-c('seqnames1','start1','end1','width1','strand1')
bl6_loopsdf<-bl6_loopsdf[,-7]
colnames(bl6_loopsdf)[c(7:9,11)]<-c('seqnames','start','end','strand')
Loops1 <- toGRanges(bl6_loopsdf)
Loops1
anno<-annotatePeakInBatch(Loops1,AnnotationData = annoData)
anno <- addGeneIDs(anno, orgAnn="org.Mm.eg.db", 
                   feature_id_type="ensembl_gene_id",
                   IDs2Add=c("symbol"))
Loops2<-data.frame(Loops1)
cat('merge anno to Loops2')
anno<-data.frame(anno)
names_test<-with(Loops2,paste(seqnames,start,end,sep=':'))
names_test
Loops2<-cbind(names_test,Loops2)
colnames(Loops2)[1]<-'merge_id'
names<-with(anno,paste(seqnames,start,end,sep=':'))
summary(names_test)
summary(names)
anno<-cbind(names,anno)
colnames(anno)[1]<-'merge_id'
anno<-anno[,-c(2:27)]
Loops2_dedup<-Loops2[!duplicated(Loops2$merge_id),]
anno_dedup<-anno[!duplicated(anno$merge_id),]
df<-merge(Loops2_dedup,anno_dedup,by='merge_id')
df<-df[,-c(1,13)]
colnames(df)[c(1:5)]<-c('seqnames2','start2','end2','width2','strand2')
df<-add_column(df, frag2_genes=df$symbol, .after = 'end2')
df<-df[,-18]
df<-df[order(df$qvalue, decreasing=FALSE),]
df<-df[,-c(6,12)]
View(df)
gc()
frag1_genes<-df$frag1_genes
frag1_genes
frag1_q_values<-df$qvalue
frag1_q_values
frag2_genes<-df$frag2_genes
frag2_genes
frag2_q_values<-df$qvalue
frag2_q_values
Gene_symbols<-c(frag1_genes,frag2_genes)
Gene_symbols
q_values<-c(frag1_q_values,frag2_q_values)
q_values
gene_annotations<-data.frame(cbind(Gene_symbols,q_values))
str(gene_annotations)
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
#write.csv(data,file = '.csv')
plotSaddle(bl6_hic)
plotSaddle(balbc_hic)
cat('redo pathway analysis on other seqnames if needed')
break 


