library(HiCool)
library(HiContacts)
library(GenomicRanges)
library(HiCExperiment)
library(BiocParallel)
library(WGCNA)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(rtracklayer)
library(InteractionSet)
library(EnsDb.Hsapiens.v86)
library(ChIPpeakAnno)
library(org.Hs.eg.db)
library(tibble)
library(biomaRt)
library(dplyr)
library(clusterProfiler)
setwd('/home/deviancedev01/hic/combined_files')
#-------------------------------------------------------------------------------
mcool_path<-file.path('/home/deviancedev01/hic/combined_files/quality_min_30.mcool')
pairs_path<-file.path('/home/deviancedev01/hic/combined_files/tmp/--quality-min=30.valid_idx_filtered.pairs')
log_path<-file.path('/home/deviancedev01/hic/combined_files/SRR26855486_BL6.log/--quality-min=30.hicstuff_20240824033519.log')
cf<-CoolFile(mcool_path,
                 pairs_path,
                 metadata = list(log=log_path),
                 resolution=5000)
cf
#--determine broad focus (~8-12M intervals)
hic<-import(cf,
                focus=c('chr16:6931971-14039342'))
plotMatrix(hic,
           caption=FALSE)
#---- refocus locus closer
hic<-refocus(hic, 'chr16:9931971-13039342')
hic
plotMatrix(hic,
           use.scores = 'balanced',
           limits = c(-4, -1),
           maxDistance=NULL,
           dpi=500,
           caption=FALSE,
           show_grid=TRUE,
           loops=NULL,
           borders=NULL)
#-----------------------------------------------------------------------
hic
hic<-normalize(hic)
hic <- detrend(hic)
hic<-autocorrelate(hic)
hic <- despeckle(hic)
head(metadata(hic))
#---call loops with chromosight
loops <- read.delim('/home/deviancedev01/hic/SRR17223321/out_10kb_loops_small.tsv') 
region1<-loops[,c(1:3)]
head(region1)
colnames(region1)<-c('chr','start','end')
region1<-GRanges(region1)
region2<-loops[,c(4:6)]
head(region2)
colnames(region2)<-c('chr','start','end')
region2<-GRanges(region2)
loops <- GInteractions(region1, region2)
loops
#--- diamond insulation 
hic<-getDiamondInsulation(hic,
                              BPPARAM = BiocParallel::bpparam(),
                              window_size = 5000)
#--- borders
hic<-getBorders(hic,weak_threshold = 0.2,strong_threshold = 0.5)
borders<-topologicalFeatures(hic,'borders')
borders
#--- compartments
phasing_track<-BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
hic<-getCompartments(hic,
                         chromosomes = 'chr18',
                         genome = phasing_track,
                         resolution=5000)
plotMatrix(hic,
           use.scores = 'ICE',
           #limits = c(-4, -1),
           maxDistance=NULL,
           dpi=500,
           caption=TRUE,
           show_grid=TRUE,
           loops=loops,
           borders=borders)
plotMatrix(hic,
           use.scores = 'expected',
           #limits = c(-4, -1),
           maxDistance=NULL,
           dpi=500,
           caption=TRUE,
           show_grid=TRUE,
           loops=loops,
           borders=borders)
plotMatrix(hic,
           use.scores = 'detrended',
           #limits = c(0, -4),
           maxDistance=NULL,
           dpi=500,
           caption=TRUE,
           show_grid=TRUE,
           loops=loops,
           borders=borders)
plotMatrix(hic,
           use.scores = 'autocorrelated',
           #limits = c(0, -4),
           maxDistance=NULL,
           dpi=500,
           caption=TRUE,
           show_grid=TRUE,
           loops=loops,
           borders=borders)
plotMatrix(hic,
           use.scores = 'balanced.despeckled',
           #limits = c(0, -4),
           maxDistance=NULL,
           dpi=500,
           caption=TRUE,
           show_grid=TRUE,
           loops=loops,
           borders=borders)
#---- distance law
ps <- distanceLaw(PairsFile('/home/deviancedev01/hic/combined_files/tmp/--quality-min=30.valid_idx_filtered.pairs'))
ps$sample <- 'Nalm6_sgNUDT21'
plotPs(ps, mapping=ggplot2::aes(x = binned_distance, y = norm_p,color=sample))
#-----annotate interaction sites------------------------------------------------
#---- THIS WILL TAKE AWHILE TO GO THROUGH. MAKE TIME. 
annoData<-toGRanges(EnsDb.Hsapiens.v86)
loopsdf<-read.delim('/home/deviancedev01/hic/SRR17223321/out_10kb_loops_small.tsv') 
loopsdf<-loopsdf[order(loopsdf$qvalue, decreasing=FALSE),]
colnames(loopsdf)[c(1:3)]<-c('seqnames','start','end')
head(loopsdf)
table(loopsdf$seqnames)
to_remove<-c('chrX','chrY','KI270442.1')
loopsdf<-loopsdf[!(loopsdf$seqnames %in% to_remove),]
table(loopsdf$seqnames)
site1_loops <- toGRanges(loopsdf)
site1_loops
site1_anno<-annotatePeakInBatch(site1_loops,
                                AnnotationData = annoData)
site1_anno
site1_anno <- addGeneIDs(site1_anno, orgAnn="org.Hs.eg.db", 
                   feature_id_type="ensembl_gene_id",
                   IDs2Add=c("symbol"))
site1_anno<-data.frame(site1_anno)
site1_anno<-add_column(site1_anno,
                       frag1_genes=site1_anno$symbol,
                       .after = 'end')
loopsdf<-site1_anno
colnames(loopsdf)[c(1:3,5:6)]<-c('seqnames1','start1','end1','width1','strand1')
loopsdf<-loopsdf[,-c(17:26)]
colnames(loopsdf)[c(7:9)]<-c('seqnames','start','end')
Loops2 <- toGRanges(loopsdf)
Loops2
site2_anno<-annotatePeakInBatch(Loops2,
                          AnnotationData = annoData)
site2_anno <- addGeneIDs(site2_anno, orgAnn="org.Hs.eg.db", 
                   feature_id_type="ensembl_gene_id",
                   IDs2Add=c("symbol"))
site2_anno
site2_anno<-data.frame(site2_anno)
#loopsdf<-add_column(loopsdf,
#                    frag2_genes=site2_anno$symbol,
#                    .after = 'end')
merge_names1<-with(site2_anno,paste(seqnames,start,end,sep=':'))
site2_anno<-cbind(merge_names1,site2_anno)
colnames(site2_anno)[1]<-'merge_id'
merge_names2<-with(loopsdf,paste(seqnames,start,end,sep=':'))
head(merge_names2)
head(merge_names1)
loopsdf<-cbind(merge_names2,loopsdf)
colnames(loopsdf)[1]<-'merge_id'
loopsdf<-loopsdf[,-c(8:14)]
colnames(loopsdf)[8:10]<-c('score1','pvalue1','qvalue1')
site2_anno<-site2_anno[,c(1:6,17:19,29)]
colnames(site2_anno) <- paste(colnames(site2_anno),"2",sep="") 
site2_anno<-add_column(site2_anno,
                    frag2_genes=site2_anno$symbol,
                    .after = 'end2')
site2_anno<-site2_anno[,-11]
colnames(site2_anno)[1]<-'merge_id'
loopsdf<-loopsdf[!duplicated(loopsdf$merge_id),]
site2_anno<-site2_anno[!duplicated(site2_anno$merge_id),]
interaction_sites<-merge(loopsdf,site2_anno,by='merge_id')
interaction_sites<-na.omit(interaction_sites)
interaction_sites[,c(2:5,11:14)]
break 
#-------------------------------------------------------------------------------
#--- site1 pathways
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useDataset("hsapiens_gene_ensembl",
                   mart = ensembl)
listAttributes(ensembl)
list1<-interaction_sites$frag1_genes
list1
list2<-interaction_sites$frag2_genes
list2
genes<-c(list1,list2)
genes<-unique(genes)
genes <-getBM(attributes = c('entrezgene_id','external_gene_name'),
              filters = 'external_gene_name',
              values = genes,
              mart = ensembl)
head(genes)
genes<-na.omit(genes)
genes$entrezgene_id
go <- enrichGO(gene = genes$entrezgene_id,
               OrgDb = org.Hs.eg.db,
               ont = "MF")
go
#--- for comparisons of files
#ps2 <- distanceLaw(PairsFile())
#ps2$sample <- ''
#plot_ps<-rbind(ps,ps2)
#plotPs(plot_ps, mapping=ggplot2::aes(x = binned_distance, y = norm_p,color=sample))
#cowplot::plot_grid(
#  plotMatrix(hic,
#             use.scores = 'balanced.despeckled',
#             scale = 'log10',
#             limits = c(-4, -1),
#             loops=loops,
#             borders=borders),
#  plotMatrix(hic2,
#             use.scores = 'balanced.despeckled',
#             scale = 'log10',
#             limits = c(-4, -1),
#             loops=loops,
#             borders=borders),
#  nrow = 1)

