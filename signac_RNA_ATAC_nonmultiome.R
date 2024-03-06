library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(caTools)
#library(car)
library(caret)
library(InformationValue)
library(pROC)
library(ROCR)
fragment<-file.path('/home/em_b/Desktop/FCCC/GSE218223_scATACseq_mouseCD8/test/fragments.sorted.bed.gz')
total_counts <- CountFragments(fragment)
head(total_counts)
summary(total_counts$frequency_count)
cutoff <- 1000 
barcodes <- total_counts[total_counts$frequency_count > cutoff, ]$CB
head(barcodes)
frags <- CreateFragmentObject(path = fragment, cells = barcodes)
frags
peaks <- CallPeaks(frags,macs2.path = '/home/em_b/anaconda3/bin/macs3')
peaks
counts <- FeatureMatrix(fragments = frags, features = peaks, cells = barcodes)
head(counts)
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = frags,
  min.cells = 10,
  min.features = 200,
  genome = 'mm10')
atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'ATAC',
  project = 'scATAC')
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
annotations
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
annotations
Annotation(atac) <- annotations
#---------------------------------------------------------------------------
atac[['ATAC']]
granges(atac)
head(atac)
summary(atac$nCount_ATAC)
VlnPlot(atac,features=c('nCount_ATAC','nFeature_ATAC'),
        pt.size=0.1,
        ncol=2,
        cols = 'skyblue')
atac<-subset(atac,subset=nCount_ATAC>1000)
VlnPlot(atac,features=c('nCount_ATAC','nFeature_ATAC'),
        pt.size=0.1,
        ncol=2,
        cols = 'skyblue')
cat('subset out all TSS scores below 4')
atac <- NucleosomeSignal(object = atac)
atac <- TSSEnrichment(atac)
DensityScatter(atac, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
head(atac)
VlnPlot(object = atac,
        features = c('nCount_ATAC',
                     'nFeature_ATAC',
                     'nucleosome_signal',
                     'nucleosome_percentile',
                     'TSS.enrichment','TSS.percentile'),
        pt.size = 0.1,
        ncol = 3,
        cols = 'skyblue')
atac$high.tss <- ifelse(atac$TSS.enrichment > 4, 'High', 'Low')
table(atac$high.tss)
#TSSPlot(data, group.by = 'high.tss')# + NoLegend()
atac<-subset(atac,subset=TSS.enrichment>4)
VlnPlot(object = atac,
        features = c('nCount_ATAC',
                     'nFeature_ATAC',
                     'nucleosome_signal',
                     'nucleosome_percentile',
                     'TSS.enrichment','TSS.percentile'),
        pt.size = 0.1,
        ncol = 3,
        cols = 'skyblue')
#-----------------------------------------------------------------------------------
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
top500.atac<-as.data.frame(FindTopFeatures(atac[['ATAC']][])[1:500,])
View(top500.atac)
atac <- RunSVD(atac) 
DepthCor(atac)
cat('remove highly correlated dimensions')
atac <- RunUMAP(object = atac, reduction = 'lsi', dims = 2:30)
atac <- FindNeighbors(object = atac, reduction = 'lsi', dims =2:30)
atac <- FindClusters(object = atac, verbose = FALSE, algorithm = 3)
cat('cleaned ATACseq data is now ok for gene activity analysis')
gene.activities <- GeneActivity(atac)
atac[['correlated_gene_activity']] <- CreateAssayObject(counts = gene.activities)
cat('transfer to gene activity layer')
DefaultAssay(atac)<-'correlated_gene_activity'
atac <- NormalizeData(
  object = atac,
  assay = 'correlated_gene_activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac$nCount_ATAC)
)
atac<-FindVariableFeatures(atac)
atac<-ScaleData(atac)
top1000.atac <- head(VariableFeatures(atac), 1000)
top1000.atac
#-------------------------------------------------------------------------------
features_path <- '/home/em_b/Desktop/FCCC/GSE218223_scATACseq_mouseCD8/geo_results/GSM6745551_IL1218_Neonate_features.tsv.gz'
barcodes_path <- '/home/em_b/Desktop/FCCC/GSE218223_scATACseq_mouseCD8/geo_results/GSM6745551_IL1218_Neonate_barcodes.tsv.gz'
matrix_path <- '/home/em_b/Desktop/FCCC/GSE218223_scATACseq_mouseCD8/geo_results/GSM6745551_IL1218_Neonate_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
rna <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'rna')
summary(rna@active.ident)
gc()
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^mt-")
VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
rna <- subset(rna, subset = nFeature_RNA > 200 & nFeature_RNA < 11000)# & percent.mt <15)
VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
rna <- NormalizeData(rna)
gc()
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
VariableFeatures(rna)
top1000.rna <- head(VariableFeatures(rna), 1000)
top1000.rna
plot1 <- VariableFeaturePlot(rna)
plot1
all.genes <- rownames(rna)
all.genes
rna <- ScaleData(rna, features = all.genes)
gc()
dim(rna)
rna <- RunPCA(rna, features = VariableFeatures(object = rna))
rna <- RunUMAP(rna, dims = 1:30)
DimHeatmap(rna, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(rna)
rna <- FindNeighbors(rna, dims = 1:30)
rna <- FindClusters(rna)
DimPlot(object = rna, label = TRUE,pt.size = 1)
VlnPlot(rna, features = c('Ptprc','Gata3','Cd3d','Cd4','Cd8a','Il12rb1','Cd19','Cd22','Ncam1','Mki67'),cols = c())
FeaturePlot(
  object = rna,
  features = c('Ptprc','Gata3','Cd3d','Cd4','Cd8a','Il12rb1','Cd19','Cd22','Ncam1','Mki67'),
  pt.size = 1,
  max.cutoff = 'q95',
  ncol = 3,
  cols = c('grey','red'))
Idents(rna)
levels(rna@meta.data$seurat_clusters)<-c('IFNhi1','IFNhi2','cMem3','Cytotoxic','KI67hi','IFNGlo1',
                                         'cMem1','cMem2','rMem','IFNlo2','IFNlo3','GZMA')
rna@meta.data$seurat_clusters
head(rna)
break 
new.cluster.ids <- c('IFNhi1','IFNhi2','cMem3','Cytotoxic','KI67hi','IFNGlo1',
                     'cMem1','cMem2','rMem','IFNlo2','IFNlo3','GZMA')
names(new.cluster.ids) <- levels(rna)
rna <- RenameIdents(rna, new.cluster.ids)
DimPlot(rna, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
head(atac)
p1 <- DimPlot(rna, label = TRUE,pt.size = 1) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(atac, group.by = "orig.ident", label = FALSE,pt.size = 1) + NoLegend() + ggtitle("ATAC")
p1 + p2
# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = VariableFeatures(object = rna),
                                        reference.assay = "RNA", query.assay = "correlated_gene_activity", reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$seurat_clusters,
                                     weight.reduction = atac[["lsi"]], dims = 2:30)
atac <- AddMetaData(atac, metadata = celltype.predictions)
p3<-DimPlot(atac, group.by = "predicted.id", label = TRUE,pt.size = 1) + ggtitle('RNA+ATAC')
p1+p2+p3
break 
FeaturePlot(atac, features = c('Cd2','Cd27','Cd28','Cd38','Cd40','Cd40lg','Sell','Cd69','Cd247','Cd274','Pdcd1','CCR4','CCR6','CCR7','CXCR3','CCL22'))
FeaturePlot(
  object = atac,
  features = c('Ifng','Nfkb2','Tnf','Lamp1','Lamp2','Gzmb','Gzma','Il1a'),
  pt.size = 1,
  max.cutoff = 'q95',
  cols = c('grey','red'))
DimPlot(object = rna, label = TRUE,pt.size = 2)
x<-FindMarkers(data, ident.1 = '1', 
               features = c(top1000),test.use='negbinom')
VlnPlot(data, features = c(row.names(x)[1:12]),cols = c())
VlnPlot(data, features = c(row.names(x)[13:24]),cols = c())
VlnPlot(data, features = c(row.names(x)[25:36]),cols = c())
FeaturePlot(
  object = atac,
  features = c('Nek6','Ksr1','C3','Tmem108','Arpp21','Reep1','Fes','Il18','Mrc1'),
  pt.size = 1,
  max.cutoff = 'q95',
  ncol = 3,
  cols = c('grey','red'))
FeaturePlot(atac, features = c('Ptprc','Cd8a','Ifng','Il4','Cd3d','Il12rb1','Ttc39c','Gas7','Map6'))
break
setwd('/home/em_b/Desktop/FCCC/GSE218223_scATACseq_mouseCD8/my_results')
saveRDS(data, file = "my.rds")
data<-readRDS('my.rds')
