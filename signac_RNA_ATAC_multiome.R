library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(caTools)
#library(car)
library(caret)
#library(InformationValue)
library(pROC)
library(ROCR)
fragment<-file.path('/home/em_b/Desktop/snRNA+ATAC_huPBMC/GSM7830528_Sample_0_atac_fragments.tsv.gz')
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
colnames(counts)
row.names(counts)
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = frags,
  min.cells = 10,
  min.features = 200,
  genome = 'hg38')
atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'ATAC',
  project = 'scATAC')
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
annotations
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
annotations
Annotation(atac) <- annotations
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
atac <- RunTFIDF(atac)
atac
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunSVD(atac) 
DepthCor(atac)
cat('remove highly correlated dimensions')
atac <- RunUMAP(object = atac, reduction = 'lsi', dims = 2:30)
atac <- FindNeighbors(object = atac, reduction = 'lsi', dims =2:30)
atac <- FindClusters(object = atac, verbose = FALSE, algorithm = 3,resolution = 0.1)
Idents(atac)
atac
DimPlot(atac,label = TRUE)
top500.atac<-FindTopFeatures(atac)[1:500]
row.names(top500.atac)
cat('cleaned ATACseq data is now ok for gene activity analysis')
gene.activities <- GeneActivity(atac)
atac[['correlated_gene_activity']] <- CreateAssayObject(counts = gene.activities)
cat('transfer to gene activity layer')
DefaultAssay(atac)<-'correlated_gene_activity'
atac <- NormalizeData(
  object = atac,
  assay = 'correlated_gene_activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac$nCount_ATAC))
atac<-FindVariableFeatures(atac)
atac<-ScaleData(atac)
atac
top1000.atac <- head(VariableFeatures(atac), 1000)
top1000.atac
VariableFeatures(atac)
DefaultAssay(atac)<-'ATAC'
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
features_path <- '/home/em_b/Desktop/snRNA+ATAC_huPBMC/GSM7830515_Sample_0_features.tsv.gz'
barcodes_path <- '/home/em_b/Desktop/snRNA+ATAC_huPBMC/GSM7830515_Sample_0_barcodes.tsv.gz'
matrix_path <- '/home/em_b/Desktop/snRNA+ATAC_huPBMC/GSM7830515_Sample_0_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
colnames(matrix)
row.names(matrix)
rna <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'rna')
summary(rna@active.ident)
gc()
#rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")
VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
rna <- subset(rna, subset = nFeature_RNA > 200 & nFeature_RNA < 11000)# & percent.mt <5)
VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
rna <- NormalizeData(rna)
gc()
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
VariableFeatures(rna)
rna
rna_df<-data.frame(rna[['RNA']]$counts)
atac_df<-data.frame(rna[['RNA']]$counts)
top5000.rna <- head(VariableFeatures(rna), 5000)
top5000.rna
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
rna <- FindClusters(rna,resolution = 0.1)
DimPlot(object = rna, label = TRUE,pt.size = 1)
VlnPlot(rna, features = c('PTPRC','CD3D','CD4','CD8A','IFNG','CD19','CD22','CD86','ITGAX'),cols = c())
FeaturePlot(
  object = rna,
  features = c('PTPRC','CD3D','CD4','CD8A','IFNG','CD19','CD22','CD86','ITGAX'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  cols = c('grey','red'))
DimPlot(object = rna, label = TRUE,pt.size = 1)
Idents(rna)
levels(rna@meta.data$seurat_clusters)<-c('cytoxic_t_cells','t_helper_cells','mono/macro','b_cells')
rna@meta.data$seurat_clusters
head(rna)
new.cluster.ids <- c('cytotoxic_t_cells','t_helper_cells','mono/macro','b_cells')
names(new.cluster.ids) <- levels(rna)
rna <- RenameIdents(rna, new.cluster.ids)
DimPlot(rna, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
head(atac)
transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = VariableFeatures(object = rna),
                                        reference.assay = "RNA", query.assay = "correlated_gene_activity", reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$seurat_clusters,
                                     weight.reduction = atac[["lsi"]], dims = 1:30)
atac <- AddMetaData(atac, metadata = celltype.predictions)
head(atac)
head(Idents(atac))
head(Idents(rna))
rna$seurat_clusters
new.cluster.ids <- c('cytotoxic_t_cells','t_helper_cells','mono/macro','b_cells')
names(new.cluster.ids) <- levels(atac)
names(new.cluster.ids)
new.cluster.ids
levels(atac)
Idents(atac)
atac <- RenameIdents(atac, new.cluster.ids)
head(Idents(atac))
head(Idents(rna))
DimPlot(object = atac, label = TRUE,pt.size = 1)
rna_umap <- DimPlot(rna, label = TRUE,pt.size = 1,label.size = 8)+ ggtitle("gene_transcripts")
atac_features<-FeaturePlot(atac,features = c('PTPRC','CD3D','CD4','ICD8A','IFNG','CD19',
                              'CD22','CD86','ITGAX'),cols = c('grey','red'),pt.size = 2)
atac_umap<-DimPlot(object = atac, label = TRUE,pt.size = 1,label.size = 8)
rna_umap+atac_umap
#-------------------------------------------------------------------
peaks <- FindMarkers(
  object = atac,
  ident.1 = 'mono/macro',
  test.use = 'negbinom',
  latent.vars = 'nCount_ATAC',
  logfc.threshold=1.5,
  min.pct=0.2)
peaks
FeaturePlot(atac,features = c(row.names(peaks)[1:12]))
VlnPlot(
  object = atac,
  features = rownames(peaks)[1:12],
  pt.size = 0.1,
  idents = c())
ClosestFeature(atac,c('chr4-8409026-8410072'))
p1<-CoveragePlot(
  object = atac,
  region = 'chr4-8409026-8410072',
  extend.upstream = 40000,
  extend.downstream = 20000)
p2 <- TilePlot(
  object = atac,
  region = 'chr4-8409026-8410072')
p3<-ExpressionPlot(
  object = rna,
  features = c("ACOX3"))
p1
p2+p3
break 
#--------------------------------s
saveRDS(rna,file = 'rna.rds')
saveRDS(atac,file='atac.rds')
rna<-readRDS('/home/em_b/Desktop/snRNA+ATAC_huPBMC/rna.rds')
atac<-readRDS('/home/em_b/Desktop/snRNA+ATAC_huPBMC/atac.rds')
