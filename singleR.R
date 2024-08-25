library(Signac)
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(htmlwidgets)
require('Seurat.utils')
library(plotly)
library(EnsDb.Mmusculus.v79)
library(celldex)
library(SingleR)
setwd('/home/deviancedev01/work_stuff/multiome_MB_tumorigenesis')
#-------------------------------------------------------------------------------
features_path <- '/home/deviancedev01/work_stuff/multiome_MB_tumorigenesis/p7_gnp/outs/filtered_feature_bc_matrix/features.tsv.gz'
barcodes_path <- '/home/deviancedev01/work_stuff/multiome_MB_tumorigenesis/p7_gnp/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
matrix_path <- '/home/deviancedev01/work_stuff/multiome_MB_tumorigenesis/p7_gnp/outs/filtered_feature_bc_matrix/matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
#head(colnames(matrix))
#head(row.names(matrix))
rna <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'rna')
summary(rna@active.ident)
#rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")
#VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
rna <- subset(rna, subset = nFeature_RNA > 200 & nFeature_RNA < 11000)# & percent.mt <5)
#VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
#FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
#head(VariableFeatures(rna))
#rna
#rna_df<-data.frame(rna[['RNA']]$counts)
#atac_df<-data.frame(atac[['ATAC']]$counts)
#top5000.rna <- head(VariableFeatures(rna), 5000)
#top5000.rna
plot1 <- VariableFeaturePlot(rna,raster = FALSE)
plot1
all.genes <- rownames(rna)
#all.genes
rna <- ScaleData(rna, features = all.genes)
#dim(rna)
rna <- RunPCA(rna, features = VariableFeatures(object = rna))
rna <- RunUMAP(rna, dims = 1:30)
DimHeatmap(rna, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(rna)
rna <- FindNeighbors(rna, dims = 1:30)
#--adjust resolution to atac umap clusters
rna <- FindClusters(rna,resolution = 0.10)
DimPlot(object = rna, label = TRUE,pt.size = 1)
#--auto annotate cell clusters
surveyReferences()
searchReferences('mouse')
ref <- fetchReference("mouse_rnaseq", "2024-02-26")
ref
#--- convert seurat object into single cell experiment object 
rna.sce <- as.SingleCellExperiment(DietSeurat(rna))
rna.sce
ref.main <- SingleR(test = rna.sce,
                       assay.type.test = 1,
                       ref = ref,
                       labels = ref$label.main)
ref.fine <- SingleR(test = rna.sce,
                    assay.type.test = 1,
                    ref = ref,
                    labels = ref$label.fine)
table(ref.main$pruned.labels)
table(ref.fine$pruned.labels)
head(rna@meta.data)
rna@meta.data$gen_celltype <- ref.main$pruned.labels
rna@meta.data$fine_celltype <- ref.fine$pruned.labels
head(rna@meta.data)
rna <- SetIdent(rna, value = "fine_celltype")
VlnPlot(rna,
           features=c('Rbfox3'))
#---have subset out all identities not NA 
identities<-c('NPCs','Neurons','OPCs','Fibroblasts activated','Astrocytes',
              'Neurons activated','Endothelial cells','aNSCs','Microglia','Ependymal',
              'Oligodendrocytes','T cells','Macrophages','Microglia activated','qNSCs',
              'Monocytes','Fibroblasts senescent','Dendritic cells','Fibroblasts',
              'Granulocytes')
rna<-subset(x = rna, idents = c(identities), invert = FALSE)
DoHeatmap(
  rna,
  features = c('Rbfox3','Ptch1','Gli2'),
  cells = NULL,
  #group.by = "cell_types",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 3.5,
  hjust = 0,
  vjust = 0,
  angle = 90,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)
table(rna@meta.data$fine_celltype)
#----downsample to see all groups
RNA<-subset(x = rna, downsample = 6)
DoHeatmap(
  RNA,
  features = c('Rbfox3','Ptch1','Gli2','Ptprc'),
  cells = NULL,
  #group.by = "cell_types",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 3.5,
  hjust = 0,
  vjust = 0,
  angle = 90,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)
