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
library(biomaRt)
library(ReactomePA)
library(clusterProfiler)
setwd('/home/deviancedev01/Downloads/GSE224679_RAW')
#-------------------------------------------------------------------------------
features_path <- 'GSM7029393_111_11_features.tsv.gz'
barcodes_path <- 'GSM7029393_111_11_barcodes.tsv.gz'
matrix_path <- 'GSM7029393_111_11_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
lesions <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'lesions')
summary(lesions@active.ident)
#-----------------------------
features_path <- 'GSM7029395_111_13_features.tsv.gz'
barcodes_path <- 'GSM7029395_111_13_barcodes.tsv.gz'
matrix_path <- 'GSM7029395_111_13_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
normal <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'normal')
summary(normal@active.ident)
lesions<-subset(lesions,
                downsample=2760)
rna<-merge(normal,
           lesions)
table(Idents(rna))
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^mt-")
VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
rna <- subset(rna, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt <30)
VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
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
rna <- FindNeighbors(rna, dims = 1:8)
#--adjust resolution to atac umap clusters
rna <- FindClusters(rna,resolution = 0.3)
DimPlot(object = rna, label = TRUE,pt.size = 1)
#--auto annotate cell clusters
surveyReferences()
searchReferences('mouse')
ref <- fetchReference("mouse_rnaseq", "2024-02-26")
ref
rna<-JoinLayers(rna)
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
        features=c('Cd14','Cd3d'))
DoHeatmap(
  rna,
  features = c('Cd14','Cd3d'),
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
DimPlot(rna,
        reduction = "umap",
        group.by = c("orig.ident", "fine_celltype"),
        label = TRUE)+NoLegend()
#----perform integration
rna
rna[["RNA"]] <- split(rna[["RNA"]],
                      f = rna$orig.ident)
rna <- IntegrateLayers(rna,
                       method = CCAIntegration,
                       orig.reduction = "pca",
                       new.reduction = "integrated.cca",
                       verbose = FALSE)
# re-join layers after integration
rna[["RNA"]] <- JoinLayers(rna[["RNA"]])
rna <- FindNeighbors(rna, reduction = "integrated.cca", dims = 1:30)
rna <- FindClusters(rna, resolution = 0.3)
rna <- RunUMAP(rna, dims = 1:30, reduction = "integrated.cca")
DimPlot(rna,
        reduction = "umap",
        group.by = c("orig.ident", "fine_celltype"),
        label = TRUE)+NoLegend()
#---comparisons between cell types
rna<-SetIdent(rna,
              value='fine_celltype')
table(Idents(rna))
rna_subset<-subset(rna,
                    idents = 'T cells')
rna_subset<-SetIdent(rna_subset,
                      value='orig.ident')
table(Idents(rna_subset))
VlnPlot(rna_subset,
        features = c('Cd3d'))
markers <- FindMarkers(rna_subset,
                       ident.1 = 'lesions',
                       ident.2 = 'normal',
                       test.use='bimod',
                       #only.pos=TRUE,
                       min.pct=0.7,
                       logfc.threshold = 1)
markers<-subset(markers,
                p_val<0.05)
markers
geneid <- row.names(markers)
tcell_subset <- RenameIdents(object = rna_subset,
                           `lesions` = "Tcell.lesions",
                           `normal`="Tcell.normal")
Idents(rna_subset)
DoHeatmap(
  tcell_subset,
  features = c(geneid),
  cells = NULL,
  #group.by = "orig.ident",
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
  angle = 50,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)
break 
#---if subset out all identities not NA 
identities<-c('NPCs','Neurons','OPCs','Fibroblasts activated','Astrocytes',
              'Neurons activated','Endothelial cells','aNSCs','Microglia','Ependymal',
              'Oligodendrocytes','T cells','Macrophages','Microglia activated','qNSCs',
              'Monocytes','Fibroblasts senescent','Dendritic cells','Fibroblasts',
              'Granulocytes')
rna<-subset(x = rna, idents = c(identities), invert = FALSE)
#----downsample to see all groups
RNA<-subset(x = rna, downsample = 6)
DoHeatmap(
  RNA,
  features = c('Cd14','Cd3d'),
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
