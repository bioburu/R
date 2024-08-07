library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
library(SeuratWrappers)
library(monocle3)
library(plotly)
cat('GSE234832. Patient-derived brain metastasis models')
setwd('/home/em_b/work_stuff/brain_metastasis/GSE234832_RAW')
barcodes_path <- 'GSM7475327_LUBMET7.barcodes.tsv.gz'
features_path <- 'GSM7475327_LUBMET7.features.tsv.gz'
matrix_path <- 'GSM7475327_LUBMET7.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data1 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'lung_cancer')
#---------------------------------------------------------------------------
setwd('/home/em_b/work_stuff/brain_metastasis/GSE234832_RAW')
barcodes_path <- 'GSM7475328_LUBMET1.barcodes.tsv.gz'
features_path <- 'GSM7475328_LUBMET1.features.tsv.gz'
matrix_path <- 'GSM7475328_LUBMET1.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data2 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'lung_cancer')
data<-merge(data1,data2)
table(data2@active.ident)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt <10)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top1000 <- head(VariableFeatures(data), 1000)
head(top1000)
gc()
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
head(all.genes)
data <- ScaleData(data, features = all.genes)
dim(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
data <- FindNeighbors(data, dims = 1:15)
data <- FindClusters(data, resolution = 0.6)
data <- RunUMAP(data, dims = 1:15)
data <-JoinLayers(data)
DimPlot(data,
        reduction = 'umap',
        label=TRUE,
        label.size = 6,
        repel = TRUE,
        label.box = FALSE,
        raster=FALSE,
        pt.size = 0.5,
        seed=1,
        cols.highlight = c('grey'),
        dims = c(1,2))
#-------------------------------------------------------------------------------
data_subset<-subset(data,ident=c(4,6,9))
data_subset <- FindNeighbors(data_subset, dims = 1:15)
data_subset <- FindClusters(data_subset, resolution = 0.7)
data_subset <- RunUMAP(data_subset, dims = 1:15)
new_idents<-c('OPC','intermediates','reactive_OPC','oligodendrocytes','astrocytes')
names(new_idents) <- levels(data_subset)
data_subset <- RenameIdents(data_subset, new_idents)
levels(data_subset)
data_subset@meta.data$cell_types<-Idents(data_subset)
table(data_subset@meta.data$cell_types)
DimPlot(data_subset,
        reduction = 'umap',
        label=TRUE,
        label.size = 6,
        repel = TRUE,
        label.box = FALSE,
        raster=FALSE,
        pt.size = 1,
        seed=1,
        cols.highlight = c('grey'),
        dims = c(1,2))
DoHeatmap(
  data_subset,
  features = c('GFAP','S100B','SLC1A3','IL6ST','OLIG2','ZNF488','MOG','VIM','CD44','LAMP2'),
  cells = NULL,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  vjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)
#-------------------------------
break 
require('Seurat.utils')
data3d <- data_subset
data3d <- RunUMAP(data3d,dims = 1:15,n.components = 3L)
head(Embeddings(data3d,reduction = "umap"))
plot3d1 <- FetchData(data3d, vars = c("umap_1", "umap_2", "umap_3", "cell_types"))
plot3d1$label <- paste(rownames(plot3d1))
fig <- plot_ly(data = plot3d1, 
               x = ~umap_1, y = ~umap_2, z = ~umap_3, 
               color = ~cell_types, 
               colors = c("lightseagreen",
                          "gray50",
                          "darkgreen",
                          "red4",
                          "red",
                          "turquoise4",
                          "black",
                          "yellow4",
                          "royalblue1",
                          "lightcyan3",
                          "peachpuff3",
                          "khaki3",
                          "gray20",
                          "orange2",
                          "royalblue4",
                          "yellow3",
                          "gray80",
                          "darkorchid1",
                          "lawngreen",
                          "plum2",
                          "darkmagenta"),
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 5, width=2),
               text=~label,
               hoverinfo="text")
fig
goi <- "GFAP"
plotting.data <- FetchData(data3d, vars = c("umap_1", "umap_2", "umap_3", "Expression"=goi), slot = 'data')
Cutoff <- quantile(plotting.data[,goi], probs = .95)
plotting.data$"ExprCutoff" <- ifelse(test = plotting.data[,goi] <Cutoff, yes = plotting.data[,goi], no = Cutoff)
plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data[,goi], sep="")
plot_ly(data = plotting.data,
        name = goi,
        x = ~umap_1, y = ~umap_2, z = ~umap_3, 
        color = ~ExprCutoff, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgrey', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 5), 
        text=~label,
        hoverinfo="text"
) %>%layout(title=goi)
#-----If subsetting use
#----------Isolate CD45+  -------------------------------------
data$CD45.groups <- 'CD45.pos'
data$CD45.groups[WhichCells(data, expression= PTPRC < 0.1)] <- 'CD45.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD45.groups')
head(data@meta.data)
data <- subset(data, subset = CD45.groups != "CD45.neg")
gc()
#-------------CD45+ CD19+  ---------------------------------------------------------
data$CD19.groups <- 'CD19.pos'
data$CD19.groups[WhichCells(data, expression= CD19 < 0.1)] <- 'CD19.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD19.groups')
head(data@meta.data)
data <- subset(data, subset = CD19.groups != "CD19.neg")
gc()
#-------------CD45+ CD19+ TRAC- -------------------------------------------
data$TRAC.groups <- 'TRAC.pos'
data$TRAC.groups[WhichCells(data, expression= TRAC < 0.1)] <- 'TRAC.neg'
DimPlot(data, reduction = 'pca',split.by = 'TRAC.groups')
head(data@meta.data)
data <- subset(data, subset = TRAC.groups != "TRAC.pos")
gc()
#-------------CD45+ CD19+ TRAC- CD14-     -------------------------------------------
data$CD14.groups <- 'CD14.pos'
data$CD14.groups[WhichCells(data, expression= CD14 < 0.1)] <- 'CD14.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD14.groups')
head(data@meta.data)
data <- subset(data, subset = CD14.groups != "CD14.pos")
gc()
#-------------CD45+ CD19+ TRAC- CD14- MKI67+    -------------------------------------------
data$MKI67.groups <- 'MKI67.pos'
data$MKI67.groups[WhichCells(data, expression= MKI67 < 0.1)] <- 'MKI67.neg'
DimPlot(data, reduction = 'pca',split.by = 'MKI67.groups')
head(data@meta.data)
data <- subset(data, subset = MKI67.groups != "MKI67.neg")
gc()
VlnPlot(data, features = c('PTPRC','CD19','TRAC','CD14','MKI67'),pt.size=0.1)

