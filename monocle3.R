sudo apt-get install libudunits2-dev
devtools::install_github('cole-trapnell-lab/monocle3')
install.packages('R.utils')
remotes::install_github('satijalab/seurat-wrappers')
#------------------------------------------------------------------------------------------
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
setwd('/home/deviancedev01/work_stuff/human_develop_timecourse/GSE217511_RAW')
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720852_9C_barcodes.tsv.gz'
features_path <- 'GSM6720852_9C_features.tsv.gz'
matrix_path <- 'GSM6720852_9C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data1 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '17wks')
data1<-subset(x = data1, downsample = 3000)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720856_26C_barcodes.tsv.gz'
features_path <- 'GSM6720856_26C_features.tsv.gz'
matrix_path <- 'GSM6720856_26C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data2 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '19wks')
data2<-subset(x = data2, downsample = 3000)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720860_56C_barcodes.tsv.gz'
features_path <- 'GSM6720860_56C_features.tsv.gz'
matrix_path <- 'GSM6720860_56C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data3 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '20wks')
data3<-subset(x = data3, downsample = 3000)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720864_30C_barcodes.tsv.gz'
features_path <- 'GSM6720864_30C_features.tsv.gz'
matrix_path <- 'GSM6720864_30C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data4 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '22wks')
data4<-subset(x = data4, downsample = 3000)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720874_3C_barcodes.tsv.gz'
features_path <- 'GSM6720874_3C_features.tsv.gz'
matrix_path <- 'GSM6720874_3C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data5 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '26wks')
data5<-subset(x = data5, downsample = 3000)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720876_60C_barcodes.tsv.gz'
features_path <- 'GSM6720876_60C_features.tsv.gz'
matrix_path <- 'GSM6720876_60C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data6 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '32wks')
data6<-subset(x = data6, downsample = 3000)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720880_6C_barcodes.tsv.gz'
features_path <- 'GSM6720880_6C_features.tsv.gz'
matrix_path <- 'GSM6720880_6C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data7 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '38wks')
data7<-subset(x = data7, downsample = 3000)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720882_31C_barcodes.tsv.gz'
features_path <- 'GSM6720882_31C_features.tsv.gz'
matrix_path <- 'GSM6720882_31C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data8 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '28yrs')
data8<-subset(x = data8, downsample = 3000)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720884_2C_barcodes.tsv.gz'
features_path <- 'GSM6720884_2C_features.tsv.gz'
matrix_path <- 'GSM6720884_2C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data9 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '45yrs')
data9<-subset(x = data9, downsample = 3000)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720886_20C_barcodes.tsv.gz'
features_path <- 'GSM6720886_20C_features.tsv.gz'
matrix_path <- 'GSM6720886_20C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data10 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '53yrs')
data10<-subset(x = data10, downsample = 3000)
data<-merge(data1,y=c(data3,data4,data5,data7,data8,data9,data10))
table(data@meta.data$orig.ident)
#---------------------------------------------------------------------------
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt <20)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
data <- RunPCA(data, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
data<-JoinLayers(data)
rm(data1,data2,data3,data4,data5,
   data6,data7,data8,data9,data10,matrix)
gc()
data <- FindNeighbors(data, dims = 1:15)
data <- FindClusters(data, resolution = 0.05)
data <- RunUMAP(data, dims = 1:15)
#-----For gene searches and labeling -------------------------------------------
DimPlot(data,
         reduction = 'umap',
         label=TRUE,
         label.size = 3,
         label.box = FALSE,
         raster=FALSE,
         pt.size = 0.1,
         seed=1,
         cols.highlight = c('black'),
         dims = c(1,2))#+NoLegend()
new_idents<-c('excitatory_neurons2','inhibitory_neurons','early_neurons','excitatory_neurons1','opc',
              'oligodendrocytes','bipolar_neurons','mature_neurons','microglia','astrocytes')
names(new_idents) <- levels(data)
data <- RenameIdents(data, new_idents)
levels(data)
DoHeatmap(
  data,
  features = c('OLIG1','OLIG2','SOX2','SOX9','MOG',
               'GFAP','S100B','SLC1A2',
               'RBFOX3','NEUROD2','CALM1',
               'AGBL4','PTPRD','MSR1','TCF4'),
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
  angle = 50,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)


table(Idents(data))
data@meta.data$cell_types<-Idents(data)
head(data@meta.data)
#----switch active ident to time course
Idents(data)<-data@meta.data$orig.ident
table(data@active.ident)
cds <- as.cell_data_set(data)
head(colData(cds))
recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
#----Edit partitions layer
names(cds@clusters@listData$UMAP$partitions) <- cds@colData@rownames
#cds
#cds@clusters@listData
#data@active.ident<-data@active.ident
cds@clusters@listData$UMAP$clusters<-data@active.ident
#----learn trajectory
cds<-learn_graph(cds,
                 use_partition = TRUE)
#head(clusters(cds))
#table(clusters(cds))
#----order cells from initial starting point
cds <- order_cells(cds, reduction_method = "UMAP",
                   root_cells = colnames(cds[,clusters(cds)=='17wks']))
#-----plots
plot_cells(cds,reduction_method = c('UMAP'),
           show_trajectory_graph = TRUE,
           trajectory_graph_color = 'red',
           group_label_size = 10,cell_size = 1,
           color_cells_by = c('pseudotime'))
cds$monocle3_pseudotime <- pseudotime(cds)
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d@colData$monocle3_pseudotime
data@meta.data$monocle_pseudotime<-cds_3d@colData$monocle3_pseudotime
head(data@meta.data)
#---- 3d plot rna 
data3d <- data
data3d <- RunUMAP(data3d,dims = 1:15,n.components = 3L)
#head(Embeddings(data3d,reduction = "umap"))
head(data3d@meta.data)
plot3d1 <- FetchData(data3d, vars = c("umap_1", "umap_2", "umap_3", "cell_types","monocle_pseudotime"))
plot3d1$label <- paste(plot3d1$cell_types)
fig <- plot_ly(data = plot3d1, 
               x = ~umap_1, y = ~umap_2, z = ~umap_3, 
               color = ~monocle_pseudotime, 
               colors = c('viridis'),
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 2, width=2),
               text=~label,
               hoverinfo="text")%>%layout(title='Monocle3')
fig
plot_cells_3d(
  cds_3d,
  dims = c(1, 2, 3),
  reduction_method = c('PCA'),
  color_cells_by = "orig.ident",
  genes = NULL,
  show_trajectory_graph = FALSE,
  trajectory_graph_color = "red",
  trajectory_graph_segment_size = 10,
  norm_method = c("log"),
  color_palette = NULL,
  color_scale = "Viridis",
  cell_size =50 ,
  alpha = 1,
  min_expr = 0.1)%>%layout(title='PCA')
