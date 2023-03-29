library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
setwd("/home/amp_prog/Downloads")
features_path <- "/home/amp_prog/Downloads/GSM6647634_HC_features.tsv.gz"
barcodes_path <- "/home/amp_prog/Downloads/GSM6647634_HC_barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Downloads/GSM6647634_HC_matrix.mtx.gz"
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
dim(matrix)
data <- CreateSeuratObject(counts = matrix, min.cells = 3, min.features = 200)
data$nFeature_RNA
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt <5)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top10 <- head(VariableFeatures(data), 10)
top10
plot1 <- VariableFeaturePlot(data)
#LabelPoints(plot = plot1, points= top10, repel = T)
plot1
all.genes <- rownames(data)
all.genes
data <- ScaleData(data, features = all.genes)
dim(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
DimPlot(data, reduction = "pca", pt.size = 7, label = FALSE)

#------------------------------------------------------------------------------------------------------------------------------------------
VlnPlot(data, features = c('CD3D','CD3G','CD3E','PTPRC','CD19','CD14'),pt.size=0.5)
data$TRAC.groups <- 'TRAC.pos'
data$TRAC.groups[WhichCells(data, expression= TRAC < 0.5)] <- 'TRAC.neg'
DimPlot(data, reduction = 'pca',split.by = 'TRAC.groups')
data@meta.data
data <- subset(data, subset = TRAC.groups != "TRAC.neg")
break
#--------------------------------monocle------------------------------------
#convert to monocle object
cds <- as.cell_data_set(data)
#get cell metadata
head(colData(cds))
#get gene data
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
head(fData(cds))
#get counts
head(counts(cds))
#retrieve cluster information from Seurat
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
#assign partitions
recreate.partitions <- as.factor(recreate.partitions)
head(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
#assign cluster information
list.cluster <- data@active.ident
head(list.cluster)
#assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- data@reductions$umap@cell.embeddings
#plot
cluster.before.traj <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
                                  group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj
#learn Trajectory
cds <- learn_graph(cds, use_partition = F)
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)
#order cells in pseudotime
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[,clusters(cds)==5]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)
#cells ordered by pseudotime
head(pseudotime(cds), 10)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
ggplot(data.pseudo, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters))+ geom_boxplot()
#reorder plot
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot()
#add pseudotime values into Seurat Object
data$pseudotime <- pseudotime(cds)
FeaturePlot(data, features = "pseudotime")
#Histograms of transcripts by cluster ID
RidgePlot(data, features = c("FLT3LG", "MILR1", "FLT3", "GZMB", "CD19"), sort = T, 
          idents = c("5", "6", "0", "7", "13"))
#---------------------------------------------------------------------------------------------
DF <- as.data.frame(as.matrix(GetAssayData(data)))
View(DF)
