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
data <- FindNeighbors(data)
data <- FindClusters(data, resolution = c(0.3,0.5,0.7,1))
head(data@meta.data)
p1 <- DimPlot(data, group.by = "RNA_snn_res.0.3", label = T)
p2 <- DimPlot(data, group.by = "RNA_snn_res.0.5", label = T)
p3 <- DimPlot(data, group.by = "RNA_snn_res.0.7", label = T)
p4 <- DimPlot(data, group.by = "RNA_snn_res.1", label = T)
p1 + p2 + p3 + p4
Idents(data) <- "RNA_snn_res.0.7"
head(data@meta.data)
data <- RunUMAP(data, dims = 1:9)
data
DimPlot(data, reduction = "umap")

VlnPlot(data, features = c("MILR1", "FLT3LG", "FLT3", "CD14","ITGAM","ITGAX","IL6",
                           "HLA-DRA","IL8","CLEC2D","CYFIP2","IFNG"))#, slot = "counts", log = TRUE)
VlnPlot(data, features = c("MILR1","IFNGR1","IL4R","IFNGR2","PTPN22","TNF","SOX4",
                           "JAK2","STAT3","NLRP3","IL1B","IL18"))#, slot = "counts", log = TRUE)
VlnPlot(data, features = c("CD3D","CD19","NCAM1","CD8A"))#, slot = "counts", log = TRUE)
FeaturePlot(data, features = c("MILR1", "FLT3LG", "FLT3", "CD14","ITGAM","ITGAX","IL6",
                               "HLA-DRA","CCL2","CLEC2D","CYFIP2","CCR5","CD209","IL4R",
                               "IFNGR2", "PTPN22"))

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
