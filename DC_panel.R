library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
setwd("/home/amp_prog/Desktop/rstudio/CRS_final/GSE176392_pDCs")
features_path <- "/home/amp_prog/Desktop/rstudio/CRS_final/GSE176392_pDCs/GSM5363809_genes.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CRS_final/GSE176392_pDCs/GSM5363809_barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CRS_final/GSE176392_pDCs/GSM5363809_matrix.mtx.gz"
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
data <- RunUMAP(data, dims = 1:8)
data
DimPlot(data, reduction = "umap")
break
VlnPlot(data, features = c("MILR1", "FLT3LG", "FLT3", "CD14","ITGAM","ITGAX","HLA-DRA",
                           "GZMA","GZMB","PRF1","LAMP1","FTH1"))#, slot = "counts", log = TRUE)
VlnPlot(data, features = c("IL1A","IL1B","IL3RA", "IL6", "IL6R", "IL6ST","CXCL8","CSF1",
                           "IL10","IL15","TGFB1","IL18"))
VlnPlot(data, features = c("IFNG","IFNGR1","TLR1","GAPDH","ACACA","IDH2","TNF",
                           "TNFAIP6","SOX4","JAK2","STAT3","CCL3"))#, slot = "counts", log = TRUE)
VlnPlot(data, features = c("CD3D", "CD4", "IL16","CD28","CD86","CD80","CD8A","IL2","CD19","CD22",
                           "NCAM1","FCGR3A"))
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
#?? plot
my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("CD19","FLT3LG",)))
