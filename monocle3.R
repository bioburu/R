sudo apt-get install libudunits2-dev
devtools::install_github('cole-trapnell-lab/monocle3')
install.packages('R.utils')
remotes::install_github('satijalab/seurat-wrappers')

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
cat('GSE234832. Patient-derived brain metastasis models')
setwd('/home/em_b/work_stuff/brain_metastasis/GSE234832_RAW')
barcodes_path <- 'GSM7475327_LUBMET7.barcodes.tsv.gz'
features_path <- 'GSM7475327_LUBMET7.features.tsv.gz'
matrix_path <- 'GSM7475327_LUBMET7.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'lung_cancer')
#---------------------------------------------------------------------------
setwd('/home/em_b/work_stuff/brain_metastasis/GSE234832_RAW')
barcodes_path <- 'GSM7475328_LUBMET1.barcodes.tsv.gz'
features_path <- 'GSM7475328_LUBMET1.features.tsv.gz'
matrix_path <- 'GSM7475328_LUBMET1.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data2 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'lung_cancer')
data<-merge(data,data2)
table(data@active.ident)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt <10)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top1000 <- head(VariableFeatures(data), 1000)
top1000
gc()
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
gc()
data <- ScaleData(data, features = all.genes)
gc()
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
new_idents<-c('OPC','astrocytes','intermediate','oligodendrocytes','react_astrocytes')
names(new_idents) <- levels(data_subset)
data_subset <- RenameIdents(data_subset, new_idents)
levels(data_subset)
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
VlnPlot(data_subset, features = c('OLIG2','SOX2','SOX9','SOX10','MOG','VIM',
                                  'GFAP','S100B','SLC1A2','SLC1A3','S100A10',
                                  'MAP2','NOTCH1','TUBB2B','GAS1','PDGFA','CLDN11',
                                  'CNP','MKI67'),
        cols = c())
#-------------------------------------------------------------------------------
cds <- as.cell_data_set(data_subset)
cds
head(colData(cds))
#-------------------------------------------------------------------------------
#----Edit partitions layer
cds@clusters@listData$UMAP$partitions
names(cds@clusters@listData$UMAP$partitions) <- cds@colData@rownames
cds
cds@clusters@listData
cds@clusters@listData$UMAP$clusters<-data_subset@active.ident
cds<-learn_graph(cds,
                 use_partition = FALSE)
head(clusters(cds))
cds <- order_cells(cds, reduction_method = "UMAP",
                   root_cells = colnames(cds[,clusters(cds)=='OPC']))
p1<-plot_cells(cds,
           reduction_method = c('UMAP'),
           show_trajectory_graph = TRUE,
           trajectory_graph_color = 'red',
           group_label_size = 10,
           cell_size = 1,
           color_cells_by = c('pseudotime'))
p2<-DimPlot(data_subset,
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
p1+p2
gene_fits<-fit_models(cds,model_formula_str = '~pseudotime')
head(gene_fits)
#----fit_coefs show genes which vary as a function of time
fit_coefs<-coefficient_table(gene_fits)
fit_coefs <- fit_coefs %>% filter(term == 'pseudotime')
fit_coefs<-subset(fit_coefs,q_value< 0.05)
fit_coefs<-fit_coefs[,c(2,6,7,13)]
#------pseudotime vs gene expression levels
cds_subset<-cds[row.names(subset(rowData(cds),
                                 row.names(rowData(cds))%in%c('OLIG2','SOX9','GFAP','VIM','S100B','MOG','SPP1','TF',
                                                              'NGFR','DPYSL5','TBC1D12','GSN','CAVIN1','ELOVL1','MARCKS'))),]
cds_subset
colData(cds_subset)$pseudotime<-pseudotime(cds_subset,reduction_method = 'UMAP')
colData(cds_subset)
plot_genes_violin(cds_subset,
                  group_cells_by="pseudotime",
                  ncol=3,
                  nrow = 5,
                  label_by_short_name = FALSE,
                  normalize = FALSE,
                  log_scale = FALSE)
#-------------------------------------------------------------------------------
cds_subset<-cds[row.names(subset(rowData(cds),
                                 row.names(rowData(cds))%in%c('HAPLN2','IER2','DNAJB1','BACE1','NRBP2',
                                                              'RAPGEF5','PMP2','MT3','CLU','GPR37L1'))),]
cds_subset
colData(cds_subset)$pseudotime<-pseudotime(cds_subset,reduction_method = 'UMAP')
colData(cds_subset)
plot_genes_violin(cds_subset,
                  group_cells_by="pseudotime",
                  ncol=2,
                  label_by_short_name = FALSE,
                  normalize = FALSE,
                  log_scale = FALSE)
cds_pr_test<-graph_test(cds,neighbor_graph = 'principal_graph',cores = 4)
cds_pr_test<-subset(cds_pr_test,q_value<0.05)
cds_pr_test<-cds_pr_test[order(cds_pr_test$morans_I, decreasing=TRUE),]
rowData(cds)$gene_short_name<-row.names(rowData(cds))
rowData(cds)
head(pseudotime(cds), 10)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(Idents(data_subset), monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot()
#----------------
data_subset <- AddMetaData(
  object = data_subset,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "astrocyte_trajectory"
)
FeaturePlot(data_subset, c("astrocyte_trajectory"), pt.size = 1) & scale_color_viridis_c()
FeaturePlot(data_subset,
            features = c(row.names(cds_pr_test)[1:20]))
VlnPlot(data_subset, features = c('MARCKS','CLU','S100B','MT3','TUBB2B','PLP1',
                                  'FOSB','GPR37L1','SLC1A3','HSPA1B',
                                  'EGR1','DNAJB1','SPP1','HAPLN2','TF',
                                  'IER2'))
# a helper function to identify the root principal points:
summary(cds@colData$seurat_clusters)
#---time bin is the earliest group assigned by user from colData metadata
get_earliest_principal_node <- function(cds, time_bin="0"){
  cell_ids <- which(colData(cds)[, "seurat_cluster"] == time_bin)
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
plot_cells(cds,
           color_cells_by = "seurat_clusters",
           label_cell_groups=TRUE,
           group_label_size = 10,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=2,
           cell_size = 1)
plot_cells_3d(
  cds_3d,
  dims = c(1, 2, 3),
  reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
  color_cells_by = "monocle3_pseudotime",
  genes = NULL,
  show_trajectory_graph = FALSE,
  trajectory_graph_color = "red",
  trajectory_graph_segment_size = 10,
  norm_method = c("log"),
  color_palette = NULL,
  color_scale = "Viridis",
  cell_size =200 ,
  alpha = 1,
  min_expr = 0.1
)
break 
