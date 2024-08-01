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
setwd('/home/em_b/work_stuff/prenatal_samples/GSE217511_RAW')
barcodes_path <- 'GSM6720852_9C_barcodes.tsv.gz'
features_path <- 'GSM6720852_9C_features.tsv.gz'
matrix_path <- 'GSM6720852_9C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'cortical_plate_17wks')
dim(data)
#data<-subset(x = data, downsample = 10000)
#---------------------------------------------------------------------------
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
data <- FindClusters(data, resolution = 0.1)
data <- RunUMAP(data, dims = 1:15)
data <-JoinLayers(data)
data_subset<-subset(data,ident=c(0,1,2,3,4))
new_idents<-c('glial','neurons','stem_like','intermediate','neural_precursor')
names(new_idents) <- levels(data_subset)
data_subset <- RenameIdents(data_subset, new_idents)
levels(data_subset)
table(Idents(data_subset))
DimPlot(data_subset,
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
VlnPlot(data_subset,features = c('OLIG2','SOX9','NES','VIM',
                          'GFAP','SOX2','S100B','SLC1A2','SLC1A3',
                          'MAP2','DIO2','RBFOX3','RORA'))
cds <- as.cell_data_set(data_subset)
cds
head(colData(cds))
#----Edit partitions layer
cds@clusters@listData$UMAP$partitions
names(cds@clusters@listData$UMAP$partitions) <- cds@colData@rownames
cds
cds@clusters@listData
cds@clusters@listData$UMAP$clusters<-data_subset@active.ident
cds<-learn_graph(cds,
                 use_partition = TRUE)
head(clusters(cds))
table(clusters(cds))
cds <- order_cells(cds, reduction_method = "UMAP",
                   root_cells = colnames(cds[,clusters(cds)=='stem_like']))
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
rowData(cds)$gene_short_name<-row.names(rowData(cds))
cds$monocle3_pseudotime <- pseudotime(cds)
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
break
saveRDS(cds_3d,file = 'cds_cortical_plate_17wks.rda')
saveRDS(data_subset,file = 'seurat_cortical_plate_17wks.rda')
plot_cells_3d(
  cds_3d,
  dims = c(1, 2, 3),
  reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
  color_cells_by = "monocle3_pseudotime",
  #genes = NULL,
  show_trajectory_graph = FALSE,
  trajectory_graph_color = "red",
  trajectory_graph_segment_size = 10,
  norm_method = c("log"),
  color_palette = NULL,
  color_scale = "Viridis",
  cell_size =50,
  alpha = 1,
  min_expr = 0.1)
colData(cds)
plot_cells_3d(
  cds_3d,
  dims = c(1, 2, 3),
  reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
  color_cells_by = "ident",
  genes = NULL,
  show_trajectory_graph = FALSE,
  trajectory_graph_color = "red",
  trajectory_graph_segment_size = 10,
  norm_method = c("log"),
  color_palette = NULL,
  color_scale = "Viridis",
  cell_size =50 ,
  alpha = 1,
  min_expr = 0.1
)
p3<-VlnPlot(data_subset,features = c('OLIG2','SOX9','RORA','MOG','VIM',
                          'GFAP','CD44','S100B','SLC1A2','SLC1A3',
                          'MAP2','DIO2','RBFOX3'))
p2+p3
wtf<-FindMarkers(data_subset,
                 ident.1 = 'astrocytes2',
                 #ident.2 = 'SHH.PBS', 
                 #features = c(),
                 logfc.threshold=1,
                 only.pos = TRUE,
                 test.use = 'bimod',
                 min.pct = 0.5)
cat(row.names(wtf))
VlnPlot(data_subset,features = c(row.names(wtf))[1:20])
break 
cds_3D<-readRDS('/home/em_b/Downloads/GSE217511_RAW/cortical_plate_17wks.rda')
break 
#-----fit model statistics 
gene_fits<-fit_models(cds,model_formula_str = '~pseudotime')
head(gene_fits)
#----fit_coefs show genes which vary as a function of time
fit_coefs<-coefficient_table(gene_fits)
fit_coefs <- fit_coefs %>% filter(term == 'pseudotime')
fit_coefs<-subset(fit_coefs,q_value< 0.05)
fit_coefs<-fit_coefs[,c(2,6,7,13)]
#-------------------------------------------------------------------------------
#-----principal graph test 
cds_pr_test<-graph_test(cds,neighbor_graph = 'principal_graph',cores = 4)
cds_pr_test<-subset(cds_pr_test,q_value<0.05)
cds_pr_test<-cds_pr_test[order(cds_pr_test$morans_I, decreasing=TRUE),]
#-------------------------------------------------------------------------------
#-----subset and graph gene by pseudotime
cds_subset<-cds[row.names(subset(rowData(cds),
                                 row.names(rowData(cds))%in%c('OLIG2','SOX9','SOX10','VIM',
                                                              'GFAP','CD44','S100B','SLC1A3',
                                                              'MAP2','DIO2','RBFOX3'))),]
cds_subset
colData(cds_subset)$pseudotime<-pseudotime(cds_subset,reduction_method = 'UMAP')
colData(cds_subset)
plot_genes_violin(cds_subset,
                  group_cells_by="pseudotime",
                  ncol=3,
                  label_by_short_name = FALSE,
                  normalize = FALSE,
                  log_scale = FALSE)
#-------------------------------------------------------------------------------
#-----plot pseudotime in ggplot 
data.pseudo <- as.data.frame(colData(cds))
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(Idents(data), monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot()
data_subset <- AddMetaData(
  object = data,
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
#summary(cds@colData$seurat_clusters)
#---time bin is the earliest group assigned by user from colData metadata
#get_earliest_principal_node <- function(cds, time_bin="5"){
#  cell_ids <- which(colData(cds)[, "seurat_cluster"] == time_bin)
#  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#  root_pr_nodes <-igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
#  root_pr_nodes
#}
#cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
#plot_cells(cds,
#           color_cells_by = "seurat_clusters",
#           label_cell_groups=TRUE,
#           group_label_size = 10,
#           label_leaves=FALSE,
#           label_branch_points=FALSE,
#           graph_label_size=2,
#           cell_size = 1)
break 
