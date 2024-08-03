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
setwd('/home/em_b/work_stuff/human_develop_timecourse/GSE217511_RAW')
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720852_9C_barcodes.tsv.gz'
features_path <- 'GSM6720852_9C_features.tsv.gz'
matrix_path <- 'GSM6720852_9C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data1 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '17wks')
dim(data1)
data1<-subset(x = data1, downsample = 3000)
dim(data1)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720856_26C_barcodes.tsv.gz'
features_path <- 'GSM6720856_26C_features.tsv.gz'
matrix_path <- 'GSM6720856_26C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data2 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '19wks')
dim(data2)
data2<-subset(x = data2, downsample = 3000)
dim(data2)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720860_56C_barcodes.tsv.gz'
features_path <- 'GSM6720860_56C_features.tsv.gz'
matrix_path <- 'GSM6720860_56C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data3 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '20wks')
dim(data3)
data3<-subset(x = data3, downsample = 3000)
dim(data3)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720864_30C_barcodes.tsv.gz'
features_path <- 'GSM6720864_30C_features.tsv.gz'
matrix_path <- 'GSM6720864_30C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data4 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '22wks')
dim(data4)
data4<-subset(x = data4, downsample = 3000)
dim(data4)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720874_3C_barcodes.tsv.gz'
features_path <- 'GSM6720874_3C_features.tsv.gz'
matrix_path <- 'GSM6720874_3C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data5 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '26wks')
dim(data5)
data5<-subset(x = data5, downsample = 3000)
dim(data5)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720876_60C_barcodes.tsv.gz'
features_path <- 'GSM6720876_60C_features.tsv.gz'
matrix_path <- 'GSM6720876_60C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data6 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '32wks')
dim(data6)
data6<-subset(x = data6, downsample = 3000)
dim(data6)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720880_6C_barcodes.tsv.gz'
features_path <- 'GSM6720880_6C_features.tsv.gz'
matrix_path <- 'GSM6720880_6C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data7 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '38wks')
dim(data7)
data7<-subset(x = data7, downsample = 3000)
dim(data7)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720882_31C_barcodes.tsv.gz'
features_path <- 'GSM6720882_31C_features.tsv.gz'
matrix_path <- 'GSM6720882_31C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data8 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '28yrs')
dim(data8)
data8<-subset(x = data8, downsample = 3000)
dim(data8)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720884_2C_barcodes.tsv.gz'
features_path <- 'GSM6720884_2C_features.tsv.gz'
matrix_path <- 'GSM6720884_2C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data9 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '45yrs')
dim(data9)
data9<-subset(x = data9, downsample = 3000)
dim(data9)
#-------------------------------------------------------------------------------
barcodes_path <- 'GSM6720886_20C_barcodes.tsv.gz'
features_path <- 'GSM6720886_20C_features.tsv.gz'
matrix_path <- 'GSM6720886_20C_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
data10 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = '53yrs')
dim(data10)
data10<-subset(x = data10, downsample = 3000)
dim(data10)
#-------------------------------------------------------------------------------
dim(data1)
dim(data2)
dim(data3)
dim(data4)
dim(data5)
dim(data6)
dim(data7)
dim(data8)
dim(data9)
dim(data10)
data<-merge(data1,y=c(data3,data4,data5,data7,data8,data9,data10))
dim(data)
gc()
#---------------------------------------------------------------------------
table(data@active.ident)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt <20)
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
break 
saveRDS(data,file = 'save_file1.rda')
data<-readRDS('save_file1.rda')
break
data <- FindNeighbors(data, dims = 1:15)
data <- FindClusters(data, resolution = 0.05)
data <- RunUMAP(data, dims = 1:15)
data <-JoinLayers(data)
head(data@meta.data)
#-----For gene searches and labeling -------------------------------------------
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
         dims = c(1,2),
         group.by = 'cell_types')#,
         #cols = c('skyblue','yellow','red','blue','orange','green','brown','grey'))
DimPlot(data,
        label = TRUE,
        label.size = 6)
VlnPlot(data,features = c('OLIG1','OLIG2','SOX2','SOX9','MOG',
                                 'GFAP','S100B','SLC1A2',
                                 'RBFOX3','NEUROD2','CALM1',
                                 'AGBL4','PTPRD','MSR1','TCF4'))
Idents(data)
new_idents<-c('excitatory_neurons2','inhibitory_neurons','early_neurons','excitatory_neurons1','opc',
              'oligodendrocytes','bipolar_neurons','mature_neurons','microglia',
              'astrocytes')
wtf<-FindMarkers(data,
                 ident.1 = '0',
                 #ident.2 = 'microglia', 
                 #features = c(),
                 logfc.threshold=0,
                 only.pos = TRUE,
                 test.use = 'bimod',
                 min.pct = 0.5)
cat(row.names(wtf))
VlnPlot(data,features = c(row.names(wtf))[12:24])
#-------------------------------------------------------------------------------
names(new_idents) <- levels(data)
data <- RenameIdents(data, new_idents)
levels(data)
table(Idents(data))
Idents(data)
head(data@meta.data)
data@meta.data$cell_types<-Idents(data)
head(data@meta.data)
break 
saveRDS(data,file = 'save_file2.rda')
data<-readRDS('save_file2.rda')
break
#----switch active ident to time course
table(data@active.ident)
head(Idents(data))
Idents(data)<-data@meta.data$orig.ident
table(data@active.ident)
head(Idents(data))
#-------make cell data set
cds <- as.cell_data_set(data)
cds
head(colData(cds))
#----Edit partitions layer
cds@clusters@listData$UMAP$partitions
names(cds@clusters@listData$UMAP$partitions) <- cds@colData@rownames
cds
cds@clusters@listData
data@active.ident<-data@active.ident
cds@clusters@listData$UMAP$clusters<-data@active.ident
#----learn trajectory
cds<-learn_graph(cds,
                 use_partition = TRUE)
head(clusters(cds))
table(clusters(cds))
#----order cells from initial starting point
cds <- order_cells(cds, reduction_method = "UMAP",
                   root_cells = colnames(cds[,clusters(cds)=='17wks']))
#-----plots
plot_cells(cds,reduction_method = c('UMAP'),
           show_trajectory_graph = TRUE,
           trajectory_graph_color = 'red',
           group_label_size = 10,cell_size = 1,
           color_cells_by = c('pseudotime'))
DoHeatmap(
  data,
  features = c('OLIG1','OLIG2','SOX2','SOX9','MOG',
               'GFAP','S100B','SLC1A2',
               'RBFOX3','NEUROD2','CALM1',
               'AGBL4','PTPRD','MSR1','TCF4'),
  cells = NULL,
  group.by = "cell_types",
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
cds$monocle3_pseudotime <- pseudotime(cds)
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
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
colData(cds_3d)
plot_cells_3d(
  cds_3d,
  dims = c(1, 2, 3),
  reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
  color_cells_by = "cell_types",
  genes = NULL,
  show_trajectory_graph = FALSE,
  trajectory_graph_color = "red",
  trajectory_graph_segment_size = 10,
  norm_method = c("log"),
  color_palette = NULL,
  color_scale = "Viridis",
  cell_size =50 ,
  alpha = 1,
  min_expr = 0.1)
plot_cells_3d(
  cds_3d,
  dims = c(1, 2, 3),
  reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
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
  min_expr = 0.1
)
break
saveRDS(cds_3d,file = 'cds_3d.rda')
saveRDS(cds,file = 'cds.rda')
break 
#-----fit model statistics 
#-----This takes awhile
gene_fits<-fit_models(cds,
                      model_formula_str = '~pseudotime')
gc()
head(gene_fits)
#----fit_coefs show genes which vary as a function of time
fit_coefs<-coefficient_table(gene_fits)
fit_coefs <- fit_coefs %>% filter(term == 'pseudotime')
fit_coefs<-subset(fit_coefs,q_value< 0.05)
fit_coefs<-fit_coefs[,c(2,6,7,13)]
summary(fit_coefs)
fit_coefs<-fit_coefs[order(fit_coefs$estimate, decreasing=FALSE),]
fit_coefs<-fit_coefs[!grepl('LINC',fit_coefs$gene_id),]
fit_coefs<-fit_coefs[!grepl('AC0',fit_coefs$gene_id),]
fit_coefs<-fit_coefs[!grepl('AP0',fit_coefs$gene_id),]
fit_coefs<-fit_coefs[!grepl('AJ0',fit_coefs$gene_id),]
fit_coefs<-fit_coefs[!grepl('.1',fit_coefs$gene_id),]
fit_coefs<-fit_coefs[!grepl('.2',fit_coefs$gene_id),]
break 
write.csv(fit_coefs,file = 'pseudotime_gene_fits.csv')
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
setwd('/home/em_b/work_stuff/human_develop_timecourse/GSE217511_RAW')
data<-readRDS('save_file2.rda')
fit_coefs<-read.csv('pseudotime_gene_fits.csv')
cds<-readRDS('cds.rda')
break
#------plot genes by pseudotime
plot(fit_coefs$estimate ~ num(1:3977),
     bty='n',
     ylab='pseudotime',
     xlab='genes',
     main='Genes_ordered_by_pseudotime')
abline(lm(fit_coefs$estimate ~ num(1:3977)), col='blue')
text(fit_coefs$estimate,
     fit_coefs$gene_id,
     pos=4,
     offset = 2,
     cex=0.7,
     col = 'red')
#-------------------------------------------------------------------------------
DoHeatmap(
  data,
  features = c('FOXN4','MKI67','SQLE','LYPD6B','SHISA7','TMX3','SLITRK4',
               'ITGB4','ZNF98'),
  cells = NULL,
  group.by = "cell_types",
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
#----stop here
break 
#-----subset and graph gene by pseudotime
cds_subset<-cds[row.names(subset(rowData(cds),
                                 row.names(rowData(cds))%in%c('FOXN4','MKI67','SQLE','LYPD6B','SHISA7','TMX3','SLITRK4',
                                                              'ITGB4','ZNF98'))),]
cds_subset
colData(cds_subset)$pseudotime<-pseudotime(cds_subset,reduction_method = 'UMAP')
colData(cds_subset)
plot_genes_violin(cds_subset,
                  group_cells_by="monocle3_pseudotime",
                  ncol=3,
                  label_by_short_name = FALSE,
                  normalize = FALSE,
                  log_scale = FALSE)
#----plot histogram of pseudotime vs gene names
head(row.names(rowData(cds)))
rowData(cds)$gene_short_name<-row.names(rowData(cds))
rowData(cds)
plot_genes_in_pseudotime(cds, color_cells_by="monocle3_pseudotime")
plot_genes_in_pseudotime(cds_subset,
                         min_expr = NULL,
                         cell_size = 0.75,
                         nrow = NULL,
                         ncol = 1,
                         panel_order = NULL,
                         color_cells_by = "pseudotime",
                         trend_formula = "~ splines::ns(pseudotime, df=3)",
                         label_by_short_name = TRUE,
                         vertical_jitter = NULL,
                         horizontal_jitter = NULL)

break
write.csv(fit_coefs,file = 'pseudotime_gene_fits.csv')
break
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
VlnPlot(data_subset,features = c('SPARC','OLIG2','MOG','MSX3','MSX1'))
