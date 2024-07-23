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
head(colData(cds))
recreate.partitions <- c(rep(1,length(cds@colData@rownames)))
head(recreate.partitions)
names(recreate.partitions) <- cds@colData@rownames
head(recreate.partitions)
recreate.partitions <- as.factor(recreate.partitions)
head(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
head(data_subset@active.ident)
list.cluster <- data_subset@active.ident
head(list.cluster)
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- data_subset@reductions$umap@cell.embeddings
cds<-learn_graph(cds,
                 use_partition = FALSE)
head(clusters(cds))
cds <- order_cells(cds, reduction_method = "UMAP",
                   root_cells = colnames(cds[,clusters(cds)=='OPC']))
plot_cells(cds,
           reduction_method = c('UMAP'),
           show_trajectory_graph = TRUE,
           trajectory_graph_color = 'red',
           group_label_size = 10,
           cell_size = 1,
           color_cells_by = c('pseudotime'))
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
my_genes <- row.names(subset(fData(cds), gene_short_name %in% c('GFAP','OLIG2','S100B','SOX9'))) 
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset,
                         cell_size = 5,
                         color_cells_by = "monocle3_pseudotime" )
FeaturePlot(data_subset, features = c('GFAP','OLIG2','S100B','SOX9'))
break 
VlnPlot(data_subset, features = c(row.names(cds_pr_test)[41:60]))
VlnPlot(data_subset, features = c(row.names(cds_pr_test)[61:80]))
VlnPlot(data_subset, features = c('MARCKS','CLU','S100B','MT3','TUBB2B','PLP1',
                                  'FOSB','GPR37L1','SLC1A3','HSPA1B',
                                  'EGR1','DNAJB1','SPP1','HAPLN2','TF',
                                  'IER2'))

