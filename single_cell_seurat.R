library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
setwd('/home/amp_prog/Desktop/hope_model/GSE208653_hpv_cancer')
features_path <- 'GSM6360680_N_HPV_NEG_1.features.tsv.gz'
barcodes_path <- 'GSM6360680_N_HPV_NEG_1.barcodes.tsv.gz'
matrix_path <- 'GSM6360680_N_HPV_NEG_1.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
no_hpv <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'no_hpv')
#---------------------------------------------------------------------------
features_path <- 'GSM6360682_N_1.features.tsv.gz'
barcodes_path <- 'GSM6360682_N_1.barcodes.tsv.gz'
matrix_path <- 'GSM6360682_N_1.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
hpv <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'hpv')
#---------------------------------------------------------------------------
features_path <- 'GSM6360684_HSIL_1.features.tsv.gz'
barcodes_path <- 'GSM6360684_HSIL_1.barcodes.tsv.gz'
matrix_path <- 'GSM6360684_HSIL_1.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
lesions <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'lesions')
#---------------------------------------------------------------------------
features_path <- 'GSM6360686_SCC_4.features.tsv.gz'
barcodes_path <- 'GSM6360686_SCC_4.barcodes.tsv.gz'
matrix_path <- 'GSM6360686_SCC_4.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
cervical_cancer <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'cervical_cancer')
data<-merge(no_hpv,y=c(hpv,lesions,cervical_cancer),project='hpv_cervical_cancer')
head(data@active.ident)
rm(no_hpv,hpv,lesions,cervical_cancer,matrix)
gc()
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt <10)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top10 <- head(VariableFeatures(data), 10)
top10
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
gc()
#----------Isolate CD45+  -------------------------------------
data$CD45.groups <- 'CD45.pos'
data$CD45.groups[WhichCells(data, expression= PTPRC < 0.1)] <- 'CD45.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD45.groups')
head(data@meta.data)
exp <- subset(data, subset = CD45.groups != "CD45.neg")
rm(data,barcodes_path,features_path,matrix_path)
gc()
#-------------CD45+ CD19-  ---------------------------------------------------------
exp$CD19.groups <- 'CD19.pos'
exp$CD19.groups[WhichCells(exp, expression= CD19 < 0.1)] <- 'CD19.neg'
DimPlot(exp, reduction = 'pca',split.by = 'CD19.groups')
head(exp@meta.data)
exp <- subset(exp, subset = CD19.groups != "CD19.pos")
gc()
#-------------CD45+ CD19+ TRAC+ -------------------------------------------
exp$TRAC.groups <- 'TRAC.pos'
exp$TRAC.groups[WhichCells(exp, expression= TRAC < 0.1)] <- 'TRAC.neg'
DimPlot(exp, reduction = 'pca',split.by = 'TRAC.groups')
head(exp@meta.data)
exp <- subset(exp, subset = TRAC.groups != "TRAC.neg")
gc()
#-------------CD45+ CD19+ TRAC+ CD14-     -------------------------------------------
exp$CD14.groups <- 'CD14.pos'
exp$CD14.groups[WhichCells(exp, expression= CD14 < 0.1)] <- 'CD14.neg'
DimPlot(exp, reduction = 'pca',split.by = 'CD14.groups')
head(exp@meta.data)
exp <- subset(exp, subset = CD14.groups != "CD14.pos")
gc()
#-------------CD45+ CD19+ TRAC+ CD14- MKI67+    -------------------------------------------
exp$MKI67.groups <- 'MKI67.pos'
exp$MKI67.groups[WhichCells(exp, expression= MKI67 < 0.1)] <- 'MKI67.neg'
DimPlot(exp, reduction = 'pca',split.by = 'MKI67.groups')
head(exp@meta.data)
exp <- subset(exp, subset = MKI67.groups != "MKI67.neg")
gc()
VlnPlot(exp, features = c('PTPRC','CD19','TRAC','CD14','MKI67'),pt.size=0.1)
break 
