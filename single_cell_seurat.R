library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
setwd('/home/amp_prog/Downloads')
features_path <- 'GSE230295_L6features.tsv.gz'
barcodes_path <- 'GSE230295_L6barcodes.tsv.gz'
matrix_path <- 'GSE230295_L6matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
x <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'one_percent')
#---------------------------------------------------------------------------
features_path <- 'GSE230295_L7features.tsv.gz'
barcodes_path <- 'GSE230295_L7barcodes.tsv.gz'
matrix_path <- 'GSE230295_L7matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
y <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'ninety_percent')
#---------------------------------------------------------------------------
data<-merge(x,y,project='B_ALL')
rm(x,y,matrix,features_path,barcodes_path,matrix_path)
head(data@active.ident)
gc()
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt <10)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
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
#-------------CD45+ CD19+ TRAC+ CD14- MKI67+    -------------------------------------------
data$MKI67.groups <- 'MKI67.pos'
data$MKI67.groups[WhichCells(data, expression= MKI67 < 0.1)] <- 'MKI67.neg'
DimPlot(data, reduction = 'pca',split.by = 'MKI67.groups')
head(data@meta.data)
data <- subset(data, subset = MKI67.groups != "MKI67.neg")
gc()
VlnPlot(data, features = c('PTPRC','CD19','CD3D','CD14','MKI67'),pt.size=0.1)
#--------------------------------------------------------------
df <- as.data.frame(as.matrix(GetAssayData(data)))
Genes <- row.names(df)
df <- cbind(Genes,df)
View(df)
break  
