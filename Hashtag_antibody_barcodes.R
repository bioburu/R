library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
setwd("/home/amp_prog/Downloads/GSE190964_RAW")
features_path <- "GSM5879885_features.tsv.gz"
barcodes_path <- "GSM5879885_barcodes.tsv.gz"
matrix_path <- "GSM5879885_matrix.mtx.gz"
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
dim(matrix)
str(matrix)
colnames(matrix)[1:10]
row.names(matrix)[1:10]
head(matrix[c('Hashtag1','Hashtag2','Hashtag6','Hashtag7','Hashtag8'),])
HTO<-matrix[c('Hashtag1','Hashtag2','Hashtag6','Hashtag7','Hashtag8'),]
HTO
#--------------------------------------------------------------------------------------------------------------
#------------Or-------------------------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
setwd('/home/amp_prog/Desktop/rstudio/CAR_T/GSE153697_+-CD19')
features_path <- "GSM4649255_mRNA_features.tsv.gz"
barcodes_path <- "GSM4649255_mRNA_barcodes.tsv.gz"
matrix_path <- "GSM4649255_mRNA_matrix.mtx.gz"
mRNA <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
features_path <- "GSM4649254_HTO_features.tsv.gz"
barcodes_path <- "GSM4649254_HTO_barcodes.tsv.gz"
matrix_path <- "GSM4649254_HTO_matrix.mtx.gz"
#--Set feature.column as needed Default=2------------------------------------------------
Hashtags <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path,feature.column = 1)
Hashtags
# Select cell barcodes detected by both RNA and HTO 
filtered <- intersect(colnames(mRNA), colnames(Hashtags))
# Subset RNA and HTO counts by joint cell barcodes
mRNA <- mRNA[, filtered]
Hashtags <- as.matrix(Hashtags[, filtered])
# Confirm that the HTO have the correct names
rownames(Hashtags)
#---------------------------------------------------------------------------------------------------------------------
data<-CreateSeuratObject(counts = matrix)
data<-NormalizeData(data)
data<-FindVariableFeatures(data)
data<-ScaleData(data,features = VariableFeatures(data))
data[['HTO']]<-CreateAssayObject(counts = HTO)
data<-NormalizeData(data, assay = 'HTO',normalization.method = 'CLR')
data<-HTODemux(data,assay = 'HTO',positive.quantile = 0.99)
table(data$HTO_classification.global)
VlnPlot(data,features = 'nCount_RNA',log = TRUE)
data<-subset(data,idents=c('Negative','Doublet'),invert=TRUE)
VlnPlot(data,features = 'nCount_RNA',log = TRUE)
head(data@meta.data)
head(data@active.ident)
head(data$nCount_RNA)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt <5)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
top10 <- head(VariableFeatures(data), 10)
top10
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
gc()
data
data <- RunPCA(data,npcs = 50, features = VariableFeatures(object = data))
#------------------------------------------------------------------
data<-RenameIdents(data,'Hashtag8'='untransduced','Hashtag2'='Transduced',
                   'Hashtag6'='Transduced','Hashtag1'='Transduced','Hashtag7'='untransduced')
head(data@meta.data)
dim(data)
cat('Get expression matrix from Seurat object')
#DF <- as.data.frame(as.matrix(GetAssayData(data)))
#Genes <- row.names(DF)
#DF <- cbind(Genes,DF)
#View(DF)
break
DimPlot(data, reduction = "pca", dims = c(1,2) ,pt.size = 3, label.box = TRUE, label.size = 2)
DimPlot(data, reduction = "pca", dims = c(3,4) ,pt.size = 3, label.box = TRUE, label.size = 2)

ProjectDim(data, reduction = "pca",dims.print = c(8))
