library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
setwd("/home/GSE190964_RAW")
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
#------------Or-------------------------------------------------------------------------------------------------
setwd('GSE153697_+-CD19')
meta<- read.delim('GSE153697_filtered_metadata_matrix.tsv')
meta$hash.ID<-as.factor(meta$hash.ID)
summary(meta)
meta<- read.delim('GSE153697_filtered_metadata_matrix.tsv')
meta$hash.ID[meta$hash.ID == "T1-CD19neg"] <- "pretreatmentCD19neg"
meta$hash.ID[meta$hash.ID == "T1-CD19pos"] <- "pretreatmentCD19pos"
meta$hash.ID[meta$hash.ID == "T2-CD19neg"] <- "treatedCD19neg"
meta$hash.ID[meta$hash.ID == "T2-CD19pos"] <- "treatedCD19pos"
summary(meta)
matrix<-read.delim('GSE153697_filtered_raw_expression_matrix.tsv')
head(matrix)[1:10]
dim(matrix)
dim(meta)
colnames(matrix) <- meta$hash.ID[match(colnames(matrix), meta$ID)]
rownames(matrix)<-matrix[,1]
matrix<-matrix[,-1]
data<-CreateSeuratObject(counts=matrix)
data
head(data@meta.data)
data<-AddMetaData(data,meta$hash.ID,col.name = 'hash.ID')
Idents(data)<-data$hash.ID
data@active.ident
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
DimPlot(data, reduction = "pca", pt.size = 1, label = FALSE)
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
DimPlot(data, reduction = "pca", dims = c(1,2) ,pt.size = 3, label.box = TRUE, label.size = 2)
DimPlot(data, reduction = "pca", dims = c(3,4) ,pt.size = 3, label.box = TRUE, label.size = 2)

ProjectDim(data, reduction = "pca",dims.print = c(8))
