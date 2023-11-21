library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(caTools)
library(car)
library(caret)
library(InformationValue)
library(pROC)
library(ROCR)
setwd('/home/deviancedev/Desktop/drive_nov2023/FCCC/alignments/finals')
matrix<-read.csv('gene.list_unordered.csv')
matrix<-matrix[order(matrix$Gene), ]
row.names(matrix)<-make.names(matrix$Gene,unique = TRUE)
matrix<-matrix[,-c(1:7)]
colnames(matrix)[2]<-'Norm_2'
#remove all none expressing genes
matrix<-filter(matrix,Norm_1>0,Norm_2>0,Norm_3>0,TAA_1>0,TAA_2>0,TAA_3>0)
#
data <- CreateSeuratObject(counts=matrix,project = 'TAA')
gc()
data@active.ident
table(data@active.ident)
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
dim(data)
data <- RunPCA(data,npcs = 5, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:5, cells = 500, balanced = T)
ElbowPlot(data)
gc()
table(data@meta.data$orig.ident)
#------------------------------------------
x<-FindMarkers(data, ident.1 = 'TAA', ident.2 = 'Norm', 
               features = c(top1000),logfc.threshold=1,min.pct1=1,
               max.pct2=0.0001,only.pos = TRUE,test.use = 'negbinom')
#x<-FindMarkers(data, ident.1 = 'TAA',ident.2 = 'Norm',
#               features = c(top1000),test.use='negbinom')
#x<-x[-c(215:997),]
x<-cbind(row.names(x),x)
colnames(x)[1]<-'Gene'
matrix2<-read.csv('gene.list_unordered.csv')
row.names(matrix2)<-make.names(matrix2$Gene,unique = TRUE)
matrix2<-matrix2[,-1]
matrix2<-cbind(row.names(matrix2),matrix2)
colnames(matrix2)[1]<-'Gene'
df<-merge(x,matrix2,by='Gene')
df<-df[order(df$avg_log2FC, decreasing=TRUE),]
break
out<-df
out$Gene
out$Gene<-sub('\\..*','',out$Gene)
out$Gene
#write.csv(out,file = 'most_sig_genes.csv')
break 
VlnPlot(data, features = c(row.names(x)[1:12]),cols = c('grey','red'))
FindMarkers(data, ident.1 = 'TAA', ident.2 = 'Norm', features = c(row.names(x)[1:12]),test.use = 'negbinom')

