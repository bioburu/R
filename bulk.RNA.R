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
library(readxl)
setwd('/home/deviancedev/Desktop/drive_nov2023/FCCC/GSE154958_GBM_dura_NSC')
matrix<-read_excel('GSE154958_cell_culture_star_counts.xlsx')
matrix<-data.frame(matrix)
str(matrix)
matrix<-data.frame(matrix %>%
  group_by(Gene) %>%
  summarize(T3_1 = sum(T3_1),
            T3_2 = sum(T3_2),
            T3_3 = sum(T3_3),
            PBS_1 = sum(PBS_1),
            PBS_2 = sum(PBS_2),
            PBS_3 = sum(PBS_3)))
str(matrix)
row.names(matrix)<-make.names(matrix$Gene.Symbol,unique = TRUE)
matrix<-matrix[,-c(1)]
data <- CreateSeuratObject(counts=matrix,project = 'gbm')
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
data <- RunPCA(data, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
gc()
table(data@meta.data$orig.ident)
#------------------------------------------
x<-FindMarkers(data, ident.1 = 'GBM', ident.2 = 'NSC', 
               features = c(top1000),logfc.threshold=2,
               only.pos = FALSE,test.use = 'DESeq2',min.pct = 0.5)
VlnPlot(data,features = c(rownames(x)[1:12]),cols = c('grey','red','grey','red'))
VlnPlot(data,features = c(rownames(x)[37:48]),cols = c('grey','red','grey','red'))

VlnPlot(data,features = c('NEUROD1','SCG2','PMP2','GFAP','MLC1','TIMP4',
                          'NOTCH1','PROM1','THRA','THRB','DIO2','DIO3','NCAM1','FOXO1','FOXG1','MKI67'),cols = c('grey','red','grey','red'))
FindMarkers(data, ident.1 = 'GBM', ident.2 = 'NSC', features = c('NEUROD1','SCG2','PMP2','GFAP','MLC1','TIMP4',
                                                                 'NOTCH1','PROM1','THRA','THRB','DIO2','DIO3','NCAM1','FOXO1','FOXG1','MKI67'),test.use = 'negbinom')
x<-subset(x,p_val_adj< 0.01)
downreg<-subset(x,avg_log2FC< -1.5)
upreg<-subset(x,avg_log2FC> 1.5)
df<-rbind(upreg,downreg)
df<-cbind(row.names(df),df)
break
write.csv(x,file = 'GSE154958_gbm_dura.csv')

