library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(caTools)
#library(caret)
#library(InformationValue)
#library(pROC)
#library(ROCR)
library(readxl)
setwd('/home/em_b/work_stuff/FCCC/T3_RNAseq/')
matrix<-read.csv('biomart.hg38.csv')
row.names(matrix)<-make.names(matrix$Gene,unique = TRUE)
matrix<-matrix[,-c(1)]
data <- CreateSeuratObject(counts=matrix,project = 'hg38')
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
data <- RunPCA(data,npcs = 11, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:11, cells = 500, balanced = T)
ElbowPlot(data)
gc()
table(data@meta.data$orig.ident)
#-----correlations
FeatureScatter(data,
               feature1 = "NEUROD1",
               feature2 = 'EZH2',
               cols = c(),
               pt.size = 5,
               shuffle = TRUE,
               seed = 1,
               jitter = TRUE)
VlnPlot(data,
            features = c('NEUROD1','EZH2'),
        pt.size = 1,
        cols = c('grey','red',
                 'grey','red'),
        layer='counts')

DimPlot(data,
        cols = c('grey','red',
                 'grey','red'),
        pt.size = 7,
        dims = c(3,6))
#-----------PCA features
TopFeatures(data[['pca']],
            dim = 3,
            balanced = TRUE)
#------ridgeplot
RidgePlot(data,
          features = c("HR","KLF9"),
          cols = c('grey','red','grey','red'))
#------------------------------------------------
break
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

