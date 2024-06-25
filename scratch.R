library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
rna<-read.csv('/home/em_b/Desktop/rna_seq_current.csv')
h3k27me3<-read.csv('/home/em_b/Desktop/h3k27me3_matrix.csv')
thra<-read.csv('/home/em_b/Desktop/thra_matrix.csv')
colnames(thra)
thra<-thra %>%
  group_by(Gene) %>%
  summarize(THRA.0hr_1 = sum(THRA.0hr_1),
            THRA.0hr_2 = sum(THRA.0hr_2),
            THRA.0hr_3 = sum(THRA.0hr_3),
            THRA.T3_1 = sum(THRA.T3_1),
            THRA.T3_2 = sum(THRA.T3_2),
            THRA.T3_3 = sum(THRA.T3_3))


test<-merge(rna,h3k27me3,by='Gene')
test<-merge(test,thra,by='Gene')
row.names(test)<-test$Gene
test<-test[,-1]

matrix<-round(test)
data <- CreateSeuratObject(counts=matrix,project = 'mm39')
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
data <- RunPCA(data,npcs = 17, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:17, cells = 500, balanced = T)
ElbowPlot(data)
gc()
table(data@meta.data$orig.ident)
summary(matrix)


x<-FindMarkers(data, ident.1 = 'RNA.T3', ident.2 = 'RNA.0hr', 
               features = c(top1000),logfc.threshold=1.5,
               only.pos = FALSE,test.use = 'DESeq2',min.pct = 0.5)
x<-row.names(x)
x<-x[!grepl('Gm',x)]
x<-x[!grepl('Rik',x)]
VlnPlot(data,
        features = c(top1000)[1:20],
        pt.size = 1,
        cols = c('grey','red',
                 'grey','red',
                 'grey','red'))#,
        #layer='counts')
VlnPlot(data,
        features = c(top1000)[21:40],
        pt.size = 1,
        cols = c('grey','red',
                 'grey','red',
                 'grey','red'))
VlnPlot(data,
        features = c(top1000)[41:60],
        pt.size = 1,
        cols = c('grey','red',
                 'grey','red',
                 'grey','red'))
VlnPlot(data,
        features = c(top1000)[61:80],
        pt.size = 1,
        cols = c('grey','red',
                 'grey','red',
                 'grey','red'))
VlnPlot(data,
        features = c('Ubash3b','Ntrk3','Fgfr2','Anxa2','Adcy8','Ccdc162','Auts2','Cntnap5b'),
        pt.size = 1,
        cols = c('grey','red',
                 'grey','red',
                 'grey','red'))

