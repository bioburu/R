library(GEOquery)
library(biomaRt)
library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
setwd('/home/amp_prog/rstudio')
filePaths = getGEOSuppFiles("GSE64346")
filePaths
matrix <- read.delim('/home/amp_prog/rstudio/GSE64346/GSE64346_Non-normalized_data.txt')
matrix <- matrix[,-c(6:11)]
#View(matrix)
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl <- useMart("ensembl")
listDatasets(ensembl)
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
filters <- listFilters(ensembl)
filters[1:10,]
attributes <- listAttributes(ensembl)
attributes[1:50,]
attributes
geneid <- matrix$Probe
head(geneid)
dim(geneid)
genes <- getBM(attributes = c('affy_primeview','ensembl_gene_id','external_gene_name','description'),
               filters = 'affy_primeview',
               values = geneid,
               mart = ensembl)
head(genes)
View(genes)
df <- matrix[matrix$Probe %in% unique(genes$affy_primeview),]
dim(df)
dim(matrix)
View(df)
genes$affy_primeview
View(genes)
colnames(genes)[1] <- 'Probe.Set.ID'
table <- merge(genes,df,by="Probe.Set.ID")
View(table)
table <- table[!(is.na(table$external_gene_name) | table$external_gene_name==""),]
dim(table)
row.names(table) <- make.names(table$external_gene_name, unique = TRUE)
table <-table[,-c(1:4)]
data <- CreateSeuratObject(counts = table)
dim(data)
data$nFeature_RNA
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
#data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt <5)
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
data <- RunPCA(data,npcs=3, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:3, cells = 500, balanced = T)
ElbowPlot(data)
data <- FindNeighbors(data)
data <- FindClusters(data, resolution = c(0.3,0.5,0.7,1))
head(data@meta.data)
p1 <- DimPlot(data, group.by = "RNA_snn_res.0.3", label = T)
p2 <- DimPlot(data, group.by = "RNA_snn_res.0.5", label = T)
p3 <- DimPlot(data, group.by = "RNA_snn_res.0.7", label = T)
p4 <- DimPlot(data, group.by = "RNA_snn_res.1", label = T)
p1 + p2 + p3 + p4
Idents(data) <- "RNA_snn_res.0.7"
head(data@meta.data)
data <- RunUMAP(data, dims = 1:9)
data
DimPlot(data, reduction = "pca")
