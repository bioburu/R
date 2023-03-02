library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
setwd("/home/amp_prog/Desktop/rstudio/CAR_M")
matrix <- read.delim('GSE120086_Raw_read_counts_challengeData.txt')
dim(matrix)
row.names(matrix) <- make.names(matrix$refGene, unique = TRUE)
matrix <- matrix[,-1]
data <- CreateSeuratObject(counts = matrix)
dim(data)
data$orig.ident
data$nFeature_RNA
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top10 <- head(VariableFeatures(data), 10)
top10
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
data <- ScaleData(data, features = all.genes)
dim(data)
data <- RunPCA(data,npcs=11, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:11, cells = 500, balanced = T)
ElbowPlot(data)
data
plot <-DimPlot(data, reduction = "pca", pt.size = 7, label = FALSE)
plot
data@meta.data
plot$data
break
#---------------------------------------------------------------------------
#------------------MACRO PANEL----------------------------------------
VlnPlot(data, features = c("CD14","FCGR3A","ITGAM","ITGAX","CD38","MILR1",
                           "CREB5","MXD3","EGR1","SPI1","TLR1","TLR2","TLR3",
                           "TLR4","TLR5","TLR6"), cols  = c('red','green','blue'))
RidgePlot(data, features = c("CD14","FCGR3A","ITGAM","ITGAX","CD38","MILR1",
                           "CREB5","MXD3","EGR1","SPI1","TLR1","TLR2","GAPDH",
                           "IDH2","TLR5","TLR6"), cols  = c('red','green','blue'))
VlnPlot(data, features = c("TLR7","TLR8","TLR8-AS1","TLR9","CCR5",
                           "HLA-A","HLA-B","HLA-DRA","HLA-DRB1","HLA-C",
                           "HLA-DMA","HLA-DMB","HLA-DPA1","HLA-DPB1","HLA-DQA1",
                           "HLA-DQB1"), cols  = c('red','green'))
RidgePlot(data, features = c("TLR7","TLR8","TLR8-AS1","TLR9","CCR5",
                           "HLA-A","HLA-B","HLA-DRA","HLA-DRB1","HLA-C",
                           "HLA-DMA","HLA-DMB","HLA-DPA1","HLA-DPB1","HLA-DQA1",
                           "HLA-DQB1"), cols  = c('red','green','blue'))
VlnPlot(data, features = c("HLA-DRB5","HLA-E","HLA-F","HLA-G","CD274","CD80",
                           "CD86","SELL","SELPLG","SPN","ITGA4",
                           "CD163","CD163L1","CD74","CCR5","CCRL2"), cols  = c('red','green'))
RidgePlot(data, features = c("HLA-DRB5","HLA-E","HLA-F","HLA-G","CD274","CD80",
                           "CD86","SELL","SELPLG","SPN","ITGA4",
                           "CD163","CD163L1","CD74","CCR5","CCRL2"), cols  = c('red','green'))
VlnPlot(data, features = c("IL1B","IL1RN","IL5RA","IL7R","CXCL8","IL10",
                           "IL17RA","IL18BP","TGFBI","TNFRSF4","CSF1",
                           "CXCL12"), cols  = c('red','green'))
RidgePlot(data, features = c("IL1B","IL1RN","IL5RA","IL7R","CXCL8","IL10",
                           "IL17RA","IL18BP","TGFBI","TNFRSF4","CSF1",
                           "CXCL12"), cols  = c('red','green'))
VlnPlot(data, features = c('CFD','C5'), cols  = c('red','green'))
VlnPlot(data, features = c('SERPING1','VSIG4','VTN'), cols  = c('red','green'))
#-----Final candidate genes ----------------------------------------------------
VlnPlot(data, features = c('CFD','C5','SERPING1','VSIG4','VTN'), cols  = c('red','green','blue'))


RidgePlot(data, features = c('CD14','ITGAM','ITGAX','MXD3','TLR2','TLR8.AS1',
                           'CD274','SPN','CD163L1','ITGA4','IL1B','CXCL8','TGFBI',
                           'CFD'), cols  = c('red','green','blue'))
RidgePlot(data, features = c('ITGAM', 'ITGAX', 'MXD3', 'SPI1', 'TLR7', 'HLA-DRA', 
                           'HLA-DRB1','HLA-DRB5', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 
                           'HLA-DQB1'), cols  = c('red','green'))
VlnPlot(data, features = c('HLA-C', 'ITGA4', 'IL1B', 'IL1RN', 'IL10', 
                           'IL18BP', 'TGFBI', 'TNFRSF4', 'CXCL12', 'CSF1', 'CFD', 'C5'), cols  = c('red','green'))
RidgePlot(data, features = c('HLA-C', 'ITGA4', 'IL1B', 'IL1RN', 'IL10', 
                           'IL18BP', 'TGFBI', 'TNFRSF4', 'CXCL12', 'CSF1', 'CFD', 'C5'), cols  = c('red','green'))
VlnPlot(data, features = c('CD74'), cols  = c('red','green'))
RidgePlot(data, features = c('CD74'), cols  = c('red','green'))

