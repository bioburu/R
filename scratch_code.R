library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
setwd('/home/amp_prog/rstudio/HLH')
matrix <- read.csv('external_gene_name.csv')
rownames(matrix) <- matrix$X
matrix <- matrix[,-1]
colnames(matrix) <- c('3.4yr_1','3.4yr_2','3.4yr_3','3.4yr_4',
            '4.3yr_1','4.3yr_2','4.3yr_3','4.3yr_4',
            '12yr_1','12yr_2','12yr_3',
            'Norm_1','Norm_2','Norm_3','Norm_4','Norm_5','Norm_6')
head(matrix)
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
data <- RunPCA(data,npcs=16, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:6, cells = 500, balanced = T)
ElbowPlot(data)
data
DimPlot(data, reduction = "pca", pt.size = 7, label.box = TRUE, label.size = 8)
break

#------------------LYMPHOCYTE PANEL----------------------------------------
VlnPlot(data, features = c("CD3D","CD19","NCAM1","CD8A","IL2", "CCR7", "LAMP1","CD3E","IL10","IGHA1",
                           "CXCR4","IL3RA","IL7R","IL6","IL18","CD14"), pt.size = 1, cols  = c('blue','green','red','orange'))

VlnPlot(data, features = c("CD14","CD22","CD8B","ITGAX","CD38","FCGR3A","ITGAM", "MILR1", "FLT3LG","IL15","IDH2","GAPDH","JAK2",
                           "CCL3","CD4","LAG3"), pt.size = 1, cols  = c('blue','green','red','orange'))
####cYTOKEINSSSSSS-------------------------
VlnPlot(data, features = c("IL1RAP","IL13RA2","IL1RN","IL12RB2","IL11",
                           "IL19","IL18","IL17RB","IL10RB","IL10","IL1RAP",
                           "IL13RA2","IL2RB","IL23A","IL27RA","IL2RA"), pt.size = 1, cols  = c('blue','green','red','orange'))
VlnPlot(data, features = c("IL22RA1","IL2RG","IL32","IL4R","IL4","IL5","IL6",
                           "IL7","CSF3","CSF1R","TNFRSF17","TNFAIP1","TNFAIP3",
                           "TNFSF10","TNFAIP6","TNFSF9"), pt.size = 1, cols  = c('blue','green','red','orange'))
VlnPlot(data, features = c("TNFSF14","TNFRSF11A","TNFAIP8","TNFRSF10B","TNFSF4",
                           "TNFSF12","IFNG","IFNAR1","IFNLR1","IFNAR2","IFNAR2.2",
                           "IL6.3"), pt.size = 1, cols  = c('blue','green','red','orange'))
#-------------Complement pathway--------------------------------------------
VlnPlot(data, features = c("CD93","CD55","CD59","C2","C5","C4A","C8G","CFP",
                           "CLU","C1QA","CR1","CR2","CFI","ITGAM","ITGAX","C1QB"), pt.size = 1, cols  = c('blue','green','red','orange'))
VlnPlot(data, features = c("C1R","VTN","CD2","C5AR2","CFH","C5AR1","LAT"), pt.size = 1, cols  = c('blue','green','red','orange'))
VariableFeatures(data)[1:200]
#--------------------------------------------------------------------
VlnPlot(data, features = c("IL6","IL6.1","IL6.2","IL6.3","IL6.4",
                           "IL6.5","IL6.7"), pt.size = 1, cols  = c('blue','green','red','orange'))
