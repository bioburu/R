library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
setwd('/home/amp_prog/rstudio/GSE120084_carM')
matrix <- read.delim('GSE120084_Raw_read_counts_PhenoMap.txt')
row.names(matrix)<- make.names(matrix$GeneID,unique = TRUE)
matrix <- matrix[,-1]
data <- CreateSeuratObject(counts = matrix)
dim(data)
data@meta.data
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
data <- RunPCA(data,npcs=18, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
data
DimPlot(data, reduction = "pca", pt.size = 7, label.box = TRUE, label.size = 8)
break
#----------------------------------------------------------
cat('General monocyte marker panel')
VlnPlot(data, features = c("CD14","FCGR3A","ITGAM","ITGAX","CD38","CD4",
                           "CD86","FCGR1A","CD68","CD163","MILR1","CD80"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
RidgePlot(data, features = c("CD14","FCGR3A","ITGAM","ITGAX","CD38","CD4",
                           "CD86","FCGR1A","CD68","CD163","MILR1","CD80"), cols  = c('blue','green','red','orange','yellow'))
#metabolic activity
VlnPlot(data, features = c("GAPDH","ACACA","IDH2","HK1","G6PD","PRDX2",
                           "ATP1A1","ATP1B1","CPT1A"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
RidgePlot(data, features = c("GAPDH","ACACA","IDH2","HK1","G6PD","PRDX2",
                           "ATP1A1","ATP1B1","CPT1A"), cols  = c('blue','green','red','orange','yellow'))
#activation/exhaustion markers
VlnPlot(data, features = c("JAK2","STAT3","IRAK1","IRAK3","MAPK1",
                           "POGLUT1","LAG3","CCL22","CCL8","SH2D1A",
                           "SLAMF7","SLAMF8"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
RidgePlot(data, features = c("JAK2","STAT3","IRAK1","IRAK3","MAPK1",
                             "POGLUT1","LAG3","CCL22","CCL8","SH2D1A",
                             "SLAMF7","SLAMF8"), cols  = c('blue','green','red','orange','yellow'))
#effector markers
VlnPlot(data, features = c("GZMA","GZMB","PRF1","PSAP","LAMP1","LAMP2",
                           "ATG7","MRC1","MSR1","NR1H3"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
RidgePlot(data, features = c("GZMA","GZMB","PRF1","PSAP","LAMP1","LAMP2",
                           "ATG7","MRC1","MSR1","NR1H3"), cols  = c('blue','green','red','orange','yellow'))
####cYTOKEINSSSSSS-------------------------
VlnPlot(data, features = c("IFNG","IFNAR1","IFNLR1","IFNAR2","IFNAR2.2",
                           "IL1B","IL1RAP","IL1RN","IL2RA","IL2RB",
                           "IL2RG","IL4I1","IL4R"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
RidgePlot(data, features = c("IFNG","IFNAR1","IFNLR1","IFNAR2","IFNAR2.2",
                           "IL1B","IL1RAP","IL1RN","IL2RA","IL2RB",
                           "IL2RG","IL4I1","IL4R"), cols  = c('blue','green','red','orange','yellow'))
VlnPlot(data, features = c("IL6","IL6R","IL6ST","IL7","IL7R","CXCL8",
                           "IL10","IL10RA","IL10RB","IL10RB.AS1","IL11",
                           "IL12RA","IL12RB1"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
RidgePlot(data, features = c("IL6","IL6R","IL6ST","IL7","IL7R","CXCL8",
                           "IL10","IL10RA","IL10RB","IL10RB.AS1","IL11",
                           "IL12RA","IL12RB1"), cols  = c('blue','green','red','orange','yellow'))
VlnPlot(data, features = c("IL12RB2","IL13","IL13RA1","IL15","IL15RA","IL16","IL17RA",
                           "IL17RB","IL17RC","IL18","IL18BP","IL21R"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
RidgePlot(data, features = c("IL12RB2","IL13","IL13RA1","IL15","IL15RA","IL16","IL17RA",
                           "IL17RB","IL17RC","IL18","IL18BP","IL21R"), cols  = c('blue','green','red','orange','yellow'))
VlnPlot(data, features = c("IL23A","IL24","IL27","IL27RA","IL32","TNF",
                           "TNFAIP1","TNFAIP3","TNFSF10","CSF1","CSF1R",
                           "CSF2","CSFRA","CSFRB"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
RidgePlot(data, features = c("IL23A","IL24","IL27","IL27RA","IL32","TNF",
                           "TNFAIP1","TNFAIP3","TNFSF10","CSF1","CSF1R",
                           "CSF2","CSFRA","CSFRB"), cols  = c('blue','green','red','orange','yellow'))
VlnPlot(data, features = c("TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7",
                           "TLR8","TLR9","TLR10"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
RidgePlot(data, features = c("TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7",
                           "TLR8","TLR9","TLR10"), cols  = c('blue','green','red','orange','yellow'))
VlnPlot(data, features = c("HLA.A","HLA.B","HLA.C","HLA.G","HLA.H","HLA.J","HLA.L",
                           "HLA.DRA","HLA.DRB1",
                           "LILRB1","RFX5","CIITA"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
RidgePlot(data, features = c("HLA.A","HLA.B","HLA.C","HLA.G","HLA.H","HLA.J","HLA.L",
                           "HLA.DRA","HLA.DRB1",
                           "LILRB1","RFX5","CIITA"), cols  = c('blue','green','red','orange','yellow'))
#-------------Complement pathway--------------------------------------------
VlnPlot(data, features = c("CD55","CD59","CD93","C3","C5","CFP","VTN","CR1",
                           "CFH","CD46"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
RidgePlot(data, features = c("CD55","CD59","CD93","C3","C5","CFP","VTN","CR1",
                           "CFH","CD46"), cols  = c('blue','green','red','orange','yellow'))
VariableFeatures(data)[1:200]
