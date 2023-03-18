library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
#-------------------------------------------------------------------------------
setwd('/home/amp_prog/Desktop/IL6/GSE134874_PRR_stims')
matrix <- read.delim('GSE134874_Human_umi_counts.tsv')
row.names(matrix) <- make.names(matrix$X,unique = TRUE)
matrix <-matrix[,-c(1)]
data <- CreateSeuratObject(counts = matrix)
#---------------------------------------------------------------------
data@meta.data
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top200 <- head(VariableFeatures(data), 200)
top200
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
data <- ScaleData(data, features = all.genes)
dim(data)
data <- RunPCA(data,npcs=50,features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
VlnPlot(data, features = 'IL6')
VlnPlot(data, features = 'CD14')
VlnPlot(data, features = 'ITGAM')
VlnPlot(data, features = 'IL1B')
VlnPlot(data, features = 'IL8')
VlnPlot(data, features = 'TNFAIP3')
VlnPlot(data, features = 'NFKB1')
VlnPlot(data, features = 'MRC1')
VlnPlot(data, features = 'IL10')
VlnPlot(data, features = 'IL18')
VlnPlot(data, features = 'IL15')
VlnPlot(data, features = 'CSF1')
VlnPlot(data, features = 'CSF2')
VlnPlot(data, features = 'IL17RA')
VlnPlot(data, features = 'IDH1')
VlnPlot(data, features = 'IDH2')
VlnPlot(data, features = 'JAK2')
VlnPlot(data, features = 'STAT3')
VlnPlot(data, features = 'MSR1')
VlnPlot(data, features = 'TLR1')
VlnPlot(data, features = 'TLR2')
VlnPlot(data, features = 'TLR4')
VlnPlot(data, features = 'TLR8')
VlnPlot(data, features = 'HLA.DRA')
VlnPlot(data, features = 'LILRB1')
VlnPlot(data, features = 'CIITA')
VlnPlot(data, features = 'CD55')
VlnPlot(data, features = 'ITGAX')
VlnPlot(data, features = 'IGFBP3')
VlnPlot(data, features = 'IGFBP4')
VlnPlot(data, features = 'PIK3CG')
VlnPlot(data, features = 'ANTXR2')
VlnPlot(data, features = 'MILR1')
VlnPlot(data, features = '')




DimPlot(data, reduction = "pca", pt.size = 7, label.box = TRUE, label.size = 8)
break
data$il6.groups <- 'il6.pos'
data$il6.groups[WhichCells(data, expression= IL6 < 0.7)] <- 'il6.neg'
Idents(data) <- 'il6.groups'
summary(data$il6.groups)
il6.markers <- FindMarkers(data, ident.1= 'il6.pos',min.pct = 0.25)
il6.markers[1:50,]
#-----------tester -------------------------------------------------------------
VlnPlot(data,features = c("IL6",'IL1B','IL8','TNFAIP3','NFKB1','MRC1','PTPRC','NLRP3','E2F1'))#, cols = c('green','red','blue'))
#-----Cell type markers---------------------------------------------------------
RidgePlot(data,features = c("CD3D","CD8A","CD8B","CD19","MS4A1",
                            "CD22","NCAM1","CD14","ITGAM","ITGAX",
                            "FCGR3A","PTPRC"), cols = c('grey','red','blue'))
#-----showing all differences in cytokines--------------------------------------
RidgePlot(data,features = c("IFNG","IFNAR1","IFNLR1","IFNAR2",
                            "IL1B","IL1RAP","IL1RN","IL2RA","IL2RB",
                            "IL2RG","IL4I1","IL4R"), cols = c('grey','red','blue'))
####cYTOKEINSSSSSS-------------------------
RidgePlot(data, features = c("IL6","IL6R","IL6ST","IL7","IL7R","CXCL8",
                           "IL10","IL10RA","IL10RB","IL11",
                           "IL12RB1"), cols  = c('grey','red','blue'))
RidgePlot(data, features = c("IL12RB2","IL13","IL13RA1","IL15","IL15RA","IL16","IL17RA",
                             "IL17RB","IL17RC","IL18","IL18BP","IL21R"), cols  = c('grey','red'))
RidgePlot(data, features = c("IL23A","IL24","IL27","IL27RA","IL32","TNF",
                           "TNFAIP1","TNFAIP3","TNFSF10","CSF1","CSF1R",
                           "CSF2","CSFRA","CSFRB"), cols  = c('grey','red'))
RidgePlot(data, features = c("GAPDH","ACACA","IDH2","HK1","G6PD","PRDX2",
                             "ATP1A1","ATP1B1","CPT1A"), cols  = c('grey','red','blue'))
RidgePlot(data, features = c("JAK2","STAT3","IRAK1","IRAK3","MAPK1",
                             "POGLUT1","LAG3","CCL22","CCL8","SH2D1A",
                             "SLAMF7","SLAMF8"), cols  = c('grey','red'))
RidgePlot(data, features = c("GZMA","GZMB","PRF1","PSAP","LAMP1","LAMP2",
                           "ATG7","MRC1","MSR1","NR1H3"), cols  = c('grey','red','blue'))
RidgePlot(data, features = c("TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7",
                             "TLR8","TLR9","TLR10"), cols  = c('grey','red','blue'))
RidgePlot(data, features = c("HLA.A","HLA.B","HLA.C","HLA.G","HLA.H","HLA.J","HLA.L",
                             "HLA.DRA","HLA.DRB1",
                             "LILRB1","RFX5","CIITA"), cols  = c('grey','red','blue'))
RidgePlot(data, features = c("CD55","CD59","CD93","C3","C5","CFP","VTN","CR1",
                             "CFH","CD46"), cols  = c('grey','red','blue'))

