library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
setwd('/home/amp_prog/rstudio/CAR_T')
matrix <- read.csv('GSE208052_Processed_data_Normalized_gene_expression.csv')
rownames(matrix) <- make.names(matrix$X, unique = TRUE)
matrix <- matrix[,-1]
colnames(matrix) <- c('4-1BB.stim.6hr_1','4-1BB.ctrl.6hr_1','CD28.stim.6hr_1','CD28.ctrl.6hr_1',
                      '4-1BB.stim.24hr_1','4-1BB.ctrl.24hr_1','CD28.stim.24hr_1','CD28.ctrl.24hr_1',
                      '4-1BB.stim.24hr_2','4-1BB.ctrl.24hr_2','CD28.stim.24hr_2','CD28.ctrl.24hr_2',
                      '4-1BB.stim.6hr_2','4-1BB.ctrl.6hr_2','CD28.stim.6hr_2','CD28.stim.6hr_2',
                      '4-1BB.stim.24hr_3','4-1BB.ctrl.24hr_3','CD28.stim.24hr_3','CD28.ctrl.24hr_3',
                      '4-1BB.stim.6hr_3','4-1BB.ctrl.6hr_3','CD28.stim.6hr_3','CD28.ctrl.6hr_3')
data <- CreateSeuratObject(counts = matrix)
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
data <- RunPCA(data,npcs=23, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
data
DimPlot(data, reduction = "pca", dims = c(1,2) ,pt.size = 7, label.box = TRUE, label.size = 8)
DimPlot(data, reduction = "pca", dims = c(3,4) ,pt.size = 7, label.box = TRUE, label.size = 8)

break
#----------------------------------------------------------
#-------General T cell markers----------------------------------------------
VlnPlot(data, features = c('CD3D','CD3G','CD3E','CD4','CD8A','CD8B'), pt.size = 1)
#-------Idunno markers-------------------------------------------------
VlnPlot(data, features = c('PTPRC','CCR7','FOXP3','LCK','CD28','LAP3'), pt.size = 1)
VlnPlot(data, features = c('ZAP70','TRAC','TRBC1','TRBC2','PLCG1','ITK'), pt.size = 1)
VlnPlot(data, features = c('FYN','ADA','RAG1','CD69','CD2','CD40LG'), pt.size = 1)
VlnPlot(data, features = c('FAS','FASLG','VAV1','CD6','FLT3LG','CTLA4'), pt.size = 1)

#---------T cell receptor alpha variables-------------------------------------
VlnPlot(data, features = c('TRAV1.1','TRAV1.2','TRAV10','TRAV12.1','TRAV12.2','TRAV12.3'), pt.size = 1)
VlnPlot(data, features = c('TRAV13.1','TRAV13.2','TRAV14DV4','TRAV16','TRAV17','TRAV18'), pt.size = 1)
VlnPlot(data, features = c('TRAV19','TRAV2','TRAV20','TRAV21','TRAV22','TRAV23DV6'), pt.size = 1)
VlnPlot(data, features = c('TRAV24','TRAV25','TRAV26.1','TRAV26.2','TRAV27','TRAV29DV5'), pt.size = 1)
VlnPlot(data, features = c('TRAV3','TRAV30','TRAV34','TRAV36DV7','TRAV38.1','TRAV38.2DV8'), pt.size = 1)
VlnPlot(data, features = c('TRAV39','TRAV4','TRAV40','TRAV41','TRAV5','TRAV6'), pt.size = 1)
VlnPlot(data, features = c('TRAV8.1','TRAV8.2','TRAV8.3','TRAV8.4','TRAV8.5','TRAV9.2'), pt.size = 1)
#--------T cell receptor beta variables-------------------------------------------
VlnPlot(data, features = c('TRBV10.1','TRBV10.2','TRBV10.3','TRBV11.1','TRBV11.2','TRBV11.3'), pt.size = 1)
VlnPlot(data, features = c('TRBV12.3','TRBV12.4','TRBV12.5','TRBV13','TRBV14','TRBV15'), pt.size = 1)



#metabolic activity
VlnPlot(data, features = c("GAPDH","ACACA","IDH2","HK1","G6PD","PRDX2"), pt.size = 1)
VlnPlot(data, features = c("PRDX2","ATP1A1","ATP1B1","CPT1A"), pt.size = 1)
#activation/exhaustion markers
VlnPlot(data, features = c("JAK2","STAT3","IRAK1","IRAK3","MAPK1",
                           "POGLUT1"), pt.size = 1)
VlnPlot(data, features = c("LAG3","CCL22","CXCL4","SH2D1A",
                           "SLAMF7","SLAMF8"), pt.size = 1)
#effector markers
VlnPlot(data, features = c("GZMA","GZMB","PRF1","PSAP","LAMP2","LAG3"), pt.size = 1)

####cYTOKEINSSSSSS-------------------------
VlnPlot(data, features = c("","IFNAR1","","IFNAR2",
                           "",""), pt.size = 1)
VlnPlot(data, features = c("IL2","IL2RA","IL3",
                           "IL4","IL4I1","IL4R"), pt.size = 1)

VlnPlot(data, features = c("IFNG","IL1RAP","IL6","IL7R","CXCL8",
                           "IL10"), pt.size = 1)
VlnPlot(data, features = c("IL10RA","IL10RB","IL10RB.AS1","IL11",
                           "IL12RB1","IL12RB2"), pt.size = 1)

VlnPlot(data, features = c("IL13","IL13RA1","IL15","IL15RA","IL16","IL17RA"), pt.size = 1)
VlnPlot(data, features = c("IL17RB","IL17RC","IL18","IL18BP","IL21R","IL23A","IL24"), pt.size = 1)

VlnPlot(data, features = c("IL27","IL27RA","IL32","TNF","TNFAIP1","TNFAIP3","TNFSF10"), pt.size = 1)
VlnPlot(data, features = c("IL31","CSF1","CSF2","LDHA","CD46","GIMAP1"), pt.size = 1)

VlnPlot(data, features = c("TLR1","TLR6","HLA.A","HLA.B","HLA.E","CD55"), pt.size = 1)

#-------------Complement pathway--------------------------------------------
VlnPlot(data, features = c("CD55","CD59","CD93","C3","C5"), pt.size = 1)
VlnPlot(data, features = c("CFP","VTN","CR1","CFH","CD46"), pt.size = 1)

VariableFeatures(data)[1:200]
