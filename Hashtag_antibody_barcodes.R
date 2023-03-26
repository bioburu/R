library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
setwd("/home/amp_prog/Downloads/GSE190964_RAW")
features_path <- "GSM5879885_features.tsv.gz"
barcodes_path <- "GSM5879885_barcodes.tsv.gz"
matrix_path <- "GSM5879885_matrix.mtx.gz"
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
dim(matrix)
str(matrix)
colnames(matrix)[1:10]
row.names(matrix)[1:10]
head(matrix[c('Hashtag1','Hashtag2','Hashtag6','Hashtag7','Hashtag8'),])
HTO<-matrix[c('Hashtag1','Hashtag2','Hashtag6','Hashtag7','Hashtag8'),]
HTO
#--------------------------------------------------------------------------------
data<-CreateSeuratObject(counts = matrix)
data<-NormalizeData(data)
data<-FindVariableFeatures(data)
data<-ScaleData(data,features = VariableFeatures(data))
data[['HTO']]<-CreateAssayObject(counts = HTO)
data<-NormalizeData(data, assay = 'HTO',normalization.method = 'CLR')
data<-HTODemux(data,assay = 'HTO',positive.quantile = 0.99)
table(data$HTO_classification.global)
VlnPlot(data,features = 'nCount_RNA',log = TRUE)
data<-subset(data,idents=c('Negative','Doublet'),invert=TRUE)
VlnPlot(data,features = 'nCount_RNA',log = TRUE)
head(data@meta.data)
head(data@active.ident)
head(data$nCount_RNA)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt <5)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
top10 <- head(VariableFeatures(data), 10)
top10
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
gc()
data
data <- RunPCA(data,npcs = 50, features = VariableFeatures(object = data))
#------------------------------------------------------------------
data<-RenameIdents(data,'Hashtag8'='untransduced','Hashtag2'='Transduced',
                   'Hashtag6'='Transduced','Hashtag1'='Transduced','Hashtag7'='untransduced')
head(data@meta.data)
dim(data)
cat('Get expression matrix from Seurat object')
#DF <- as.data.frame(as.matrix(GetAssayData(data)))
#Genes <- row.names(DF)
#DF <- cbind(Genes,DF)
#View(DF)
break
DimPlot(data, reduction = "pca", dims = c(1,2) ,pt.size = 3, label.box = TRUE, label.size = 2)
DimPlot(data, reduction = "pca", dims = c(3,4) ,pt.size = 3, label.box = TRUE, label.size = 2)

ProjectDim(data, reduction = "pca",dims.print = c(8))
data@meta.data
VlnPlot(data, features = c('HMGB1','HMGB2','HMGB3','TP53','HSP90AA1','HSPB1','HSPA1A','HSF1','HSPA4'))
VlnPlot(data, features = c('HSPA5','HSP90AB1','HSPA1B','HSPD1','HSPD1','HSPH1')) 
VlnPlot(data, features = c('MRPL1','SAAL1','CARD8','CARD9','SAP130','SCARB1')) 
VlnPlot(data, features = c('SEC14L1','TRIM25','TFAM')) 

#---damage associated molecular pattern (alarmins)-----------------------
VlnPlot(data, features = c('HSP90AA1','HSPB1','HSPA1A','HSF1','HSF2','HSPA4')) 
VlnPlot(data, features = c('HSPA5','HSP90AB1','HSPA1B','HSPD1','HSPD1','HSPH1')) 
VlnPlot(data, features = c('S100A1','S100A10','S100A11','S100A13','S100A2','S100A3')) 
VlnPlot(data, features = c('S100A4','S100A5','S100A6','S100A7A','S100A8','S100A9')) 
VlnPlot(data, features = c('S100B','S100G','S100P','S100PBP','S100Z','MRPL1')) 
VlnPlot(data, features = c('SAAL1','AGER','SYK','CARD6','CARD8','CARD9')) 
VlnPlot(data, features = c('CARD11','CARD14','LGALS1','SAP130','TREM1','TREML2')) 
VlnPlot(data, features = c('IFIH1','DDX58','DDX12P','LY96','PTAFR','SCARB1')) 
VlnPlot(data, features = c('SCARB2','MAVS','RNF135','SEC14L1','TRIM25','TFAM')) 
VlnPlot(data, features = c(),pt.size = 2)  
#------------finals----------------------------------------------------------------

VlnPlot(data, features = c('FOS', 'TNFAIP3', 'CXCR4', 'DUSP1', 'TSC22D3', 'PPP1R15A'),pt.size = 2)
VlnPlot(data, features = c('TOP2A', 'PTTG1','UBE2T', 'H2AX','CDKN3', 'CDK1'),pt.size = 2)
VlnPlot(data, features = c('CENPF', 'BIRC5', 'CCNA2', 'TK1','UBE2S','PDE6G'),pt.size = 2)
VlnPlot(data, features = c('LRRN3','IGFBP2', 'CD248','ITM2C','DBN1','CD38'),pt.size = 2)
VlnPlot(data, features = c('MMP25','IL9R', 'NDFIP2'),pt.size = 2)
#--------------PC1 negative component-------------------------------------------
VlnPlot(data, features = c('TYMS', 'STMN1', 'MKI67', 'RRM2', 'UBE2C','TUBA1B'))
VlnPlot(data, features = c('CKS2', 'TUBB', 'TOP2A', 'PTTG1','UBE2T', 'H2AX'))
VlnPlot(data, features = c('CDKN3', 'CDK1', 'CENPF', 'BIRC5', 'CCNA2', 'TK1'))
VlnPlot(data, features = c('HMGB2', 'UBE2S'))
#--------------PC8 positive component-------------------------------------------
VlnPlot(data, features = c('PDE6G', 'LRRN3', 'ACTB', 'IGFBP2', 'CD248', 'GZMA'))
VlnPlot(data, features = c('CD8A', 'GZMB', 'COTL1', 'ITM2C','IFITM2','DBN1'))
VlnPlot(data, features = c('CD38', 'PKM', 'IFITM3', 'C9orf135', 'MMP25', 'ZYX'))
VlnPlot(data, features = c('IL9R', 'NDFIP2')) 
#----------------------------------------------------------
#-------General T cell markers----------------------------------------------
VlnPlot(data, features = c('CD3D','CD3G','CD3E','CD4','CD8A','CD8B'))
#-------Idunno markers-------------------------------------------------
VlnPlot(data, features = c('PTPRC','CCR7','FOXP3','LCK','CD28','LAP3'))
VlnPlot(data, features = c('ZAP70','TRAC','TRBC1','TRBC2','PLCG1','ITK'))
VlnPlot(data, features = c('FYN','ADA','','CD69','CD2','CD40LG'))
VlnPlot(data, features = c('FAS','FASLG','VAV1','CD6','FLT3LG','CTLA4'))
#---------T cell receptor alpha variables-------------------------------------
VlnPlot(data, features = c('TRAV1.1','TRAV1-2','TRAV10','TRAV12-1','TRAV12-2','TRAV12-3'))
VlnPlot(data, features = c('TRAV13-1','TRAV13-2','TRAV14DV4','TRAV16','TRAV17','TRAV18'))


VlnPlot(data, features = c(), pt.size = 1)
VlnPlot(data, features = c('TRAV19','TRAV2','TRAV20','TRAV21','TRAV22','TRAV23DV6'), pt.size = 1)
VlnPlot(data, features = c('TRAV24','TRAV25','TRAV26.1','TRAV26.2','TRAV27','TRAV29DV5'), pt.size = 1)
VlnPlot(data, features = c('TRAV3','TRAV30','TRAV34','TRAV36DV7','TRAV38.1','TRAV38.2DV8'), pt.size = 1)
VlnPlot(data, features = c('TRAV39','TRAV4','TRAV40','TRAV41','TRAV5','TRAV6'), pt.size = 1)
VlnPlot(data, features = c('TRAV8.1','TRAV8.2','TRAV8.3','TRAV8.4','TRAV8.5','TRAV9.2'), pt.size = 1)
#--------T cell receptor beta variables-------------------------------------------
VlnPlot(data, features = c('TRBV10.1','TRBV10.2','TRBV10.3','TRBV11.1','TRBV11.2','TRBV11.3'), pt.size = 1)
VlnPlot(data, features = c('TRBV12.3','TRBV12.4','TRBV12.5','TRBV13','TRBV14','TRBV15'), pt.size = 1)
#metabolic activity
VlnPlot(data, features = c("GAPDH","ACACA","IDH2","HK1","G6PD","PRDX2"))
VlnPlot(data, features = c("PRDX2","ATP1A1","ATP1B1","CPT1A"))
#activation/exhaustion markers
VlnPlot(data, features = c("JAK2","STAT3","IRAK1","IRAK3","MAPK1",
                           "POGLUT1"))
VlnPlot(data, features = c("LAG3","CCL22","CXCL4","SH2D1A",
                           "SLAMF7","SLAMF8"))
#effector markers
VlnPlot(data, features = c("GZMA","GZMB","PRF1","PSAP","LAMP2","LAG3"))

####cYTOKEINSSSSSS-------------------------
VlnPlot(data, features = c("IL2","IL2RA","IL3","IL4","IL4I1","IL4R"))
VlnPlot(data, features = c("IFNG","IL1RAP","IL6","IL7R","CXCL8","IL10"))
VlnPlot(data, features = c("IL10RA","IL10RB","IL10RB.AS1","IL11","IL12RB1","IL12RB2"))
VlnPlot(data, features = c("IL13","IL13RA1","IL15","IL15RA","IL16","IL17RA"))
VlnPlot(data, features = c("IL17RB","IL17RC","IL18","IL18BP","IL21R","IL23A",""))
VlnPlot(data, features = c("","IL27RA","IL32","TNF","TNFAIP1","TNFAIP3","TNFSF10"))
VlnPlot(data, features = c("IL31","CSF1","CSF2","LDHA","CD46","GIMAP1"))

VlnPlot(data, features = c("TLR4","TLR6","HLA-A","HLA-B","HLA-E","CD55"))

VariableFeatures(data)[1:50]
