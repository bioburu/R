library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
setwd("/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts")
cat('Yescarta/CD28z anti-CD19 CAR-T cells')
#-------------Baseline samples into objects--------------------------------------------------
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient6-Baseline/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient6-Baseline/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient6-Baseline/matrix.mtx.gz"
b6 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
B6 <- CreateSeuratObject(counts = b6,project = 'Untransduced') 
gc()
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient8-Baseline/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient8-Baseline/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient8-Baseline/matrix.mtx.gz"
b8 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
B8 <- CreateSeuratObject(counts = b8,project = 'Untransduced') 
gc()
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient9-Baseline/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient9-Baseline/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient9-Baseline/matrix.mtx.gz"
b9 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
B9 <- CreateSeuratObject(counts = b9,project = 'Untransduced') 
gc()
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient10-Baseline/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient10-Baseline/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient10-Baseline/matrix.mtx.gz"
b10 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
B10 <- CreateSeuratObject(counts = b10,project = 'Untransduced') 
gc()
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient11-Baseline/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient11-Baseline/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient11-Baseline/matrix.mtx.gz"
b11 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
B11 <- CreateSeuratObject(counts = b11,project = 'Untransduced') 
gc()
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient12-Baseline/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient12-Baseline/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient12-Baseline/matrix.mtx.gz"
b12 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
B12 <- CreateSeuratObject(counts = b12,project = 'Untransduced') 
gc()
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient13-Baseline/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient13-Baseline/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient13-Baseline/matrix.mtx.gz"
b13 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
B13 <- CreateSeuratObject(counts = b13,project = 'Untransduced') 
gc()
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient14-Baseline/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient14-Baseline/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient14-Baseline/matrix.mtx.gz"
b14 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
B14 <- CreateSeuratObject(counts = b14,project = 'Untransduced') 
gc()
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient15-Baseline/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient15-Baseline/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient15-Baseline/matrix.mtx.gz"
b15 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
B15 <- CreateSeuratObject(counts = b15,project = 'Untransduced') 
gc()
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient17-Baseline/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient17-Baseline/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient17-Baseline/matrix.mtx.gz"
b17 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
B17 <- CreateSeuratObject(counts = b17,project = 'Untransduced') 
rm(b6,b8,b9,b10,b11,b12,b13,b14,b15,b17)
gc()
#-------------Infusion product into objects---------------------------------------
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient1-Infusion/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient1-Infusion/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient1-Infusion/matrix.mtx.gz"
inf1 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
INF1 <- CreateSeuratObject(counts = inf1,project = 'Transduced')
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient2-Infusion/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient2-Infusion/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient2-Infusion/matrix.mtx.gz"
inf2 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
INF2 <- CreateSeuratObject(counts = inf2,project = 'Transduced')
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient3-Infusion/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient3-Infusion/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient3-Infusion/matrix.mtx.gz"
inf3 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
INF3 <- CreateSeuratObject(counts = inf3,project = 'Transduced')
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient4-Infusion/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient4-Infusion/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient4-Infusion/matrix.mtx.gz"
inf4 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
INF4 <- CreateSeuratObject(counts = inf4,project = 'Transduced')
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient5-Infusion/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient5-Infusion/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient5-Infusion/matrix.mtx.gz"
inf5 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
INF5 <- CreateSeuratObject(counts = inf5,project = 'Transduced')
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient7-Infusion/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient7-Infusion/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient7-Infusion/matrix.mtx.gz"
inf7 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
INF7 <- CreateSeuratObject(counts = inf7,project = 'Transduced')
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient8-Infusion/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient8-Infusion/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient8-Infusion/matrix.mtx.gz"
inf8 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
INF8 <- CreateSeuratObject(counts = inf8,project = 'Transduced')
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient9-Infusion/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient9-Infusion/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient9-Infusion/matrix.mtx.gz"
inf9 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
INF9 <- CreateSeuratObject(counts = inf9,project = 'Transduced')
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient10-Infusion/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient10-Infusion/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient10-Infusion/matrix.mtx.gz"
inf10 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
INF10 <- CreateSeuratObject(counts = inf10,project = 'Transduced')
features_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient11-Infusion/features.tsv.gz"
barcodes_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient11-Infusion/barcodes.tsv.gz"
matrix_path <- "/home/amp_prog/Desktop/rstudio/CAR_T/GSE197268_3timepts/Patient11-Infusion/matrix.mtx.gz"
inf11 <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
INF11 <- CreateSeuratObject(counts = inf11,project = 'Transduced')
rm(inf1,inf2,inf3,inf4,inf5,inf7,inf8,inf9,inf10,inf11)
gc()

#----------------Merge all objects---------------------------------------------

data <- merge(B6, y = c(B8,B9,B10,B11,B12,INF1,INF2,INF3,INF4,INF5,INF7),
                                  add.cell.ids = c('Untransduced','Untransduced','Untransduced','Untransduced','Untransduced','Untransduced',
                                                  # 'Untransduced','Untransduced','Untransduced','Untransduced','Untransduced',
                                                   'Transduced','Transduced','Transduced','Transduced','Transduced','Transduced'))
                                                   #'Transduced','Transduced','Transduced','Transduced','Transduced'), project = "Tcell_test")
rm(B6,B8,B9,B10,B11,B12,B13,B14,B15,B17,INF1,INF2,INF3,INF4,INF5,INF7,INF8,INF9,INF10,INF11)
gc()
head(colnames(data))
table(data$orig.ident)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
RidgePlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt <5)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
gc()
top10 <- head(VariableFeatures(data), 10)
top10
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
gc()
data <- ScaleData(data, features = all.genes)
dim(data)
gc()
data <- RunPCA(data,npcs = 50, features = VariableFeatures(object = data))
#------------------------------------------------------------------
dim(data)
data$TRAC.groups <- 'TRAC.pos'
data$TRAC.groups[WhichCells(data, expression= TRAC < 0.5)] <- 'TRAC.neg'
DimPlot(data, reduction = 'pca',split.by = 'TRAC.groups')
data@meta.data
data <- subset(data, subset = TRAC.groups != "TRAC.neg")
dim(data)
#---------------------------------------------------------------------------------
data@meta.data
data <- RunPCA(data,npcs = 50, features = VariableFeatures(object = data))
Loadings(data[['pca']])
VariableFeatures(data)
cat('Get expression matrix from Seurat object')
#DF <- as.data.frame(as.matrix(GetAssayData(data)))
#Genes <- row.names(DF)
#DF <- cbind(Genes,DF)
#View(DF)
break
DimPlot(data, reduction = "pca", dims = c(1,4) ,pt.size = 0.5, label.box = TRUE, label.size = 3)
DimPlot(data, reduction = "pca", dims = c(7,8) ,pt.size = 0.5, label.box = TRUE, label.size = 3)
DimPlot(data, reduction = "pca", dims = c(13,14) ,pt.size = 0.5, label.box = TRUE, label.size = 3)
ProjectDim(data, reduction = "pca",dims.print = c(8))
#---damage associated molecular pattern (alarmins)-----------------------
data@meta.data
VlnPlot(data, features = c('HMGB1','HMGB2','VIM','TP53')) 
VlnPlot(data, features = c('HSP90AA1','HSPB1','HSPA1A','HSF1','HSF2','HSPA4')) 
VlnPlot(data, features = c('HSPA5','HSP90AB1','HSPA1B','HSPD1','HSPD1','HSPH1')) 
#------------finals----------------------------------------------------------------
VlnPlot(data, features = c('TYMS', 'STMN1', 'MKI67', 'RRM2', 'UBE2C','CKS2'))
VlnPlot(data, features = c('TOP2A', 'PTTG1','UBE2T', 'H2AX','CDKN3', 'CDK1'))
VlnPlot(data, features = c('CENPF', 'BIRC5', 'CCNA2', 'TK1','UBE2S','PDE6G'))
VlnPlot(data, features = c('LRRN3','IGFBP2', 'CD248','ITM2C','DBN1','CD38'))
VlnPlot(data, features = c('MMP25','IL9R', 'NDFIP2'))
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

break
#----------------------------------------------------------
#-------General T cell markers----------------------------------------------
VlnPlot(data, features = c('CD3D','CD3G','CD3E','CD4','CD8A','CD8B'))
#-------Idunno markers-------------------------------------------------
VlnPlot(data, features = c('PTPRC','CCR7','FOXP3','LCK','CD28','LAP3'))
VlnPlot(data, features = c('ZAP70','TRAC','TRBC1','TRBC2','PLCG1','ITK'))
VlnPlot(data, features = c('FYN','ADA','','CD69','CD2','CD40LG'))
VlnPlot(data, features = c('FAS','FASLG','VAV1','CD6','FLT3LG','CTLA4'))

#---------T cell receptor alpha variables-------------------------------------
VlnPlot(data, features = c('TRAV1-1','TRAV1-2','TRAV10','TRAV12-1','TRAV12-2','TRAV12-3'))
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
