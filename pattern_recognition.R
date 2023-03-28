library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
setwd("/home/amp_prog/Desktop/rstudio/CAR_T/GSE179249_PAMP")
cat('DAMP_CARs')
#-------------RN7SL1+CARs---------------------------------------------------
features_path1 <- "GSM5412146_M5_7SL_1_features.tsv.gz"
barcodes_path1 <- "GSM5412146_M5_7SL_1_barcodes.tsv.gz"
matrix_path1 <- "GSM5412146_M5_7SL_1_matrix.mtx.gz"
exp1 <- ReadMtx(mtx= matrix_path1, features = features_path1, cells= barcodes_path1)
M5BBz_7SL_1 <- CreateSeuratObject(counts = exp1, project = 'RN7SL1+_TME') 
rm(exp1)
gc()
features_path2 <- "GSM5412147_M5_7SL_2_features.tsv.gz"
barcodes_path2 <- "GSM5412147_M5_7SL_2_barcodes.tsv.gz"
matrix_path2 <- "GSM5412147_M5_7SL_2_matrix.mtx.gz"
exp2 <- ReadMtx(mtx= matrix_path2, features = features_path2, cells= barcodes_path2)
M5BBz_7SL_2 <- CreateSeuratObject(counts = exp2, project = 'RN7SL1+_TME') 
rm(exp2)
gc()
#----------------------------------------------------------------------------------
#-------------control CARs---------------------------------------------------------------
features_path3 <- "GSM5412148_M5BBz_1_features.tsv.gz"
barcodes_path3 <- "GSM5412148_M5BBz_1_barcodes.tsv.gz"
matrix_path3 <- "GSM5412148_M5BBz_1_matrix.mtx.gz"
ctrl1 <- ReadMtx(mtx= matrix_path3, features = features_path3, cells= barcodes_path3)
M5BBz_1 <- CreateSeuratObject(counts = ctrl1, project = 'RN7SL1-_TME') 
rm(ctrl1)
gc()
features_path4 <- "GSM5412149_M5BBz_2_features.tsv.gz"
barcodes_path4 <- "GSM5412149_M5BBz_2_barcodes.tsv.gz"
matrix_path4 <- "GSM5412149_M5BBz_2_matrix.mtx.gz"
ctrl2 <- ReadMtx(mtx= matrix_path4, features = features_path4, cells= barcodes_path4)
M5BBz_2 <- CreateSeuratObject(counts = ctrl2, project = 'RN7SL1-_TME') 
rm(ctrl2)
gc()
#----------------Merge all objects---------------------------------------------
data <- merge(M5BBz_1, y = c(M5BBz_2,M5BBz_7SL_1,M5BBz_7SL_2),
                                  add.cell.ids = c('M5BBz','M5BBz','M5BBz_7SL','M5BBz_7SL'), 
                                  project = 'DAMP_CARs')
head(colnames(data))
table(data$orig.ident)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt <5)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
gc()
top20 <- head(VariableFeatures(data), 20)
top20
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
data <- ScaleData(data, features = all.genes)
dim(data)
gc()
data <- RunPCA(data, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
#---------------------------------------------------------------------------------
VlnPlot(data, features = c('CD14'),pt.size = 2)
dim(data)
data$CD14.groups <- 'CD14.pos'
data$CD14.groups[WhichCells(data, expression= CD14 < 0.5)] <- 'CD14.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD14.groups')
data@meta.data
data <- subset(data, subset = CD14.groups != "CD14.neg")
dim(data)
#---------------------------------------------------------------------------------
data@meta.data
data <- RunPCA(data,npcs = 8, features = VariableFeatures(object = data))
Loadings(data[['pca']])
VariableFeatures(data)
cat('Get expression matrix from Seurat object')
DF <- as.data.frame(as.matrix(GetAssayData(data)))
Genes <- row.names(DF)
DF <- cbind(Genes,DF)
View(DF)
break
#-----------------------------------------------------------------------------
DimPlot(data, reduction = "pca", dims = c(1,3) ,pt.size = 2, label.box = TRUE, label.size = 3,cols = c('green','red'))
DimPlot(data, reduction = "pca", dims = c(3,7) ,pt.size = 2, label.box = TRUE, label.size = 3,cols = c('green','red'))
#--------------TME CD14+ genes--------------------------------------------
#-------------Cytoprotective autophagy phenotype---------------------- 
VlnPlot(data, features = c('VAMP8','CD59','HLA-A','HMGB1','HSPB1','CD68','DDX58','IFIH1','MILR1'),pt.size = 3)
VlnPlot(data, features = c(),pt.size = 3)


VlnPlot(data,features = c('CD55','CD59','CD68','NR1H3','IL7','HLA-A'), pt.size = 2)
VlnPlot(data, features = c('HLA-B','HLA-C','MKI67','NFKB1','IL18'))
#PC_3 positives------Cell division/mitosis DAVID--------------------------------------------------------------
VlnPlot(data, features = c('PRODH2', 'PBK', 'NEK2', 'GPD1', 'HELLPAR', 'POLQ'),pt.size = 0.5)
VlnPlot(data, features = c('NREP', 'SERPINA1', 'IL6R', 'PKMYT1', 'EFCAB2', 'NCAPH'),pt.size = 0.5)
VlnPlot(data, features = c('NR2F1', 'WASF3', 'AC010359.3', 'CCDC18', 'KIF18A', 'HIRIP3'),pt.size = 0.5)
VlnPlot(data, features = c('BIRC5', 'CDCA3', 'AC007364.1', 'GOLT1A', 'GSDMB', 'DNAJC27'),pt.size = 0.5)
VlnPlot(data, features = c('CDK1', 'GALT', 'IFI44L', 'AKR1C1', 'CGN', 'VIL1'),pt.size = 0.5)
#-------------Variable genes----------------------------------------------------------
VlnPlot(data, features = c('CD14','CD68','MILR1','GAPDH','HK1','PRDX2',
                             'ATP1A1','ATP1B1','CPT1A','STAT1','IRAK1','PSAP'),pt.size = 0.5)
VlnPlot(data, features = c('LAMP1','LAMP2','NR1H3','IFNAR1','IFNAR2','IL2RG',
                             'IL4R','IL6R','IL7','IL10RB','IL12RB1','IL13RA1'),pt.size = 0.5)
VlnPlot(data, features = c('IL17RC','IL18','IL18BP','IL32','TNFSF10','HLA-A',
                             'HLA-B','HLA-C','CD55','CD59','CD46','PTPRC'),pt.size = 0.5)
#-------------monocyte markers---------------------------------------------------
VlnPlot(data, features = c("CD14","FCGR3A","ITGAM","ITGAX","CD38","CD4",
                           "CD86","FCGR1A","CD68","CD163","MILR1","CD80"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
#metabolic activity
VlnPlot(data, features = c("GAPDH","ACACA","IDH2","HK1","G6PD","PRDX2",
                           "ATP1A1","ATP1B1","CPT1A"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
#activation/exhaustion markers
VlnPlot(data, features = c("JAK2","STAT3","IRAK1","IRAK3","MAPK1",
                           "POGLUT1","LAG3","CCL22","CCL8","SH2D1A",
                           "SLAMF7","SLAMF8"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
#effector markers
VlnPlot(data, features = c("GZMA","GZMB","PRF1","PSAP","LAMP1","LAMP2",
                           "ATG7","MRC1","MSR1","NR1H3"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
####cYTOKEINSSSSSS-------------------------
VlnPlot(data, features = c("IFNG","IFNAR1","IFNLR1","IFNAR2","IFNAR2.2",
                           "IL1B","IL1RAP","IL1RN","IL2RA","IL2RB",
                           "IL2RG","IL4I1","IL4R"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
VlnPlot(data, features = c("IL6","IL6R","IL6ST","IL7","IL7R","CXCL8",
                           "IL10","IL10RA","IL10RB","IL10RB.AS1","IL11",
                           "IL12RA","IL12RB1"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
VlnPlot(data, features = c("IL12RB2","IL13","IL13RA1","IL15","IL15RA","IL16","IL17RA",
                           "IL17RB","IL17RC","IL18","IL18BP","IL21R"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
VlnPlot(data, features = c("IL23A","IL24","IL27","IL27RA","IL32","TNF",
                           "TNFAIP1","TNFAIP3","TNFSF10","CSF1","CSF1R",
                           "CSF2","CSFRA","CSFRB"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
VlnPlot(data, features = c("TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7",
                           "TLR8","TLR9","TLR10"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
VlnPlot(data, features = c("HLA-A","HLA-B","HLA-C","HLA-G","HLA-H","HLA-J","HLA-L",
                           "HLA-DRA","HLA-DRB1",
                           "LILRB1","RFX5","CIITA"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))
#-------------Complement pathway--------------------------------------------
VlnPlot(data, features = c("CD55","CD59","CD93","C3","C5","CFP","VTN","CR1",
                           "CFH","CD46"), pt.size = 1, cols  = c('blue','green','red','orange','yellow'))

#---damage associated molecular pattern (alarmins)-----------------------
data@meta.data
VlnPlot(data, features = c('HMGB1','HMGB2','HMGB3','VIM','TP53'),pt.size = 2) 
VlnPlot(data, features = c('HSP90AA1','HSPB1','HSPA1A','HSF1','HSF2','HSPA4'),pt.size = 2) 
VlnPlot(data, features = c('HSPA5','HSP90AB1','HSPA1B','HSPD1','HSPD1','HSPH1'),pt.size = 2) 
VlnPlot(data, features = c('S100A1','S100A10','S100A11','S100A13','S100A2','S100A3'),pt.size = 2) 
VlnPlot(data, features = c('S100A4','S100A5','S100A6','S100A7A','S100A8','S100A9'),pt.size = 2) 
VlnPlot(data, features = c('S100B','S100G','S100P','S100PBP','S100Z','MRPL1'),pt.size = 2) 
VlnPlot(data, features = c('SAAL1','AGER','SYK','CARD6','CARD8','CARD9'),pt.size = 2) 
VlnPlot(data, features = c('CARD11','CARD14','LGALS1','SAP130','TREM1','TREML2'),pt.size = 2) 
VlnPlot(data, features = c('IFIH1','DDX58','DDX12P','LY96','PTAFR','SCARB1'),pt.size = 2) 
VlnPlot(data, features = c('SCARB2','MAVS','RNF135','SEC14L1','TRIM25','TFAM'),pt.size = 2) 
VlnPlot(data, features = c(),pt.size = 2) 
#-------------phagosome genes----------------------------------
VlnPlot(data, features = c('VAMP3',"VAMP5","VAMP8",'MR1','CORO1C','TLR8'),pt.size = 2) 
VlnPlot(data, features = c('SEC22B','FCAR','FCGRT','STX18','DCSTAMP','MARCO'),pt.size = 2) 
VlnPlot(data, features = c('CD36','RAB10','RAB5A','LAMP3','LAMP5','M6PR'),pt.size = 2) 
VlnPlot(data, features=c('PIKFYVE','TAP1','TAP2','SEC61A1','SEC61B','SEC61G'),pt.size = 2)
#----------------RIG1--------------------------------------------
VlnPlot(data, features = c('DDX58','IFIH1','IRF7','TRIM25','CCL8','TAP1'),pt.size = 2) 
#-----------------il6---------------------------
VlnPlot(data, features = c('IL6','CD55','LAMP3','VAMP8','MRC1'),pt.size = 2) 
VlnPlot(data, features = c('CD14','MR1'),pt.size = 2) 
