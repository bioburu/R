library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(caTools)
library(car)
library(caret)
library(InformationValue)
library(pROC)
library(ROCR)
setwd('/home/deviancedev/Desktop/drive_nov2023/FCCC/GSE207921_h9.nsc_scRNAseq')
features_path <- 'GSM6322990_Stage_I_Exp1_features.tsv.gz'
barcodes_path <- 'GSM6322990_Stage_I_Exp1_barcodes.tsv.gz'
matrix_path <- 'GSM6322990_Stage_I_Exp1_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
x <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'Day13')
summary(x@active.ident)
x<-subset(x = x, downsample = 3160)
summary(x@active.ident)
#---------------------------------------------------------------------------
setwd('/home/deviancedev/Desktop/drive_nov2023/FCCC/GSE207921_h9.nsc_scRNAseq')
features_path <- 'GSM6322991_Stage_II_Exp1_features.tsv.gz'
barcodes_path <- 'GSM6322991_Stage_II_Exp1_barcodes.tsv.gz'
matrix_path <- 'GSM6322991_Stage_II_Exp1_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
x1 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'Day21')
summary(x1@active.ident)
x1<-subset(x = x1, downsample = 3160)
summary(x1@active.ident)
#---------------------------------------------------------------------------
setwd('/home/deviancedev/Desktop/drive_nov2023/FCCC/GSE207921_h9.nsc_scRNAseq')
features_path <- 'GSM6322993_Stage_III_Exp2_features.tsv.gz'
barcodes_path <- 'GSM6322993_Stage_III_Exp2_barcodes.tsv.gz'
matrix_path <- 'GSM6322993_Stage_III_Exp2_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
x2 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'Day30')
summary(x2@active.ident)
x2<-subset(x = x2, downsample = 3160)
summary(x2@active.ident)
#-------------------------------------------------------------------------------
setwd('/home/deviancedev/Desktop/drive_nov2023/FCCC/GSE207921_h9.nsc_scRNAseq')
features_path <- 'GSM6322995_Stage_IV_Exp3_features.tsv.gz'
barcodes_path <- 'GSM6322995_Stage_IV_Exp3_barcodes.tsv.gz'
matrix_path <- 'GSM6322995_Stage_IV_Exp3_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
x3 <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'Day45')
summary(x3@active.ident)
x3<-subset(x = x3, downsample = 3160)
summary(x3@active.ident)
#-------------------------------------------------------------------------------
data<-merge(x,y=c(x1,x2,x3),project='h9.esc_timepoints')
table(data@meta.data$orig.ident)
head(data@active.ident)
rm(x,x1,x2,x3,matrix)
gc()
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 11000 & percent.mt <15)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
gc()
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top1000 <- head(VariableFeatures(data), 1000)
top1000
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
data <- ScaleData(data, features = all.genes)
gc()
dim(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
data <- RunUMAP(data, dims = 1:30)
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
gc()
DimPlot(data, reduction = "pca",dims = c(1,2))
break
FeatureScatter(data, feature1 = "CMTM8", feature2 = "PCDH9",cols = c('grey','black','skyblue','red'))
DimPlot(data, reduction = 'umap',cols = c('grey','black','skyblue','red'))
cat('One day after passaging, the culture medium of hESCs was changed to the neural induction medium (NIM) that containing DMEM/F12, 1×N2 supplement (Gibco), and 1×NEAA, with the addition of SHH (C25II, R&D Systems, 500 ng/ml), CHIR99021 (Tocris, 0.4μM), DMH-1 (Tocris, 2 μM), and SB431542 (Stemgent, 2 μM) for 8 days (Day1 to Day9). On day 9, singular colonies were blown off gently and transferred to 6-well dish coated with fresh MEFs in NIM with the addition of SHH (100 ng/ml), SAG (Millipore, 1 μM), and CHIR99021 (0.4μM) for 4 days until stage I (Day 13). At stage I, singular colonies were blown off gently and transferred to non-adherent 25 ml dish cultured in suspension medium containing NIM with the addition of SHH (20 ng/ml), SAG (Millipore, 0.5 μM), and FGF8b (PeproTech, 100 ng/ml) for 8 days until stage II (Day21). From stage II to stage III (Day 30), the neuro-spheres were allowed to continue differentiating and proliferating in the suspension medium containing NIM with the addition of SHH (20 ng/ml), and FGF8b (20 ng/ml). At stage III, the neuro-spheres were dissociated by Accutase (Innovative Cell Technologies) at 37°C for 6 minutes and then replated onto the 24-well dish coated with Matrigel (BD Biosciences). From stage III to stage IV (Day 45), the differentiated cells were fed on neural differentiation medium (NDM) containing neurobasal medium, 1×N2 supplement (Gibco), and 1×B27 (Life Technologies) with the addition of brain-derived neurotrophic factor (BDNF, Peprotech, 10 ng/ml), glial-derived neurotrophic factor (GDNF, peprotech, 10 ng/ml), transforming growth factorβ3 (TGFβ3, R&D Systems, 1 ng/ml), ascorbic acid (AA, Sigma-Aldrich, 200 μM),cAMP(Sigma-Aldrich, 1 μM), and Compound E (Calbiochem, 1 μM).')
FeaturePlot(data, features = c('PRTG','CMTM8','MKI67','PCDH9'),reduction = 'umap',cols = c('grey','red'))
VlnPlot(data, features = c('PRTG','CMTM8','MKI67','PCDH9'),cols = c('grey','black','skyblue','red'))
RidgePlot(data, feature = c('PRTG','CMTM8','MKI67','PCDH9'),cols = c('grey','black','skyblue','red'))
DoHeatmap(data, features = c('PRTG','CMTM8','MKI67','PCDH9'))
data<-JoinLayers(data)
day13 <- FindMarkers(data, ident.1 = 'Day13',logfc.threshold = 1,min.pct = 0.4,test.use = 'bimod')
View(day13)
#write.csv(day13,file = 'GSE207921_d13.csv')
day21 <- FindMarkers(data, ident.1 = 'Day21',logfc.threshold = 1,min.pct = 0.4,test.use = 'bimod')
View(day21)
#write.csv(day21,file = 'GSE207921_d21.csv')
day30 <- FindMarkers(data, ident.1 = 'Day30',logfc.threshold = 1,min.pct = 0.4,test.use = 'bimod')
View(day30)
#write.csv(day30,file = 'GSE207921_d30.csv')
day45 <- FindMarkers(data, ident.1 = 'Day45',logfc.threshold = 1,min.pct = 0.4,test.use = 'bimod')
View(day45)
#write.csv(day45,file = 'GSE207921_d45.csv')
#-------------------------------------------------------------------------------
#---If clustering 
#devtools::install_github('immunogenomics/presto')
data <- FindNeighbors(data, dims = 1:30)
data <- FindClusters(data)
cluster1.markers <- FindMarkers(data, ident.1 = 1)
head(cluster1.markers,n=200)
#-------------------------------------------------------------------------------
#-------Paired cluster comparisons
x<-FindMarkers(data, ident.1 = '2', ident.2 = '1', 
               features = c(top1000),logfc.threshold=1,min.pct1=1,
               max.pct2=0.0001,only.pos = TRUE)
VlnPlot(data, features = c(row.names(x)[1:12]),cols = c('grey','red'),idents = c(1,2))
break 
#----------Isolate CD45+  -------------------------------------
data$CD45.groups <- 'CD45.pos'
data$CD45.groups[WhichCells(data, expression= PTPRC < 0.1)] <- 'CD45.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD45.groups')
head(data@meta.data)
data <- subset(data, subset = CD45.groups != "CD45.neg")
gc()
table(data@meta.data$orig.ident)
#-------------CD45+ CD19-  ---------------------------------------------------------
data$CD19.groups <- 'CD19.pos'
data$CD19.groups[WhichCells(data, expression= CD19 < 0.1)] <- 'CD19.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD19.groups')
head(data@meta.data)
data <- subset(data, subset = CD19.groups != "CD19.pos")
gc()
table(data@meta.data$orig.ident)
#-------------CD45+ CD19- TRAC- -------------------------------------------
data$TRAC.groups <- 'TRAC.pos'
data$TRAC.groups[WhichCells(data, expression= TRAC < 0.1)] <- 'TRAC.neg'
DimPlot(data, reduction = 'pca',split.by = 'TRAC.groups')
head(data@meta.data)
data <- subset(data, subset = TRAC.groups != "TRAC.pos")
gc()
table(data@meta.data$orig.ident)
#-------------CD45+ CD19- TRAC- CD14+     -------------------------------------------
data$CD14.groups <- 'CD14.pos'
data$CD14.groups[WhichCells(data, expression= CD14 < 0.1)] <- 'CD14.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD14.groups')
head(data@meta.data)
data <- subset(data, subset = CD14.groups != "CD14.neg")
gc()
table(data@meta.data$orig.ident)
#----------Isolate MILR1+  -------------------------------------
data$MILR1.groups <- 'MILR1.pos'
data$MILR1.groups[WhichCells(data, expression= MILR1 < 0.1)] <- 'MILR1.neg'
DimPlot(data, reduction = 'pca',split.by = 'MILR1.groups')
head(data@meta.data)
data <- subset(data, subset = MILR1.groups != "MILR1.neg")
gc()
table(data@meta.data$orig.ident)
VlnPlot(data, features = c('MILR1'),cols = c())
#---------------------------------------------------------------------------
DimPlot(data,dims = c(1,2),reduction = 'pca',cols = c(),pt.size = 0.5)
table(data@meta.data$orig.ident)
#----------------------------------------------------
setwd()
list<-read.csv('gene_list.csv')
gene_list<-list$x
#----split data
reg<-FetchData(data,vars = c('ident',gene_list),slot = 'counts')
table(reg$ident)
reg$ident<-ifelse(reg$ident=='HCC', 1, 0)
table(reg$ident)
#------Even out group numbers and shuffle
edit<-reg[-c(76:246),]
table(edit$ident)
reg<-edit[sample(1:nrow(edit)),]
table(reg$ident)
test<-reg
#--------------------------------------------------
set.seed(14)
#----------run model
setwd()
model<-read_rds('hcc5_model.rda')
vif(model)
summary(model)
varImp(model)
logLik(model)
#------Displaying variance inflation factors
vif(model)
#------Displaying variable importances
varImp(model)
newdata = test
summary(newdata)
dim(newdata)
dim(test)
summary(predict(model, newdata, type = 'response'))
predicted<-predict(model, newdata, type = 'response')
pred_factor <- predicted
pred_factor<- round(pred_factor)
pred_factor<-as.factor(pred_factor)
summary(pred_factor)
table(test$ident)
actual<-as.factor(test$ident)
caret::confusionMatrix(pred_factor, actual,positive='1')
table(actual)
table(pred_factor)
actuals<-test$ident
df<-cbind(actuals,predicted)
df<-data.frame(df,check.names = FALSE)
optCutoff<-optimalCutoff(actuals = actuals,
                         predictedScores = predicted,
                         optimiseFor = "Ones",
                         returnDiagnostics = TRUE)
head(optCutoff)
auc(actuals,predicted)
str(df$predicted)
ROC_pred<-prediction(df$predicted,df$actuals)
ROC_perf<-performance(ROC_pred,'tpr','fpr')
plot(ROC_perf,colorize=TRUE,print.cutoffs.at=seq(0.1,by=0.1))
summary(model)
logLik(model)
caret::confusionMatrix(pred_factor, actual,positive='1')
table(data@meta.data$orig.ident)
table(test$ident)
break
