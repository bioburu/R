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
setwd('/home/amp_prog/Desktop/hope_model/GSE208653_hpv_cancer')
features_path <- 'GSM6360680_N_HPV_NEG_1.features.tsv.gz'
barcodes_path <- 'GSM6360680_N_HPV_NEG_1.barcodes.tsv.gz'
matrix_path <- 'GSM6360680_N_HPV_NEG_1.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
no_hpv <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'no_hpv')
#---------------------------------------------------------------------------
features_path <- 'GSM6360682_N_1.features.tsv.gz'
barcodes_path <- 'GSM6360682_N_1.barcodes.tsv.gz'
matrix_path <- 'GSM6360682_N_1.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
hpv <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'hpv')
#---------------------------------------------------------------------------
features_path <- 'GSM6360684_HSIL_1.features.tsv.gz'
barcodes_path <- 'GSM6360684_HSIL_1.barcodes.tsv.gz'
matrix_path <- 'GSM6360684_HSIL_1.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
lesions <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'lesions')
#---------------------------------------------------------------------------
features_path <- 'GSM6360686_SCC_4.features.tsv.gz'
barcodes_path <- 'GSM6360686_SCC_4.barcodes.tsv.gz'
matrix_path <- 'GSM6360686_SCC_4.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
cervical_cancer <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'cervical_cancer')
data<-merge(no_hpv,y=c(hpv,lesions,cervical_cancer),project='hpv_cervical_cancer')
head(data@active.ident)
rm(no_hpv,hpv,lesions,cervical_cancer,matrix)
gc()
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt <10)
data <- NormalizeData(data)
gc()
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top10 <- head(VariableFeatures(data), 10)
top10
gc()
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
gc()
gc()
gc()
data <- ScaleData(data, features = all.genes)
gc()
dim(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
gc()
table(data@meta.data$orig.ident)
#----------Isolate CD45+  -------------------------------------
data$CD45.groups <- 'CD45.pos'
data$CD45.groups[WhichCells(data, expression= PTPRC < 0.1)] <- 'CD45.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD45.groups')
head(data@meta.data)
data <- subset(data, subset = CD45.groups != "CD45.neg")
gc()
#-------------CD45+ CD19-  ---------------------------------------------------------
data$CD19.groups <- 'CD19.pos'
data$CD19.groups[WhichCells(data, expression= CD19 < 0.1)] <- 'CD19.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD19.groups')
head(data@meta.data)
data <- subset(data, subset = CD19.groups != "CD19.pos")
gc()
#-------------CD45+ CD19- TRAC- -------------------------------------------
data$TRAC.groups <- 'TRAC.pos'
data$TRAC.groups[WhichCells(data, expression= TRAC < 0.1)] <- 'TRAC.neg'
DimPlot(data, reduction = 'pca',split.by = 'TRAC.groups')
head(data@meta.data)
data <- subset(data, subset = TRAC.groups != "TRAC.pos")
gc()
#-------------CD45+ CD19- TRAC- CD14+     -------------------------------------------
data$CD14.groups <- 'CD14.pos'
data$CD14.groups[WhichCells(data, expression= CD14 < 0.1)] <- 'CD14.neg'
DimPlot(data, reduction = 'pca',split.by = 'CD14.groups')
head(data@meta.data)
data <- subset(data, subset = CD14.groups != "CD14.neg")
gc()
#----------Isolate MILR1+  -------------------------------------
data$MILR1.groups <- 'MILR1.pos'
data$MILR1.groups[WhichCells(data, expression= MILR1 < 0.1)] <- 'MILR1.neg'
DimPlot(data, reduction = 'pca',split.by = 'MILR1.groups')
head(data@meta.data)
data <- subset(data, subset = MILR1.groups != "MILR1.neg")
gc()
VlnPlot(data, features = c('MRC1','CD14','DDX58','IFIH1','NFKB1','NFKB2','CSF1','MYC','MILR1'),cols = c())
FindMarkers(data, ident.1 = 'cervical_cancer', ident.2 = 'no_hpv', features = c('MRC1','CD14','DDX58','IFIH1','NFKB1','NFKB2','CSF1','MYC','MILR1'))
#-----regression modeling
reg<-FetchData(data,vars = c('ident','MRC1','CD14'),slot = 'counts')
reg$ident<-ifelse(reg$ident=='cervical_cancer', 1, 0)
set.seed(3)
sample <- sample(c(TRUE, FALSE), nrow(reg), replace=TRUE, prob=c(0.7,0.3))
train  <- reg[sample, ]
test   <- reg[!sample, ]
#----------run model
model<-glm(ident~MRC1+CD14,data = train, family = binomial)
summary(model)
logLik(model)
#----------Mcfadden's pseudo R squared
null<-model$null.deviance/-2
resdDEV<-model$deviance/-2
pR2<-(null-resdDEV)/null
print(pR2)
#------Displaying variance inflation factors
#vif(model)
#------Displaying variable importance factors
varImp(model)
newdata = data.frame(MRC1=test$MRC1,CD14=test$CD14)
summary(predict(model, newdata, type = 'response'))
predicted<-predict(model, newdata, type = 'response')
pred_factor <- predicted
pred_factor<- round(pred_factor)
pred_factor<-as.factor(pred_factor)
actual<-as.factor(test$ident)
caret::confusionMatrix(actual, pred_factor)
actuals<-test$ident
df<-cbind(actuals,predicted)
df<-data.frame(df,check.names = FALSE)
optCutoff<-optimalCutoff(actuals = actuals,
                         predictedScores = predicted,
                         optimiseFor = "Ones",
                         returnDiagnostics = TRUE)
head(optCutoff)
auc(actuals,predicted)
ROC_pred<-prediction(df$predicted,df$actual)
ROC_perf<-performance(ROC_pred,'tpr','fpr')
plot(ROC_perf,colorize=TRUE,print.cutoffs.at=seq(0.1,by=0.1))
summary(model)
logLik(model)
caret::confusionMatrix(actual, pred_factor)
break 


