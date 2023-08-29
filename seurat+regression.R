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
setwd('/home/amp_prog/Desktop/TAM_manuscript/datasets/GSE221575_CRC')
#---------------------------------------------------------------------------
features_path <- 'GSM6886539_Rectum-03T_genes.tsv.gz'
barcodes_path <- 'GSM6886539_Rectum-03T_barcodes.tsv.gz'
matrix_path <- 'GSM6886539_Rectum-03T_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
x <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'tumor')
summary(x@active.ident)
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
features_path <- 'GSM6886540_Liver_metastasis_genes.tsv.gz'
barcodes_path <- 'GSM6886540_Liver_metastasis_barcodes.tsv.gz'
matrix_path <- 'GSM6886540_Liver_metastasis_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
y <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'meta_tumor')
summary(y@active.ident)
data<-merge(x,y=c(y),project='')
table(data@meta.data$orig.ident)
head(data@active.ident)
rm(x,y,z)
gc()
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt <25)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
gc()
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
VariableFeatures(data)
top10 <- head(VariableFeatures(data), 10)
top10
plot1 <- VariableFeaturePlot(data)
plot1
all.genes <- rownames(data)
all.genes
data <- ScaleData(data, features = all.genes)
gc()
dim(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
gc()
table(data@meta.data$orig.ident)
DimPlot(data,dims = c(1,2),reduction = 'pca',cols = c(),pt.size = 4)
VlnPlot(data, features = c('PTPRC','CD19','TRAC','CD14','MILR1','MKI67','FCGR3A','MRC1'),cols = c())
FindMarkers(data, ident.1 = 'meta_tumor', ident.2 = 'tumor', features = c('PTPRC','TRAC','CD14','MILR1','MKI67','FCGR3A','MRC1'))
#-----regression modeling
reg<-FetchData(data,vars = c('ident','CD14','FCGR3A','MILR1','MRC1','PTPRC'),slot = 'counts')
table(reg$ident)
reg$ident<-ifelse(reg$ident=='meta_tumor', 1, 0)
table(reg$ident)
set.seed(387)
sample <- sample(c(TRUE, FALSE), nrow(reg), replace=TRUE, prob=c(0.9,0.1))
train  <- reg[sample, ]
table(train$ident)
test   <- reg[!sample, ]
table(test$ident)
#----------run model
model<-glm(ident~MRC1+FCGR3A+CD14+MILR1+PTPRC,data = train, family = binomial)
summary(model)
logLik(model)
#----------Mcfadden's pseudo R squared
null<-model$null.deviance/-2
resdDEV<-model$deviance/-2
pR2<-(null-resdDEV)/null
print(pR2)
#------Displaying variance inflation factors
vif(model)
#------Displaying variable importance factors
varImp(model)
newdata = data.frame(PTPRC=test$PTPRC,FCGR3A=test$FCGR3A,MILR1=test$MILR1,CD14=test$CD14,MRC1=test$MRC1)
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
ROC_pred<-prediction(df$predicted,df$actual)
ROC_perf<-performance(ROC_pred,'tpr','fpr')
plot(ROC_perf,colorize=TRUE,print.cutoffs.at=seq(0.1,by=0.1))
summary(model)
logLik(model)
caret::confusionMatrix(pred_factor, actual,positive='1')
table(data@meta.data$orig.ident)
table(test$ident)
break 
VlnPlot(data, features = c('PTPRC','MRC1','MILR1','CD14','FCGR3A'),cols = c('red','grey','grey'))

