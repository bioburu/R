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
library(scTenifoldNet)
setwd('/home/amp_prog/Desktop/TAM_manuscript/datasets/GSE228512_GBM.sera')
matrix<-read.csv('GSE228512_hiseq_counts.csv')
row.names(matrix)<-make.names(matrix$Gene,unique = TRUE)
matrix<-matrix[,-1]
data <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'normal')
gc()
data <- RenameIdents(object = data, `preop` = "GBM.EV")
table(data@active.ident)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
data <- NormalizeData(data)
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
data <- ScaleData(data, features = all.genes)
gc()
dim(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
DimHeatmap(data, dims = 1:15, cells = 500, balanced = T)
ElbowPlot(data)
gc()
table(data@meta.data$orig.ident)
#----split data
reg<-FetchData(data,vars = c('ident','GTF2IP1','RNA5SP30','SNORD89','GOLGA6L23P','MAGED4B'),slot = 'counts')
table(reg$ident)
reg$ident<-ifelse(reg$ident=='preop', 1, 0)
table(reg$ident)
set.seed(14)
sample <- sample(c(TRUE, FALSE), nrow(reg), replace=TRUE, prob=c(0.5,0.5))
train  <- reg[sample, ]
table(train$ident)
test   <- reg[!sample, ]
table(test$ident)
#----------run model
model<-glm(ident~GTF2IP1+RNA5SP30+SNORD89+GOLGA6L23P+MAGED4B,
           data = train, family = binomial)
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
VlnPlot(data, features = c('GTF2IP1','RNA5SP30','SNORD89','GOLGA6L23P','MAGED4B'),cols = c('grey','red'))
FindMarkers(data, ident.1 = 'preop', ident.2 = 'control', features = c('GTF2IP1','RNA5SP30','SNORD89','GOLGA6L23P','MAGED4B'))
VlnPlot(data, features = c('MILR1','FLT3','MRC1','IFIH1','CCL8','TAP1','NFKB1','CD274','PDCD1LG2','IL4','IL13','IFNA1','IFNB1','IFNG','TNF','IL12B'),cols = c('grey','red'))
FindMarkers(data, ident.1 = 'preop', ident.2 = 'control', features = c('MILR1','MRC1','DHX58','IFIH1','CCL8','TAP1','NFKB1','CD274','PDCD1LG2','IL4','IL13','IFNA1','IFNB1','IFNG','TNF','IL12B'))

