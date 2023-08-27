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
setwd()
data<-read_rds('GSE237779_integrated_seurat_object_CD45_scRNAseq.rds')
#-----Seurat 
VlnPlot(data, features = c(),pt.size=0.1)
FindMarkers(data, ident.1 = '', ident.2 = '', features = c())
#-----regression modeling
table(data@meta.data$orig.ident)
reg<-FetchData(data,vars = c(),slot = 'counts')
reg$ident<-ifelse(reg$ident=='', 1, 0)
set.seed(3)
sample <- sample(c(TRUE, FALSE), nrow(reg), replace=TRUE, prob=c(0.7,0.3))
train  <- reg[sample, ]
test   <- reg[!sample, ]
#----------run model
model<-glm(ident~,data = train, family = binomial)
summary(model)
logLik(model)
table(data$ident)
#----------Mcfadden's pseudo R squared
null<-model$null.deviance/-2
resdDEV<-model$deviance/-2
pR2<-(null-resdDEV)/null
print(pR2)
#------Displaying variance inflation factors
vif(model)
#------Displaying variable importance factors
varImp(model)
dim(test)
dim(train)
newdata = data.frame(MRC1=30, CD14=0)
newdata
summary(predict(model, newdata, type = 'response'))
break 
predicted<- predict(model,newdata , type='response')
print(head(predicted))
length(predicted)
pred_factor <- predicted
pred_factor<- round(pred_factor)
pred_factor<-as.factor(pred_factor)
print(head(pred_factor))
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
