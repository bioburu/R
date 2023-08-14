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
setwd('/home/amp_prog/Desktop/hope_model/GSE237779_immune.cells_GBM')
data<-read_rds('GSE237779_integrated_seurat_object_CD45_scRNAseq.rds')
#-----Seurat 
VlnPlot(data, features = c('MILR1','PSAP','CD68','MRC1','ITGAM'),pt.size=0.1)
FindMarkers(data, ident.1 = 'Mono/Macro', ident.2 = 'B', features = c('MILR1','PSAP','CD68','MRC1','ITGAM'))
#-----regression modeling
data<-FetchData(data,vars = c('ident','MILR1','PSAP','CD68','MRC1','ITGAM'),slot = 'counts')
data$ident<-ifelse(data$ident=='Mono/Macro', 1, 0)
set.seed(3)
subsample<-caTools::sample.split(data, SplitRatio=0.7)
train<-subset(data, subsample==TRUE)
train<-na.omit(train)
dim(train)
test<-subset(data,subsample==FALSE)
test<-na.omit(test)
dim(test)
#----------run model
model<-glm(ident~MILR1+PSAP+CD68+MRC1+ITGAM,data = train, family = binomial)
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
predicted<- predict(model, data.frame(MILR1=test$MILR1,PSAP=test$PSAP,CD68=test$CD68,MRC1=test$MRC1,ITGAM=test$ITGAM), type='response')
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
