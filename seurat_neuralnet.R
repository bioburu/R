library(ISLR2)
library(keras3)
library(ggplot2)
library(InformationValue)
library(neuralnet)
library(pROC)
library(ROCR)
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggridges)
library(caTools)
library(car)
library(caret)
setwd('/home/em_b/work_stuff/human_develop_timecourse')
#-------extract matrix from seurat object 
data<-readRDS('/home/em_b/work_stuff/human_develop_timecourse/GSE217511_RAW/save_file2.rda')
gene_list<-row.names(data@assays$RNA$counts)
table(data@meta.data$cell_types)
#------subset group types
#data<-subset(data,subset=cell_types=='opc')
#data@meta.data
#table(data@meta.data$cell_types)
#head(data@meta.data)
#-------extract data
Data<-subset(data,orig.ident=='17wks')
Data2<-subset(data,orig.ident=='20wks')
data<-merge(Data,Data2)
table(data@meta.data$orig.ident)
head(data@meta.data)
matrix<-FetchData(data,vars = c('orig.ident',gene_list),layer = 'counts')
table(matrix$orig.ident)
#matrix<-subset(matrix,
#       orig.ident%in%c('17wks','20wks'))
#Idents(data)<-data@meta.data$orig.ident
#------find genes
Idents(data)<-data@meta.data$orig.ident
head(Idents(data))
data<-JoinLayers(data)
x<-FindMarkers(data,
            ident.1 = '20wks',
            ident.2 = '17wks',
            logfc.threshold = 3,
            only.pos=TRUE,
            test.use = 'bimod')
VlnPlot(data,
        features = c(row.names(x))[1:12],
        layer = 'data')
#----model candidates
cat(row.names(x)[1:40],sep='+')
matrix$orig.ident<-ifelse(matrix$orig.ident=='20wks', 1, 0)
table(matrix$orig.ident)
#------Even out group numbers and shuffle
matrix<-sample(matrix)
table(matrix$orig.ident)
sample <- sample(c(TRUE, FALSE), nrow(matrix), replace=TRUE, prob=c(0.5,0.5))
train  <- matrix[sample, ]
table(train$orig.ident)
test   <- matrix[!sample, ]
table(test$orig.ident)
set.seed(13)
str(matrix)
#------------logistic regression
model<-glm(orig.ident~PLCG2+PDZD2+SYN3+GAS7+CELF4+
             CACNA1A+LINC00854+RYR3+SGSM1+ASIC2+
             RYR2+GRIK4+GABRG3+DLC1+NGEF+
             PCP4+SPTBN4+MTUS2+GSG1L+MICAL2+
             RAP1GAP2+PLEKHH2,
           data = train, family = binomial)
vif(model)
summary(model)
varImp(model)
logLik(model)
summary(predict(model, test, type = 'response'))
predicted<-predict(model, test, type = 'response')
pred_factor <- predicted
pred_factor<- round(pred_factor)
pred_factor<-as.factor(pred_factor)
summary(pred_factor)
table(test$orig.ident)
actual<-as.factor(test$orig.ident)
confusion_matrix<-caret::confusionMatrix(pred_factor, actual,positive='1')
actuals<-test$orig.ident
df<-cbind(actuals,predicted)
df<-data.frame(df,check.names = FALSE)
optCutoff<-optimalCutoff(actuals = actuals,
                         predictedScores = predicted,
                         optimiseFor = "Ones",
                         returnDiagnostics = TRUE)
head(optCutoff)
auc(actuals,predicted)
ROC_pred<-prediction(df$predicted,df$actuals)
ROC_perf<-performance(ROC_pred,'tpr','fpr')
confusion_matrix<-caret::confusionMatrix(pred_factor, actual,positive='1')
confusion_matrix
fourfoldplot(as.table(confusion_matrix),
             color = c('grey','skyblue'),
             std='all.max',
             main='10wks=0 20wks=1')
plot.roc(actuals, predicted,
         percent = TRUE,
         main = 'AUC: 17wks>20wks',
         add =  FALSE,
         asp = NA,
         print.auc = TRUE,
         print.auc.col='red',
         print.auc.cex=1,
         grid=TRUE,
         grid.col='skyblue',
         identity.col='blue',
         identity.lty=8,
         col='red',
         print.thres=TRUE,
         print.thres.pch=5,
         print.thres.col='black',
         ci=TRUE)
list<-c('PLCG2','PDZD2','SYN3','GAS7','CELF4',
        'CACNA1A','LINC00854','RYR3','SGSM1','ASIC2',
        'RYR2','GRIK4','GABRG3','DLC1','NGEF',
        'PCP4','SPTBN4','MTUS2','GSG1L','MICAL2',
        'RAP1GAP2','PLEKHH2')
DoHeatmap(
  data,
  features = list,
  cells = NULL,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  vjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)
#-------------neural network
nn=neuralnet(orig.ident~PLCG2+PDZD2+SYN3+GAS7+CELF4+
               CACNA1A+LINC00854+RYR3+SGSM1+ASIC2+
               RYR2+GRIK4+GABRG3+DLC1+NGEF+
               PCP4+SPTBN4+MTUS2+GSG1L+MICAL2+
               RAP1GAP2+PLEKHH2,data=train, 
             hidden=c(20,12),act.fct = "logistic",
             linear.output = FALSE)
plot(nn)
dev.print(pdf, 'nn_plot.pdf') 
#-------------------------------------------------------------------------------
summary(predict(nn, test, type = 'response'))
predicted2<-predict(nn, test, type = 'response')
pred_factor2 <- predicted2
pred_factor2<- round(pred_factor2)
pred_factor2<-as.factor(pred_factor2)
summary(pred_factor2)
table(test$orig.ident)
actual2<-as.factor(test$orig.ident)
confusion_matrix2<-caret::confusionMatrix(pred_factor2, actual2,positive='1')
confusion_matrix2
actuals2<-test$orig.ident
df2<-cbind(actuals2,predicted2)
df2<-data.frame(df2,check.names = FALSE)
optCutoff2<-optimalCutoff(actuals = actuals2,
                         predictedScores = predicted2,
                         optimiseFor = "Ones",
                         returnDiagnostics = TRUE)
head(optCutoff2)
auc(actuals2,predicted2)
#str(df$predicted)
ROC_pred<-prediction(df2$V2,df2$actuals2)
ROC_perf<-performance(ROC_pred,'tpr','fpr')
plot(ROC_perf,colorize=TRUE,print.cutoffs.at=seq(0.1,by=0.1))
summary(nn)
caret::confusionMatrix(pred_factor2, actual2,positive='1')
fourfoldplot(as.table(confusion_matrix2),
             color = c('grey','skyblue'),
             std='all.max',
             main='17wks=0 20wks=1')
plot.roc(actuals2, predicted2,
         percent = TRUE,
         main = 'AUC:_17wks>20wks_neural_network',
         add =  FALSE,
         asp = NA,
         print.auc = TRUE,
         print.auc.col='red',
         print.auc.cex=1,
         grid=TRUE,
         grid.col='skyblue',
         identity.col='blue',
         identity.lty=8,
         col='red',
         print.thres=TRUE,
         print.thres.pch=5,
         print.thres.col='black',
         ci=TRUE)
break
