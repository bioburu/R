library(caTools)
library(car)
library(caret)
library(InformationValue)
library(pROC)
library(ROCR)
library(Seurat)
#-------------Set working directory
setwd('/home/amp_prog/Desktop/hope_model/GSE208653_hpv_cancer')
#-------------Add files for cancer group                
features_path <- 'GSM6360686_SCC_4.features.tsv.gz'
barcodes_path <- 'GSM6360686_SCC_4.barcodes.tsv.gz'
matrix_path <- 'GSM6360686_SCC_4.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
matrix<-as.data.frame(matrix[c('CD274','EBNA1BP2'),])
factor<-rep(c(1),times=ncol(matrix))
df<-rbind(factor,matrix)
rownames(df)[1]<-'disease_status'
#-------------Add files for lesions group
features_path <- 'GSM6360684_HSIL_1.features.tsv.gz'
barcodes_path <- 'GSM6360684_HSIL_1.barcodes.tsv.gz'
matrix_path <- 'GSM6360684_HSIL_1.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
matrix<-as.data.frame(matrix[c('CD274','EBNA1BP2'),])
factor<-rep(c(0),times=ncol(matrix))
matrix<-rbind(factor,matrix)
df<-cbind(df,matrix)
#--------------Add files for hpv group
features_path <- 'GSM6360682_N_1.features.tsv.gz'
barcodes_path <- 'GSM6360682_N_1.barcodes.tsv.gz'
matrix_path <- 'GSM6360682_N_1.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
matrix<-as.data.frame(matrix[c('CD274','EBNA1BP2'),])
factor<-rep(c(0),times=ncol(matrix))
matrix<-rbind(factor,matrix)
df<-cbind(df,matrix)
#-------------no hpv
features_path <- 'GSM6360680_N_HPV_NEG_1.features.tsv.gz'
barcodes_path <- 'GSM6360680_N_HPV_NEG_1.barcodes.tsv.gz'
matrix_path <- 'GSM6360680_N_HPV_NEG_1.matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
matrix<-as.data.frame(matrix[c('CD274','EBNA1BP2'),])
factor<-rep(c(0),times=ncol(matrix))
matrix<-rbind(factor,matrix)
df<-cbind(df,matrix)
#-----------Subsample groups
set.seed(3)
df<-as.data.frame(t(df))
summary(df)
str(df)
subsample<-caTools::sample.split(df, SplitRatio=0.7)
train<-subset(df, subsample==TRUE)
dim(train)
train<-na.omit(train)
dim(train)
test<-subset(df,subsample==FALSE)
dim(test)
test<-na.omit(test)
dim(test)
break 
#----------run model
model<-glm(disease_status~TNFAIP3,data = train, family = binomial)
summary(model)
model<-glm(disease_status~NR1H3,data = train, family = binomial)
summary(model)
model<-glm(disease_status~CD14,data = train, family = binomial)
summary(model)
model<-glm(disease_status~EBNA1BP2,data = train, family = binomial)
summary(model)
model<-glm(disease_status~CD274+EBNA1BP2,data = train, family = binomial)
summary(model)
logLik(model)
table(df$disease_status)
print("Calculating Mcfadden's pseudo R squared")
null<-model$null.deviance/-2
resdDEV<-model$deviance/-2
pR2<-(null-resdDEV)/null
print(pR2)
print("Displaying variance inflation factors")
vif(model)
print("Displaying variable importance factors")
varImp(model)
dim(test)
dim(train)
predict<- predict(model, data.frame(EBNA1BP2=test$EBNA1BP2,CD274=test$CD274), type="response")
print(head(predict))
length(predict)
pred_factor <- predict
print(head(pred_factor))
pred_factor<- round(pred_factor)
print(head(pred_factor))
pred_factor<-as.factor(pred_factor)
print(head(pred_factor))
actual<-as.factor(test$disease_status)
caret::confusionMatrix(actual, pred_factor)
predicted.scores<-predict
head(predicted.scores)
actual.scores<-test$disease_status
head(actual.scores)
df2<-cbind(actual.scores,predicted.scores)
df2<-data.frame(df2,check.names = FALSE)
str(df2)
optCutoff<-optimalCutoff(actuals = actual.scores,
                         predictedScores = predicted.scores,
                         optimiseFor = "Ones",
                         returnDiagnostics = TRUE)
optCutoff
auc(actual.scores,predicted.scores)
ROC_pred<-prediction(df2$predicted.scores,df2$actual.scores)
ROC_perf<-performance(ROC_pred,'tpr','fpr')
plot(ROC_perf,colorize=TRUE,print.cutoffs.at=seq(0.1,by=0.1))
table(df$disease_status)
summary(model)
logLik(model)
caret::confusionMatrix(actual, pred_factor)

