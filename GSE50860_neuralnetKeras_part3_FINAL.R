#https://hastie.su.domains/ISLR2/Labs/Rmarkdown_Notebooks/Ch10-deeplearning-lab-keras.html
library(ISLR2)
library(keras)
library(ggplot2)
library(InformationValue)
library(neuralnet)
library(pROC)
library(ROCR)
print("Uploading file and editing")
df<-read.csv(file = file.choose())
row.names(df)<-df$X
colnames(df)[2]<-'Genotype'
df$Genotype<-ifelse(df$Genotype=='mutated',1,0)
df<-df[,-1]
print("Splitting into train and test sets")
print("Adjust split ratio here as needed:")
subsample<-caTools::sample.split(df, SplitRatio=0.3)
train<-subset(df, subsample==TRUE)
dim(train)
train<-na.omit(train)
dim(train)
print("Assigning testing data set")
test<-subset(df,subsample==FALSE)
dim(test)
test<-na.omit(test)
dim(test)
str(df)
print("Fitting with Neuralnet")
nn=neuralnet(Genotype~
              CD34+ATP10A+TRH+POU4F1+CTAG2+LAMP5+SERPINF1+CSRP2+DNTT+MAN1A1+
              HOXA4+MREG+MEIS1+HOXA6+HOXA7+PKLR+PBX3+CTSE+ROBO1+ITM2C+
              HPGDS+DPP4+HOXA5+SPON1+BAALC+GNG7+NPR3+PKLR+ARHGAP22+HOXA1+
              PPP1R26+APP+CD200+SIDT1+KIF17+ASB9+S100B+POLE+GYPC+FBLN5+
              SPON1+HOMER2+FZD6+CDKN2C+SPIB+FECH+GPA33+GALC+GJA1+ITM2A
            , data=train, 
            hidden=c(20,12),act.fct = "logistic",
            linear.output = FALSE)
plot(nn)
#-------------------------------------------------------------------------
print("Prepping test data for keras tensorflow")
id<-nrow(test)
total_rows <- nrow(df)
testid <- sample(1:total_rows, id)
testid
print("Convert data frame to x and y tensors for Keras")
print("x matrix is without intercept (-1) and scaled with mean of 0 with std of 1")
x <- scale(model.matrix(Genotype ~ . - 1, data = df))
print("y matrix is response integer")
y <- df$Genotype
print("Make model with 2 hidden layers with same hidden units as Neuralnet and RedLU activation function")
print("Drop out layer rate of 40% per iteration")
print("output layer at 1 for binary response")
keras <- keras_model_sequential() %>%
  layer_dense(input_shape = ncol(x),
              units = 20, activation = "relu") %>%
  layer_dense(units = 12, activation = "softmax") %>%
  layer_dropout(rate = 0.6)%>%
  layer_dense(units=1)
#
print("Compile function to output mean absolute error")
keras %>% compile(loss = "mse", optimizer = optimizer_rmsprop(),metrics = list("mean_absolute_error"))
#keras %>% compile(loss = 'categorical_crossentropy', optimizer = 'adam', metrics = c('accuracy'))

print("Enter data frame into model now")
history <- keras %>% fit(x[-testid, ], y[-testid], epochs = 400, batch_size = 50,
                      validation_data = list(x[testid, ], y[testid])
)
keras %>%evaluate(x[testid, ], y[testid])
#
print("Confusion matrix for Neural Net model")
predict<-compute(nn,test)
predict<- round(predict$net.result)
predict<-as.factor(predict)
act<-as.factor(test$Genotype)
detach(package:neuralnet,unload = T)
caret::confusionMatrix(act, predict)
print("Standard error rate")
head(nn$result.matrix)
table(df$Genotype)
predicted.scores<-predict
actual.scores<-test$Genotype
df2<-cbind(actual.scores,predicted.scores)
df2<-data.frame(df2,check.names = FALSE)
print("Displaying AUC score of Neural Net model")
auc(df2$actual.scores,df2$predicted.scores)
print("Plotting ROC curve with Neural Net model")
ROC_pred<-prediction(df2$predicted.scores,df2$actual.scores)
ROC_perf<-performance(ROC_pred,'tpr','fpr')
#plot(ROC_perf,colorize=TRUE,print.cutoffs.at=seq(0.1,by=0.1))
print("Confusion matrix for keras model")
predicted <- predict(keras, x[testid, ])
predicted[predicted < 0] <- 0     
predicted[predicted >1] <- 1 
predicted<-round(predicted)
predicted<-as.factor(predicted)
predicted
testid
actuals<-df$Genotype[testid]
actuals<-as.factor(actuals)
caret::confusionMatrix(actuals, predicted)
print("Compare tensorflow with neuralnet")
print("Tensorflow error:")
keras %>%evaluate(x[testid, ], y[testid])
print("Neuralnet error:")
head(nn$result.matrix)
print("Did Keras beat me?")
plot(history)
#plot(nn)
break









