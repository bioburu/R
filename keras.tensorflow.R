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
library(plotly)
library(Signac)
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(htmlwidgets)
require('Seurat.utils')
library(plotly)
library(EnsDb.Mmusculus.v79)
library(celldex)
library(SingleR)
library(biomaRt)
library(ReactomePA)
library(clusterProfiler)
setwd('/home/deviancedev01/Desktop/interesting_exps/GSE224679_RAW')
#-------------------------------------------------------------------------------
features_path <- 'GSM7029393_111_11_features.tsv.gz'
barcodes_path <- 'GSM7029393_111_11_barcodes.tsv.gz'
matrix_path <- 'GSM7029393_111_11_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
lesions <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'lesions')
summary(lesions@active.ident)
#-----------------------------
features_path <- 'GSM7029395_111_13_features.tsv.gz'
barcodes_path <- 'GSM7029395_111_13_barcodes.tsv.gz'
matrix_path <- 'GSM7029395_111_13_matrix.mtx.gz'
matrix <- ReadMtx(mtx= matrix_path, features = features_path, cells= barcodes_path)
normal <- CreateSeuratObject(counts=matrix,min.cells=20,min.features=200,project = 'normal')
summary(normal@active.ident)
lesions<-subset(lesions,
                downsample=2760)
data<-merge(normal,
           lesions)
table(Idents(data))
data<-JoinLayers(data)
gene_list<-row.names(data@assays$RNA)
table(data@meta.data$orig.ident)
#Data<-subset(data,orig.ident=='17wks')
#Data2<-subset(data,orig.ident=='20wks')
#data<-merge(Data,Data2)
matrix<-FetchData(data,vars = c('orig.ident',gene_list),layer = 'counts')
#Idents(data)<-data@meta.data$orig.ident
#data<-JoinLayers(data)
matrix$orig.ident<-ifelse(matrix$orig.ident=='lesions', 1, 0)
#------Even out group numbers and shuffle
set.seed(1435)
matrix<-matrix[sample(1:nrow(matrix)), ]
sample <- sample(c(TRUE, FALSE), nrow(matrix), replace=TRUE, prob=c(0.5,0.5))
train  <- matrix[sample, ]
test   <- matrix[!sample, ]
#-------------------------------------------------------------------------
print("Prepping test data for keras tensorflow")
#------use the split groups to subsample from entire dataset
number.cells.train<-nrow(train)
number.cells.all <- nrow(matrix)
#------extract test number of observations and sample
testid <- sample(1:number.cells.all,
                 number.cells.train)
testid
#-------------------------
print("Convert data frame to x and y tensors for Keras")
print("x matrix is without intercept (-1) and scaled with mean of 0 with std of 1")
#---make a matrix where the reference control is the orig.ident
#----build model contrast with independent variables
x <- scale(model.matrix(orig.ident ~ . - 1,
                        data = matrix))
print("y matrix is response integer")
#---assign response variable
y <- matrix$orig.ident
y
print("Make model with 2 hidden layers with same hidden units as Neuralnet and RedLU activation function")
print("Drop out layer rate of 40% per iteration")
print("output layer at 1 for binomial response")
#--sudo apt-get install python3-venv
library(tensorflow)
#install_tensorflow(envname = "r-tensorflow")
keras <- keras_model_sequential() %>%
  layer_dense(input_shape = ncol(x),
              units = 20, activation = "relu") %>%
  layer_dense(units = 12, activation = "softmax") %>%
  layer_dropout(rate = 0.6)%>%
  layer_dense(units=1)
keras
#
keras %>% compile(
  optimizer = 'adam',
  loss = 'binary_crossentropy',
  metrics = c('accuracy')
)
#keras %>% compile(loss = 'categorical_crossentropy', optimizer = 'adam', metrics = c('accuracy'))
print("Enter data frame into model now")
history <- keras %>% fit(x[-testid, ], y[-testid], epochs = 100, batch_size = 50,
                         validation_data = list(x[testid, ], y[testid])
)
keras %>%evaluate(x[testid, ], y[testid])
#
print("Confusion matrix for keras model")
predicted <- predict(keras, x[testid, ])
predicted
predicted[predicted < 0] <- 0     
predicted[predicted >1] <- 1 
predicted<-round(predicted)
predicted<-as.factor(predicted)
predicted
actuals<-matrix$orig.ident[testid]
actuals<-as.factor(actuals)
caret::confusionMatrix(actuals, predicted)
plot(history)
print(history)
class(keras)
break 
#model_type.keras.models.Sequential <- function(keras) {
#  "classification"}

#predict_model.keras.models.Sequential <- function (x, newdata, type, ...) {
#  pred <- predict(object = x, x = as.matrix(newdata))
#  data.frame (Positive = pred, Negative = 1 - pred) }

#predict_model.keras.models.Sequential <- function (keras, matrix, type, ...) {
#  pred <- predict(object = keras, x = as.matrix(matrix))
#  data.frame (Positive = pred, Negative = 1 - pred) }
#predict_model (x       = keras, 
#               newdata = matrix, 
#               type    = 'raw') %>%
#  tibble::as_tibble()

model<-as_classifier(keras, labels = NULL)
model_type(model)

explainer <- lime::lime (
  x              = matrix, 
  model          = model, 
  bin_continuous = FALSE)
explanation <- lime::explain (
  matrix[1:10, ], # Just to show first 10 cases
  explainer    = explainer, 
  n_labels     = 1, # explaining a `single class`(Polarity)
  n_features   = 4, # returns top four features critical to each case
  kernel_width = 0.5)

