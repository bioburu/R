library(keras3)
library(lime)
library(tidyquant)
library(rsample)
library(recipes)
library(yardstick)
library(corrr)
library(tidyverse)
library(tensorflow)
library(reticulate)
library(Seurat)
library(h2o)
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(org.Hs.eg.db)
library(SeuratWrappers)
library(plotly)
library(htmlwidgets)
library(caTools)
library(car)
library(caret)
library(InformationValue)
library(pROC)
library(ROCR)
library(celldex)
library(SingleR)
library(plotly)
library(data.tree)
library(DiagrammeR)
library(rpart)
library(rpart.plot)
library(h2o)
tf$config$list_physical_devices("GPU")
#use_condaenv('/home/deviancedev01/.virtualenvs/r-tensorflow')
tf$constant("Hello TensorFlow!")
tensorflow::set_random_seed(42)
setwd('/home/deviancedev01/Desktop/GSE145370_escc')
#-------------------------------------
validation<-readRDS('escc_hsc_validation.rds')
VlnPlot(validation,features=c('CD274','TP53'))
validation <- SetIdent(validation, value = "orig.ident")
table(Idents(validation))
validation<-subset(validation,downsample=102)
table(Idents(validation))
all.genes<-row.names(validation)
markers<-FindMarkers(validation,features=c(all.genes),
                     ident.1 = 'escc.tumor',
                     logfc.threshold=0.5,
                     min.pct=0.5,
                     only.pos=FALSE)
markers<-row.names(markers)
cat(markers)
DoHeatmap(validation,
          features=c(markers))+NoLegend()
df<-FetchData(validation,vars = c('ident',markers),layer = 'counts')
set.seed(1232)
levels(df$ident)
table(df$ident)
df$ident<-ifelse(df$ident=='escc.tumor', 1, 0)
table(df$ident)
df<-arrange(df,ident)
df<-df[sample(1:nrow(df)),]
table(df$ident)
test<-df

training<-readRDS('escc_hsc_training.rds')
VlnPlot(training,features=c('CD274','TP53'))
training <- SetIdent(training, value = "orig.ident")
table(Idents(training))
training<-subset(training,downsample=426)
table(Idents(training))
all.genes<-row.names(training)
#-------all genes
df<-FetchData(training,vars = c('ident',markers),layer = 'counts')
set.seed(1232)
levels(df$ident)
table(df$ident)
df$ident<-ifelse(df$ident=='escc.tumor', 1, 0)
table(df$ident)
df<-arrange(df,ident)
df<-df[sample(1:nrow(df)),]
table(df$ident)
train <-df

x_train <- train %>% dplyr::select(-ident)
str(x_train)

x_test  <- test %>% dplyr::select(-ident)
str(x_test)

y_train <- train %>% dplyr::select(ident)
y_test  <-test %>% dplyr::select(ident)

model_keras <- keras_model_sequential()
model_keras %>% layer_dense(
    units              = 16, 
    kernel_initializer = "uniform", 
    activation         = "leaky_relu", 
    input_shape        = ncol(x_train)) %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(
    units              = 16, 
    kernel_initializer = "uniform", 
    activation         = "sigmoid") %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(
    units              = 1, 
    kernel_initializer = "uniform", 
    activation         = "sigmoid") %>%
  compile(
    optimizer = 'adam',
    loss      = 'binary_crossentropy',
    metrics   = c('accuracy')
  )

history <- fit(
  object           = model_keras, 
  x                = as.matrix(x_train), 
  y                = y_train,
  batch_size       = 10, 
  epochs           = 50,
  validation_split = 0.50
)
print(history)
plot(history) 
# Predicted Class
yhat_keras_class_vec <- predict(object = model_keras, x = as.matrix(x_test)) %>%
  as.vector()

# Predicted Class Probability
yhat_keras_prob_vec  <- predict(object = model_keras, x = as.matrix(x_test)) %>%
  as.vector()


# Format test data and predictions for yardstick metrics
estimates_keras_tbl <- tibble(
  truth      = as.factor(y_test$ident) %>% fct_recode(yes = "1", no = "0"),
  estimate   = as.factor(ifelse(yhat_keras_class_vec>0.5,'yes','no')) %>% fct_recode(yes = "1", no = "0"),
  class_prob = yhat_keras_prob_vec
)
str(estimates_keras_tbl)
# Confusion Table
estimates_keras_tbl %>% conf_mat(truth, estimate)
confusion_matrix<-caret::confusionMatrix(estimates_keras_tbl$truth, estimates_keras_tbl$estimate,
                                         positive='yes')
fourfoldplot(as.table(confusion_matrix),
             color = c('grey','skyblue'),
             std='all.max',
             main='Keras tensorflow model predictions')
confusion_matrix

# Accuracy
estimates_keras_tbl %>% metrics(truth, estimate)

# Setup lime::model_type() function for keras
model_type.keras.src.models.sequential.Sequential <- function(x, ...) {
  "classification"
}


# Setup lime::predict_model() function for keras
predict_model.keras.src.models.sequential.Sequential <- function(x, newdata, type, ...) {
  pred <- predict(object = x, x = as.matrix(newdata))
  data.frame(Yes = pred, No = 1 - pred)
}

# Test our predict_model() function
predict_model(x = model_keras, newdata = x_test, type = 'raw') %>%
  tibble::as_tibble()

explainer <- lime(x_train, model_keras)

explanation <- lime::explain(x_test[1:6,], explainer, n_labels = 1, n_features = 6)

plot_features(explanation) +
  labs(title = "LIME Feature Importance Visualization",
       subtitle = "Hold Out (Test) Set, First 10 Cases Shown")

final<-unique(explanation$feature)
markers<-final

#-----final markers
df<-FetchData(validation,vars = c('ident',markers),layer = 'counts')
levels(df$ident)
table(df$ident)
df$ident<-ifelse(df$ident=='escc.tumor', 1, 0)
table(df$ident)
df<-arrange(df,ident)
df<-df[sample(1:nrow(df)),]
table(df$ident)
test<-df
DoHeatmap(validation,features=c(final))
table(validation@meta.data$gen_celltype)

df<-FetchData(training,vars = c('ident',markers),layer = 'counts')
levels(df$ident)
table(df$ident)
df$ident<-ifelse(df$ident=='escc.tumor', 1, 0)
table(df$ident)
df<-arrange(df,ident)
df<-df[sample(1:nrow(df)),]
table(df$ident)
train <-df
DoHeatmap(training,features=c(final))
table(training@meta.data$gen_celltype)

x_train <- train %>% dplyr::select(-ident)
str(x_train)

x_test  <- test %>% dplyr::select(-ident)
str(x_test)

y_train <- train %>% dplyr::select(ident)
y_test  <-test %>% dplyr::select(ident)

model_keras <- keras_model_sequential()
model_keras %>% layer_dense(
  units              = 16, 
  kernel_initializer = "uniform", 
  activation         = "leaky_relu", 
  input_shape        = ncol(x_train)) %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(
    units              = 16, 
    kernel_initializer = "uniform", 
    activation         = "sigmoid") %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(
    units              = 1, 
    kernel_initializer = "uniform", 
    activation         = "sigmoid") %>%
  compile(
    optimizer = 'adam',
    loss      = 'binary_crossentropy',
    metrics   = c('accuracy')
  )

history <- fit(
  object           = model_keras, 
  x                = as.matrix(x_train), 
  y                = y_train,
  batch_size       = 20, 
  epochs           = 50,
  validation_split = 0.50
)
print(history)
plot(history) 
# Predicted Class
yhat_keras_class_vec <- predict(object = model_keras, x = as.matrix(x_test)) %>%
  as.vector()

# Predicted Class Probability
yhat_keras_prob_vec  <- predict(object = model_keras, x = as.matrix(x_test)) %>%
  as.vector()


# Format test data and predictions 
estimates_keras_tbl <- tibble(
  truth      = as.factor(y_test$ident) %>% fct_recode(yes = "1", no = "0"),
  estimate   = as.factor(ifelse(yhat_keras_class_vec>0.5,'yes','no')) %>% fct_recode(yes = "1", no = "0"),
  class_prob = yhat_keras_prob_vec
)
str(estimates_keras_tbl)

# Confusion Table
#estimates_keras_tbl %>% conf_mat(truth, estimate)
confusion_matrix<-caret::confusionMatrix(estimates_keras_tbl$truth, estimates_keras_tbl$estimate,
                                         positive='yes')
fourfoldplot(as.table(confusion_matrix),
             color = c('grey','skyblue'),
             std='all.max',
             main='Keras tensorflow model predictions')
confusion_matrix

# Accuracy
estimates_keras_tbl %>% metrics(truth, estimate)

# Setup lime::model_type() function for keras
model_type.keras.src.models.sequential.Sequential <- function(x, ...) {
  "classification"
}


# Setup lime::predict_model() function for keras
predict_model.keras.src.models.sequential.Sequential <- function(x, newdata, type, ...) {
  pred <- predict(object = x, x = as.matrix(newdata))
  data.frame(Yes = pred, No = 1 - pred)
}

# Test our predict_model() function
predict_model(x = model_keras, newdata = x_test, type = 'raw') %>%
  tibble::as_tibble()

explainer <- lime(x_train, model_keras)

explanation <- lime::explain(x_test[1:6,], explainer, n_labels = 1, n_features = 9)

plot_features(explanation) +
  labs(title = "LIME Feature Importance Visualization Final Genes/Cell",
       subtitle = "Hold Out (Test) Set, First 6 Cases Shown")
markers
break 
