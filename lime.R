# required packages
# install vip from github repo: devtools::install_github("koalaverse/vip")
library(lime)       # ML local interpretation
library(vip)        # ML global interpretation
library(pdp)        # ML global interpretation
library(ggplot2)    # visualization pkg leveraged by above packages
library(caret)      # ML model building
library(h2o)  
library(rsample)# ML model building
# initialize h2o
h2o.init()
#h2o.no_progress()
# create data sets
library(modeldata)
data("attrition", package = "modeldata")
df<-attrition
#---unorder ordered factors and change attrition to factor class w levels
df <- df %>% 
  dplyr::mutate_if(is.ordered, factor, ordered = FALSE) %>%
  dplyr::mutate(Attrition = factor(Attrition, levels = c("Yes", "No")))
#---set amount of sample size to save for validation observations
index <- 1:5
#----subset training observations
train_obs <- df[-index, ]
#--subset validation observations
local_obs <- df[index, ]
# create h2o objects for modeling
#------split variables by columns
y <- "Attrition"
y
x <- setdiff(names(train_obs), y)
x
train_obs.h2o <- as.h2o(train_obs)
train_obs.h2o
local_obs.h2o <- as.h2o(local_obs)
local_obs.h2o
#-------------------------------------------------------------------------------
# Create Random Forest model with ranger via caret
fit.caret <- train(
  Attrition ~ ., 
  data = train_obs, 
  method = 'ranger',
  trControl = trainControl(method = "cv", number = 5, classProbs = TRUE),
  tuneLength = 1,
  importance = 'impurity'
)
fit.caret

explainer_caret <- lime(train_obs, fit.caret, n_bins = 5)
explainer_caret
class(explainer_caret)
summary(explainer_caret)
local_obs$Attrition

explanation_caret <-explain(
  local_obs,
  explainer_caret,
  labels = 'Yes',
  n_labels = NULL,
  n_features=10,
  n_permutations = 5000,
  feature_select = "highest_weights",
  dist_fun = "gower",
  kernel_width = .75,
  gower_pow = 1)
tibble::glimpse(explanation_caret)
plot_features(explanation_caret)
plot_explanations(explanation_caret)
explanation_caret <-explain(
  local_obs,
  explainer_caret,
  labels = 'Yes',
  n_labels = NULL,
  n_features=10,
  n_permutations = 5000,
  feature_select = "lasso_path",
  dist_fun = "gower",
  kernel_width = 3,
  gower_pow = 1)
plot_features(explanation_caret)
#-------------------------------------------------------------------------------
# create h2o models
h2o_rf <- h2o.randomForest(x, y, training_frame = train_obs.h2o)
h2o_rf
h2o_glm <- h2o.glm(x, y, training_frame = train_obs.h2o, family = "binomial")
h2o_glm
h2o_gbm <- h2o.gbm(x, y, training_frame = train_obs.h2o)
h2o.gbm
# built-in PDP support in H2O
h2o.partialPlot(h2o_rf, data = train_obs.h2o, cols = "MonthlyIncome")

#-------------------------------------------------------------------------------
# ranger model --> model type not built in to LIME
fit.ranger <- ranger::ranger(
  Attrition ~ ., 
  data = train_obs, 
  importance = 'impurity',
  probability = TRUE
)
fit.ranger
vip(fit.ranger) + ggtitle("ranger: RF")
fit.ranger %>%
  partial(pred.var = "MonthlyIncome", grid.resolution = 25, ice = TRUE) %>%
  autoplot(rug = TRUE, train = train_obs, alpha = 0.1, center = TRUE)
#-------------------------------------------------------------------------------
break
explainer_h2o_rf  <- lime(train_obs, h2o_rf, n_bins = 5)
explainer_h2o_glm <- lime(train_obs, h2o_glm, n_bins = 5)
explainer_h2o_gbm <- lime(train_obs, h2o_gbm, n_bins = 5)

explanation_rf <- explain(local_obs, explainer_h2o_rf, n_features = 5, labels = "Yes", kernel_width = .1, feature_select = "highest_weights")
explanation_glm <- explain(local_obs, explainer_h2o_glm, n_features = 5, labels = "Yes", kernel_width = .1, feature_select = "highest_weights")
explanation_gbm <- explain(local_obs, explainer_h2o_gbm, n_features = 5, labels = "Yes", kernel_width = .1, feature_select = "highest_weights")

p1 <- plot_features(explanation_rf, ncol = 1) + ggtitle("rf")
p2 <- plot_features(explanation_glm, ncol = 1) + ggtitle("glm")
p3 <- plot_features(explanation_gbm, ncol = 1) + ggtitle("gbm")
gridExtra::grid.arrange(p1, p2, p3, nrow = 1)

#----make new model type to use with lime
explainer_ranger <- lime(train, fit.ranger, n_bins = 5)
## Error in UseMethod("lime", x): no applicable method for 'lime' applied to an object of class "function"
# get the model class
class(fit.ranger)
## [1] "ranger"

# need to create custom model_type function
model_type.ranger <- function(x, ...) {
  # Function tells lime() what model type we are dealing with
  # 'classification', 'regression', 'survival', 'clustering', 'multilabel', etc
  
  return("classification")
}

model_type(fit.ranger)
## [1] "classification"
# need to create custom predict_model function
predict_model.ranger <- function(x, newdata, ...) {
  # Function performs prediction and returns data frame with Response
  pred <- predict(x, newdata)
  return(as.data.frame(pred$predictions))
}

predict_model(fit.ranger, newdata = local_obs)
##          Yes        No
## 1 0.27451508 0.7254849
## 2 0.08705952 0.9129405
## 3 0.44530397 0.5546960
## 4 0.32226270 0.6777373
## 5 0.23780397 0.7621960
explainer_ranger <- lime(train_obs, fit.ranger, n_bins = 5)
explanation_ranger <- explain(local_obs, explainer_ranger, n_features = 5, n_labels = 2, kernel_width = .1)
plot_features(explanation_ranger, ncol = 2) + ggtitle("ranger")
