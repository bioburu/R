# required packages
# install vip from github repo: devtools::install_github("koalaverse/vip")
library(lime)       # ML local interpretation
library(vip)        # ML global interpretation
library(pdp)        # ML global interpretation
library(ggplot2)    # visualization pkg leveraged by above packages
library(caret)      # ML model building
library(h2o)  
library(rsample)# ML model building
library(dplyr)
# initialize h2o
h2o.init()
#h2o.no_progress()
# create data sets
df<-read.csv('/home/deviancedev01/Desktop/interesting_exps/model_explaining/WA_Fn-UseC_-Telco-Customer-Churn.csv',
             stringsAsFactors = TRUE)%>%mutate_if(is.numeric,round)%>%mutate_if(is.numeric, as.integer)
df$customerID<-NULL
str(df)
#---unorder ordered factors and change attrition to factor class w levels
#df <- df %>% 
#  dplyr::mutate_if(is.ordered, factor, ordered = FALSE) %>%
#  dplyr::mutate(Attrition = factor(Attrition, levels = c("Yes", "No")))
#str(df)
#---set amount of sample size to save for validation observations
index <- 1:200
#----subset training observations
train_obs <- df[-index, ]
str(train_obs)
#--subset validation observations
local_obs <- df[index, ]
str(local_obs)
# create h2o objects for modeling
#------split variables by columns
y <- "Churn"
y
x <- setdiff(names(train_obs), y)
x
train_obs.h2o <- as.h2o(train_obs)
train_obs.h2o
local_obs.h2o <- as.h2o(local_obs)
local_obs.h2o
#----randomforest
h2o_rf <- h2o.randomForest(x, y, training_frame = train_obs.h2o)
h2o_rf
#----logistic regression
h2o_glm <- h2o.glm(x, y, training_frame = train_obs.h2o, family = "binomial")
h2o_glm
#-----gradient boosting tree model
h2o_gbm <- h2o.gbm(x, y, training_frame = train_obs.h2o)
h2o.gbm
#---partial dependency plot for random forest
h2o.partialPlot(h2o_rf, data = train_obs.h2o, cols = "MonthlyCharges")
#---partial dependency plot for regression
h2o.partialPlot(h2o_glm, data = train_obs.h2o, cols = "MonthlyCharges")
#---partial dependency plot for gradient boosting tree model
h2o.partialPlot(h2o_gbm, data = train_obs.h2o, cols = "MonthlyCharges")

explainer_h2o_rf  <- lime::lime(train_obs, h2o_rf, n_bins = 5)
explainer_h2o_glm <- lime::lime(train_obs, h2o_glm, n_bins = 5)
explainer_h2o_gbm <- lime::lime(train_obs, h2o_gbm, n_bins = 5)

explanation_rf <- lime::explain(local_obs[1:2,], explainer_h2o_rf, n_features = 5, labels = "Yes", kernel_width = .1, feature_select = "highest_weights")
explanation_glm <- lime::explain(local_obs[1:2,], explainer_h2o_glm, n_features = 5, labels = "Yes", kernel_width = .1, feature_select = "highest_weights")
explanation_gbm <- lime::explain(local_obs[1:2,], explainer_h2o_gbm, n_features = 5, labels = "Yes", kernel_width = .1, feature_select = "highest_weights")

p1 <- plot_features(explanation_rf, ncol = 1) + ggtitle("rf")
p2 <- plot_features(explanation_glm, ncol = 1) + ggtitle("glm")
p3 <- plot_features(explanation_gbm, ncol = 1) + ggtitle("gbm")
gridExtra::grid.arrange(p1, p2, p3, nrow = 1)
