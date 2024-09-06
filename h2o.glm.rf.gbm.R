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
df$SeniorCitizen<-NULL
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
#----logistic regression
h2o_glm <- h2o.glm(x, y, training_frame = train_obs.h2o, family = "binomial")
h2o_glm
#---partial dependency plot for regression
h2o.partialPlot(h2o_glm,
                newdata = train_obs.h2o,
                cols = "MonthlyCharges",
                plot=TRUE,
                plot_stddev=TRUE)
explainer_h2o_glm <- lime(train_obs,
                  h2o_glm,
                  bin_continuous = TRUE,
                  quantile_bins = FALSE)
explanation_glm <- lime::explain(local_obs[1:6,],
                                 explainer_h2o_glm,
                             n_labels = 1,
                             n_features = 3)
h2o.confusionMatrix(h2o_glm, local_obs.h2o, valid = FALSE, xval = FALSE)
plot_features(explanation_glm)
plot_explanations(explanation_glm)
#----randomforest
h2o_rf <- h2o.randomForest(x, y,
                           nfolds = 5,
                           training_frame = train_obs.h2o)
h2o_rf
h2o.partialPlot(h2o_rf,
                newdata = train_obs.h2o,
                cols = "MonthlyCharges")
explainer_h2o_rf  <- lime(train_obs,
                          h2o_rf,
                          bin_continuous = TRUE,
                          quantile_bins = FALSE)
explanation_rf <- lime::explain(local_obs[1:6,],
                                explainer_h2o_rf,
                                n_labels = 1,
                                n_features = 3)
h2o.confusionMatrix(h2o_rf,local_obs.h2o, valid = FALSE, xval = FALSE)
plot_features(explanation_rf)
plot_explanations(explanation_rf)

#-----gradient boosting tree model
h2o_gbm <- h2o.gbm(x, y,
                   training_frame = train_obs.h2o,
                   nfolds = 5)
head(h2o.gbm)
h2o.partialPlot(h2o_gbm,
                data = train_obs.h2o,
                cols = "MonthlyCharges")
explainer_h2o_gbm  <- lime(train_obs,
                          h2o_gbm,
                          bin_continuous = TRUE,
                          quantile_bins = FALSE)
explanation_gbm <- lime::explain(local_obs[1:6,],
                                 explainer_h2o_gbm,
                                 n_labels = 1,
                                 n_features = 3)
plot_features(explanation_gbm)
plot_explanations(explanation_gbm)

p1 <- plot_features(explanation_rf[4:6,], ncol = 1) + ggtitle("rf")
p2 <- plot_features(explanation_glm[4:6,], ncol = 1) + ggtitle("glm")
p3 <- plot_features(explanation_gbm[4:6,], ncol = 1) + ggtitle("gbm")
gridExtra::grid.arrange(p1, p2, p3, nrow = 1)
