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
# create data sets
df<-read.csv('/home/deviancedev01/Desktop/interesting_exps/model_explaining/WA_Fn-UseC_-Telco-Customer-Churn.csv',
             stringsAsFactors = TRUE)%>%mutate_if(is.numeric,round)%>%mutate_if(is.numeric, as.integer)
df$customerID<-NULL
str(df)
#df<-na.omit(df)
#---unorder ordered factors and change attrition to factor class w levels
#df <- df %>% 
#  dplyr::mutate_if(is.ordered, factor, ordered = FALSE) %>%
#  dplyr::mutate(Churn = factor(Churn, levels = c("Yes", "No")))
#str(df)
#---set amount of sample size to save for validation observations
index <- 1:80
#----subset training observations
train_obs <- df[-index, ]
str(train_obs)
#--subset validation observations
local_obs <- df[index, ]
str(local_obs)
#-------------------------------------------------------------------------------
# Create Random Forest model with ranger via caret
fit.caret <- train(
  Churn ~ ., 
  data = train_obs, 
  method = 'ranger',
  trControl = trainControl(method = "cv", number = 5, classProbs = TRUE),
  tuneLength = 1,
  importance = 'impurity',
  na.action=na.exclude)
fit.caret
class(fit.caret)
model_type(fit.caret)
predict(fit.caret, local_obs)
explainer <- lime(train_obs,
                  fit.caret,
                  bin_continuous = TRUE,
                  quantile_bins = FALSE)
explanation <- lime::explain(local_obs[1:8,],
                       explainer,
                       n_labels = 1,
                       n_features = 9)
plot_features(explanation)
plot_explanations(explanation)
