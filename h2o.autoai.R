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
aml <- h2o.automl(x, y,
                  training_frame = train_obs.h2o,
                  max_models = 20,
                  seed = 1)
lb <- h2o.get_leaderboard(object = aml, extra_columns = "ALL")
lb<-as.data.frame(lb)
head(lb)
best.model <- h2o.getModel("GBM_3_AutoML_1_20240930_104436")
best.model
break 
#----plot tree
Tree = h2o.getModelTree(model = best.model, tree_number = 1)
createDataTree <- function(h2oTree) {
  
  h2oTreeRoot = h2oTree@root_node
  
  dataTree = Node$new(h2oTreeRoot@split_feature)
  dataTree$type = 'split'
  
  addChildren(dataTree, h2oTreeRoot)
  
  return(dataTree)
}

addChildren <- function(dtree, node) {
  
  if(class(node)[1] != 'H2OSplitNode') return(TRUE)
  
  feature = node@split_feature
  id = node@id
  na_direction = node@na_direction
  
  if(is.na(node@threshold)) {
    leftEdgeLabel = printValues(node@left_levels, na_direction=='LEFT', 4)
    rightEdgeLabel = printValues(node@right_levels, na_direction=='RIGHT', 4)
  }else {
    leftEdgeLabel = paste("<", node@threshold, ifelse(na_direction=='LEFT',',NA',''))
    rightEdgeLabel = paste(">=", node@threshold, ifelse(na_direction=='RIGHT',',NA',''))
  }
  
  left_node = node@left_child
  right_node = node@right_child
  
  if(class(left_node)[[1]] == 'H2OLeafNode')
    leftLabel = paste("prediction:", left_node@prediction)
  else
    leftLabel = left_node@split_feature
  
  if(class(right_node)[[1]] == 'H2OLeafNode')
    rightLabel = paste("prediction:", right_node@prediction)
  else
    rightLabel = right_node@split_feature
  
  if(leftLabel == rightLabel) {
    leftLabel = paste(leftLabel, "(L)")
    rightLabel = paste(rightLabel, "(R)")
  }
  
  dtreeLeft = dtree$AddChild(leftLabel)
  dtreeLeft$edgeLabel = leftEdgeLabel
  dtreeLeft$type = ifelse(class(left_node)[1] == 'H2OSplitNode', 'split', 'leaf')
  
  dtreeRight = dtree$AddChild(rightLabel)
  dtreeRight$edgeLabel = rightEdgeLabel
  dtreeRight$type = ifelse(class(right_node)[1] == 'H2OSplitNode', 'split', 'leaf')
  
  addChildren(dtreeLeft, left_node)
  addChildren(dtreeRight, right_node)
  
  return(FALSE)
}

printValues <- function(values, is_na_direction, n=4) {
  l = length(values)
  
  if(l == 0)
    value_string = ifelse(is_na_direction, "NA", "")
  else
    value_string = paste0(paste0(values[1:min(n,l)], collapse = ', '),
                          ifelse(l > n, ",...", ""),
                          ifelse(is_na_direction, ", NA", ""))
  
  return(value_string)
}
dtree = createDataTree(Tree)
GetEdgeLabel <- function(node) {return (node$edgeLabel)}
GetNodeShape <- function(node) {switch(node$type, split = "diamond", leaf = "oval")}
SetEdgeStyle(dtree, fontname = 'Palatino', label = GetEdgeLabel, labelfloat = TRUE,fontsize=60,color='blue')
SetNodeStyle(dtree, fontname = 'Palatino', shape = GetNodeShape,fontsize=60,color='red')
SetGraphStyle(dtree, rankdir = "LR", dpi=70.)
plot(dtree, output="graph")
