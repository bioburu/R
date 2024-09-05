library(MASS)
library(lime)
data(biopsy)
biopsy$ID <- NULL
biopsy <- na.omit(biopsy)
names(biopsy) <- c('clump thickness', 'uniformity of cell size', 
                   'uniformity of cell shape', 'marginal adhesion',
                   'single epithelial cell size', 'bare nuclei', 
                   'bland chromatin', 'normal nucleoli', 'mitoses',
                   'class')
str(biopsy)
set.seed(4)
#-----sample out 
test_set <- sample(seq_len(nrow(biopsy)), 100)
test_set
#--get all dependent variables
prediction <- biopsy$class
prediction
#--remove all dependent variables from dataframe
biopsy$class <- NULL
#---model- enter in dataframes
model <- lda(biopsy[-test_set, ], prediction[-test_set])
model
str(biopsy)
predict(model, biopsy[test_set, ])
explainer <- lime(biopsy[-test_set,],
                  model,
                  bin_continuous = TRUE,
                  quantile_bins = FALSE)
val_set<-biopsy[test_set,]
explanation <- lime::explain(val_set[2:8,],
                       explainer,
                       n_labels = 1,
                       n_features = 9)
plot_features(explanation)
plot_explanations(explanation)
