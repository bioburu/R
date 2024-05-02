library(tidyverse)
library(hrbrthemes)
library(viridis)
setwd('/home/em_b/work_stuff/FCCC/T3_RNAseq/')
matrix<-read.csv('biomart.hg38.csv')
row.names(matrix)<-make.names(matrix$Gene,unique = TRUE)
matrix<-matrix[,-c(1)]
matrix<-data.frame(t(matrix))
classes<-c('SHH_PBS','SHH_PBS','SHH_PBS',
           'SHH_T3','SHH_T3','SHH_T3',
           'G3_PBS','G3_PBS','G3_PBS',
           'G3_T3','G3_T3','G3_T3')

classes
NEUROD1<-matrix$NEUROD1
EZH2<-matrix$EZH2
boxplot<-data.frame(cbind(classes,NEUROD1,EZH2))
str(boxplot)
boxplot$NEUROD1<-as.numeric(boxplot$NEUROD1)
boxplot$EZH2<-as.numeric(boxplot$EZH2)
str(boxplot)
#------Neurod1
boxplot %>%
  ggplot(aes(x=classes, y=NEUROD1, fill=classes)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(
    legend.position="",
    plot.title = element_text(size=11)
  ) +
  ggtitle('Human PDX medulloblastoma tumor cells (n=3/condition)') +
  xlab('Subgroups')
#--------Ezh2
boxplot %>%
  ggplot(aes(x=classes, y=EZH2, fill=classes)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(
    legend.position="",
    plot.title = element_text(size=11)
  ) +
  ggtitle('Human PDX medulloblastoma tumor cells (n=3/condition)') +
  xlab('Subgroups')
