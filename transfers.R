library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(caTools)
library(car)
library(caret)
library(InformationValue)
library(pROC)
library(ROCR)
library(readxl)
#----extract all genes expressed in tumor tissues
tumors<-read.csv('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/bulk_tissue/csv_results/bulk_tissues_biomart.csv')
tumors<-tumors[,-c(2:12)]
colnames(tumors)<-c('Gene','tumor_1','tumor_2','tumor_3')
tumors<-filter(tumors,tumor_1>2000,tumor_2>2000,tumor_3>2000)
#----extract all genes expressed in tumor-tissue exosomes
t_exosomes<-read.csv('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/exosomes/csv_results/exosomes_biomart.csv')
t_exosomes<-t_exosomes[,-c(2:21)]
colnames(t_exosomes)<-c('Gene','tumor_1','tumor_2','tumor_3')
t_exosomes<-filter(t_exosomes,tumor_1>2000,tumor_2>2000,tumor_3>2000)
#---test overlaps
x<-semi_join(tumors,t_exosomes,by=c('Gene'))
library(ggVennDiagram)
venn<-list(tumors$Gene,t_exosomes$Gene)
ggVennDiagram(venn,category.names = c(' ',' '),
              label = 'both',set_color = 'red',edge_size = 1,
              label_size = 4,label_alpha = 0)+scale_fill_distiller(palette='Set3')+
  labs(title='')
#-----find overlapping genes
df<-merge(tumors,t_exosomes,by=c('Gene'))
df$Gene<-sub('\\..*','',df$Gene)
setwd('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/exosomes/csv_results')
#write.csv(df,file='global_transcript_comparison.csv')
break 
