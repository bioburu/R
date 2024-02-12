library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggridges)
library(caTools)
library(car)
library(caret)
library(pROC)
library(ROCR)
library(readxl)
library(ggVennDiagram)
#----extract all genes expressed in tumor tissues
bulk_tissues<-read.csv('/Users/burudpc/Desktop/exosome_exp/data/bulk_tissues_biomart.csv',row.names = 1)
bulk_tissues<-bulk_tissues[,-c(2:8)]
colnames(bulk_tissues)<-c('Gene','norm_1','norm_2','norm_3','tumor_1','tumor_2','tumor_3')
#--gapdh lowest count is 168
bulk_tissues<-filter(bulk_tissues,norm_1>168,norm_2>168,norm_3>168,tumor_1>168,tumor_2>168,tumor_3>168)
#----extract all genes expressed in tumor-tissue exosomes
tissue_exos<-read.csv('/Users/burudpc/Desktop/exosome_exp/data/exosomes_biomart.csv',row.names = 1)
tissue_exos<-tissue_exos[,-c(2:8,12:20)]
colnames(tissue_exos)[1]<-c('Gene')
tissue_exos<-filter(tissue_exos,A.cereb_1>1,A.cereb_2>1,A.cereb_3>1,
                   T.tiss_1>1,T.tiss_2>1,T.tiss_3>1)
#-----extract all genes expressed in plasma
plasma_exos<-read.csv('/Users/burudpc/Desktop/exosome_exp/data/exosomes_biomart.csv',row.names = 1)
plasma_exos<-plasma_exos[,-c(2:11,15:17,21:23)]
colnames(plasma_exos)[1]<-'Gene'
plasma_exos<-filter(plasma_exos,A.pla_1>1,A.pla_2>1,A.pla_3>1,T.pla_1>1,T.pla_2>1,T.pla_3>1)
#---test overlaps
#---------------
bulk_tissues$Gene<-sub('\\..*','',bulk_tissues$Gene)
setwd('/Users/burudpc/Desktop/exosome_exp/data')
#write.csv(bulk_tissues,file = 'bulk_tiss_all_genes.csv')
tissue_exos$Gene<-sub('\\..*','',tissue_exos$Gene)
setwd('/Users/burudpc/Desktop/exosome_exp/data')
#write.csv(bulk_tissues,file = 'tiss_exos_all_genes.csv')
plasma_exos$Gene<-sub('\\..*','',plasma_exos$Gene)
setwd('/Users/burudpc/Desktop/exosome_exp/data')
#write.csv(plasma_exos,file = 'plasma_exos_all_genes.csv')

#-----------extract all genes again for venn
bulk_tissues<-read.csv('/Users/burudpc/Desktop/exosome_exp/data/bulk_tissues_biomart.csv',row.names = 1)
bulk_tissues<-bulk_tissues[,-c(2:8)]
colnames(bulk_tissues)<-c('Gene','norm_1','norm_2','norm_3','tumor_1','tumor_2','tumor_3')
#--gapdh lowest count is 168
bulk_tissues<-filter(bulk_tissues,norm_1>168,norm_2>168,norm_3>168,tumor_1>168,tumor_2>168,tumor_3>168)
#----extract all genes expressed in tumor-tissue exosomes
tissue_exos<-read.csv('/Users/burudpc/Desktop/exosome_exp/data/exosomes_biomart.csv',row.names = 1)
tissue_exos<-tissue_exos[,-c(2:8,12:20)]
colnames(tissue_exos)[1]<-c('Gene')
tissue_exos<-filter(tissue_exos,A.cereb_1>1,A.cereb_2>1,A.cereb_3>1,
                    T.tiss_1>1,T.tiss_2>1,T.tiss_3>1)
#-----extract all genes expressed in plasma
plasma_exos<-read.csv('/Users/burudpc/Desktop/exosome_exp/data/exosomes_biomart.csv',row.names = 1)
plasma_exos<-plasma_exos[,-c(2:11,15:17,21:23)]
colnames(plasma_exos)[1]<-'Gene'
plasma_exos<-filter(plasma_exos,A.pla_1>1,A.pla_2>1,A.pla_3>1,T.pla_1>1,T.pla_2>1,T.pla_3>1)

global_venn<-list(bulk_tissues$Gene,tissue_exos$Gene,plasma_exos$Gene)

ggVennDiagram(global_venn,category.names = c('bulk','tiss.exo','plas.exo'),
              label = 'both',set_color = c('red'),edge_size = 1,
              label_size = 4,label_alpha = 0)+scale_fill_distiller(palette='Set3')+
  labs(title='')
break
#------pathway and GO enrichement comparisons
bulk_enr<-read_excel('/Users/burudpc/Desktop/exosome_exp/data/global_gene_enrich_MB_tissue_only.xlsx')
tiss.exo_enr<-read_excel('/Users/burudpc/Desktop/exosome_exp/data/global_gene_enrich_MB_tissue_exosomes_only.xlsx')
plas.exo_enr<-read_excel('/Users/burudpc/Desktop/exosome_exp/data/global_gene_enrich_MB_plasma_exosomes.xlsx')

enr_3617<-bulk_enr[!bulk_enr$Term %in% tiss.exo_enr$Term,]
enr_3617<-enr_3617[!enr_3617$Term %in% plas.exo_enr$Term,]

enr_2387<-merge(bulk_enr,tiss.exo_enr,by=c('Term'))
enr_2387<-enr_2387[!enr_2387$Term %in% plas.exo_enr$Term,]

enr_679<-merge(bulk_enr,plas.exo_enr,by=c('Term'))
enr_679<-enr_679[!enr_679$Term %in% tiss.exo_enr$Term,]

enr_3378<-merge(bulk_enr,tiss.exo_enr,by=c('Term'))
enr_3378<-merge(enr_3378,plas.exo_enr,by=c('Term'))

enr_599<-merge(tiss.exo_enr,plas.exo_enr,by=c('Term'))
enr_599<-enr_599[!enr_599$Term %in% bulk_enr$Term,]

enr_1934<-tiss.exo_enr[!tiss.exo_enr$Term %in% bulk_enr$Term,]
enr_1934<-enr_1934[!enr_1934$Term %in% plas.exo_enr$Term,]

enr_1066<-plas.exo_enr[!plas.exo_enr$Term %in% bulk_enr$Term,]
enr_1066<-enr_1066[!enr_1066$Term %in% tiss.exo_enr$Term,]
