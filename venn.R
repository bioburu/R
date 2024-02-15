library(dplyr)
#-----------------
bulk_tissue<- read.csv('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/bulk_tissue/csv_results/bulk_tissues_deseq.edgeR.csv')
tissue_exosome<- read.csv('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/exosomes/csv_results/exosomes_tissues_deseq.edgeR.csv')
plasma_exosome<-read.csv('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/exosomes/csv_results/exosomes_plasma_deseq.edgeR.csv')
#--------------------
library(ggVennDiagram)
venn<-list(bulk_tissue$R_trackingID,tissue_exosome$R_trackingID,plasma_exosome$R_trackingID)
#View(venn)
library(ggplot2)
ggVennDiagram(venn,category.names = c('T.cereb','T.cereb_exo','T.pla_exo'),
              label = 'both',set_color = 'red',label_alpha = 0)+scale_fill_distiller(palette='Set3')+
  labs(title='')
#--------
no_3<-semi_join(bulk_tissue,tissue_exosome,by='R_trackingID')
no_3<-semi_join(no_3,plasma_exosome,by='R_trackingID')
no_3$R_trackingID
#---------
x<-semi_join(bulk_tissue,tissue_exosome,by ='R_trackingID')
y<-semi_join(bulk_tissue,plasma_exosome,by ='R_trackingID')
no_8727<-bulk_tissue[!bulk_tissue$R_trackingID %in% x$R_trackingID,]
no_8727<-no_8727[!no_8727$R_trackingID %in% y$R_trackingID,]
#write.csv(no_8727,file='Venn_8727.csv')
#----------
no_749<-semi_join(x,no_3,by ='R_trackingID')
no_749<-x[!x$R_trackingID %in% no_3$R_trackingID,]
#write.csv(no_749,file='Venn_749.csv')
#----------
no_4<-y[!y$R_trackingID %in% no_3$R_trackingID,]
#write.csv(no_4,file='Venn_4.csv')
#---------
x<-semi_join(tissue_exosome,plasma_exosome,by='R_trackingID')
y<-semi_join(bulk_tissue,tissue_exosome,by ='R_trackingID')
no_936<-tissue_exosome[!tissue_exosome$R_trackingID %in% x$R_trackingID,]
no_936<-no_936[!no_936$R_trackingID %in% y$R_trackingID,]
#write.csv(no_936,file='Venn_936.csv')
x<-semi_join(tissue_exosome,plasma_exosome,by='R_trackingID')
no_8<-x[!x$R_trackingID %in% no_3$R_trackingID,]
#write.csv(no_8,file='Venn_8.csv')
#------------------------
y<-semi_join(bulk_tissue,plasma_exosome,by='R_trackingID')
no_14<-plasma_exosome[!plasma_exosome$R_trackingID %in% x$R_trackingID,]
no_14<-no_14[!no_14$R_trackingID %in% y$R_trackingID,]
#write.csv(no_14,file='Venn_14.csv')

break 
write.csv(final_join,file='venn_results.csv')
