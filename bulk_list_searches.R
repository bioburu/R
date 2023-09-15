library(dplyr)
setwd('/home/amp_prog/Desktop/models/GSE242889_liver_cancer')
list1<-read.csv('pt1.csv')
list2<-read.csv('pt2.csv')
list3<-read.csv('pt3.csv')
list4<-read.csv('pt4.csv')
list5<-read.csv('pt5.csv')
#--------------------------------------------------------------------------------
x<-semi_join(list1, list2, by = c('x'))
x<-semi_join(x, list3, by = c('x'))
x<-semi_join(x, list4, by = c('x'))
x<-semi_join(x, list5, by = c('x'))
write.csv(x,file = 'gene_list.csv')

