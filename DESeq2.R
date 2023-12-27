library(DESeq2)
library(pheatmap)
library(dplyr)
setwd('/home/deviancedev/Desktop/drive_nov2023/FCCC/In_house_astrocyte_rnaseq/results')
matrix<-read.csv('biomaRt_out.csv')
colnames(matrix)[1]<-'Gene'
row.names(matrix)<-make.names(matrix$Gene,unique=TRUE)
matrix<-matrix[,-c(1:8)]
colnames(matrix)[2]<-'Norm_2'
head(matrix)
str(matrix)
matrix<-round(matrix)
head(matrix)
#---remove all zeros
matrix<-filter(matrix,Norm_1>0,Norm_2>0,Norm_3>0,TAA_1>0,TAA_2>0,TAA_3>0)
head(matrix)
#-------------------------------------------------------------------------------
names<-colnames(matrix)
condition<-c('A','A','A','B','B','B')
type<-c('paired','paired','paired','paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
row.names(coldata)<-coldata$names
coldata<-coldata[,-1]
head(coldata)
#-----------------------------------------------------------------------------
deseq<-DESeqDataSetFromMatrix(countData = matrix,colData = coldata,design = ~ condition)
deseq
DE<-DESeq(deseq)
plotMA(DE,ylim=c(-2,2))
plotDispEsts(DE)
#-------------------------------------------------------------------------------
results<-results(DE)
results
x<-data.frame(results)
x<-x[order(x$padj, decreasing=FALSE),]
x<-subset(x,padj< 0.01)
summary(x)
downreg<-subset(x,log2FoldChange< -1.5)
upreg<-subset(x,log2FoldChange> 1.5)
df<-rbind(upreg,downreg)
df<-cbind(row.names(df),df)
df$`row.names(df)`<-sub('\\..*','',df$`row.names(df)`)
colnames(df)[1]<-'Gene'
df<-cbind(row.names(df),df)
colnames(df)[1]<-'DESeq2_id'
biomart<-read.csv('biomaRt_out.csv')
df <- merge(df, biomart,by="DESeq2_id")
colnames(df)[2]<-'Gene'
colnames(df)[17]<-'Norm_2'
df<-df[,-9]
head(df)
break 
write.csv(df,file = 'In_house_DESeq2.csv')
#-----Wnt genes 
up<-c('Lgr5','Senp2.3','Tle2.13','Ccnd1','Rock2.2','Camk2d.14','Camk2b',
      'Ppp3ca.6','Ppp3r1.1')
up_wnt<-pheatmap(matrix[up,],scale = 'row',border_color = 'black',fontsize = 12,color = c('yellow','orange','red'))
down<-c('Sfrp1','Frzb','Lgr6','Fzd2','Dvl2','Nkd1.1','Tcf7l1',
        'Ctnnbip1','Daam2.2','Mapk9.6')
down_wnt<-pheatmap(matrix[down,],scale = 'row',border_color = 'black',fontsize = 8,color = c('yellow','orange','red'))
both<-c('Lgr5','Senp2.3','Tle2.13','Ccnd1','Rock2.2','Camk2d.14','Camk2b',
        'Ppp3ca.6','Ppp3r1.1','Sfrp1','Frzb','Lgr6','Fzd2','Dvl2','Nkd1.1','Tcf7l1',
        'Ctnnbip1','Daam2.2','Mapk9.6')
wnt<-pheatmap(matrix[both,],scale = 'row',border_color = 'black',fontsize = 8,color = c('yellow','orange','red'))
