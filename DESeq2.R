library(DESeq2)
setwd('/home/deviancedev/Desktop/drive_nov2023/FCCC/alignments/finals')
matrix<-read.csv('gene.list_ordered.csv')
matrix<-matrix[,-c(2:7)]
row.names(matrix)<-make.names(matrix$Gene,unique=TRUE)
matrix<-matrix[,-1]
colnames(matrix)[2]<-'Norm_2'
head(matrix)
str(matrix)
matrix<-round(matrix)
#-------------------------------------------------------------------------------
names<-c('Norm_1','Norm_2','Norm_3','TAA_1','TAA_2','TAA_3')
#-----set conditions to A and B. A will be the reference condition. 
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
#-------------------------------------------------------------------------------
results<-results(DE)
results
df<-data.frame(results)
df<-df[order(df$pvalue, decreasing=FALSE),]
df<-cbind(row.names(df),df)
df$`row.names(df)`<-sub('\\..*','',df$`row.names(df)`)
colnames(df)[1]<-'Gene'
df<-df[-c(20328:115163),]
write.csv(df,file = 'DESeq2_out.csv')
