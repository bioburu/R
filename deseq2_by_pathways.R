library(DESeq2)
library(pheatmap)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(dplyr)
library(RColorBrewer)
setwd('/home/em_b/work_stuff/FCCC/T3_RNAseq')
matrix<-read.csv('biomart.mm39.csv')
go_neural_diff<-read.delim('/home/em_b/Downloads/GO_term_summary_20240502_214014.txt',
                 row.names = NULL)
list<-go_neural_diff$MGI.Gene.Marker.ID
list
colnames(go_neural_diff)
matrix<-subset(matrix, subset = Gene %in% list)
str(matrix)
row.names(matrix)<-make.names(matrix$Gene,unique = TRUE)
matrix<-matrix[,-1]
head(matrix)
tester<-matrix[,-c(7:27)]
colnames(tester)
colnames(tester)<-c('PBS','PBS.1','PBS.2',
                    'T3','T3.1','T3.2')
names<-colnames(tester)
names
condition<-c('A','A','A','B','B','B')
type<-c('paired','paired','paired','paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
row.names(coldata)<-make.names(coldata$names,
                               unique=TRUE)
coldata<-coldata[,-1]
coldata
#-----------------------------------------------------------------------------
tester<-round(tester)
colnames(tester)
row.names(coldata)
deseq<-DESeqDataSetFromMatrix(countData = tester,colData = coldata,design = ~ condition)
deseq
#----if using standard negbinomial walds test  
DE<-DESeq(deseq)
plotMA(DE,ylim=c(-5,5))
plotDispEsts(DE)
#-------------------------------------------------------------------------------
results<-results(DE)
results
x<-data.frame(results)
x<-x[order(x$pvalue, decreasing=FALSE),]
x<-subset(x,pvalue< 0.05)
x<-subset(x,log2FoldChange> 0)
x
x
names<-row.names(x)
str(names)
names<-append(names, 'Neurod1', after = 0)
pheatmap(tester[names,],
         scale = 'row',
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'grey',
         fontsize = 17,
         color = colorRampPalette(rev(brewer.pal(n=7,name='RdYlBu')))(100))
break 
#----if using LRT---------------------------------------------- 
deseq<-estimateSizeFactors(deseq)
deseq<-estimateDispersions(deseq)
test <- DESeq(deseq, test="LRT", reduced= ~ 1)
resLRT <- results(test)
resLRT<-data.frame(resLRT)
resLRT
resLRT<-resLRT[order(resLRT$pvalue, decreasing=FALSE),]
resLRT<-subset(resLRT,pvalue< 0.05)
resLRT
