library(DESeq2)
library(pheatmap)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(dplyr)
library(RColorBrewer)
setwd('/home/em_b/work_stuff/FCCC/T3_RNAseq/R')
cat('GO list taken from https://www.informatics.jax.org/go/term/GO:0030182')
matrix<-read.csv('/home/em_b/work_stuff/FCCC/T3_RNAseq/biomart.mm39.csv')
go<-read.delim('/home/em_b/work_stuff/FCCC/T3_RNAseq/R/gene_ontology_list/GO:0030182_mmus_neuron_differentiation.txt',
                 row.names = NULL)
colnames(go)
list<-go$MGI.Gene.Marker.ID
list
colnames(matrix)
dim(matrix)
matrix<-subset(matrix, subset = Gene %in% list)
dim(matrix)
row.names(matrix)<-make.names(matrix$Gene,unique = TRUE)
matrix<-matrix[,-1]
head(matrix)
mb.pbs.t3<-matrix[,-c(7:27)]
colnames(mb.pbs.t3)
names<-colnames(mb.pbs.t3)
names
condition<-c('A','A','A','B','B','B')
type<-c('paired','paired','paired','paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
row.names(coldata)<-make.names(coldata$names,
                               unique=TRUE)
coldata<-coldata[,-1]
coldata
#-----------------------------------------------------------------------------
mb.pbs.t3<-round(mb.pbs.t3)
colnames(mb.pbs.t3)
row.names(coldata)
deseq<-DESeqDataSetFromMatrix(countData = mb.pbs.t3,colData = coldata,design = ~ condition)
deseq
#----if using standard negbinomial walds test  
DE<-DESeq(deseq)
plotMA(DE,ylim=c(-5,5))
plotDispEsts(DE)
#-------------------------------------------------------------------------------
results<-results(DE)
results
x<-data.frame(results)
htmp<-x[order(x$pvalue, decreasing=FALSE),]
htmp<-subset(htmp,pvalue< 0.09)
htmp<-subset(htmp,log2FoldChange> 0)
head(htmp)
go_genes<-row.names(htmp)
str(go_genes)
go_genes<-append(go_genes, 'Neurod1', after = 0)
go_genes
colnames(mb.pbs.t3)<-c('PBS','PBS','PBS',
                    'T3','T3','T3')
pheatmap(mb.pbs.t3[go_genes,],
         scale = 'row',
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'grey',
         fontsize = 9,
         color = colorRampPalette(rev(brewer.pal(n=7,name='RdYlBu')))(100))
nd1<-subset(x,subset=row.names(x) %in% 'Neurod1')
nd1
final<-rbind(nd1,htmp)
final<-final[,-c(1,3,4,6)]
final
write.csv(final,file='R1Q5_neural_differentiation_genes.csv')
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
