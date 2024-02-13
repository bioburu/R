library(DESeq2)
library(pheatmap)
library(dplyr)
matrix<-read.csv('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/bulk_tissue/csv_results/all_conditions_counts.csv',row.names = 1)
tissues<-matrix[,-c(1:6,10:15)]
names<-colnames(tissues)
names
condition<-c('A','A','A','B','B','B')
type<-c('paired','paired','paired','paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
row.names(coldata)<-coldata$names
coldata<-coldata[,-1]
head(coldata)
#-----------------------------------------------------------------------------
deseq<-DESeqDataSetFromMatrix(countData = tissues,colData = coldata,design = ~ condition)
deseq
DE<-DESeq(deseq)
plotMA(DE,ylim=c(-2,2))
plotDispEsts(DE)
#-------------------------------------------------------------------------------
results<-results(DE)
results
x<-data.frame(results)
summary(x)
x<-na.omit(x)
#----only upregulated genes 
downreg<-subset(x,log2FoldChange< -1.5)
upreg<-subset(x,log2FoldChange> 1.5)
df<-rbind(upreg,downreg)
df<-cbind(row.names(df),df)
df$`row.names(df)`<-sub('\\..*','',df$`row.names(df)`)
colnames(df)[1]<-'external_gene_name'
df<-cbind(row.names(df),df)
colnames(df)[1]<-'R_trackingID'
biomart_output<-read.csv('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/bulk_tissue/csv_results/bulk_tissues_biomart.csv',row.names = 1)
df <- merge(df, biomart_output,by="R_trackingID")
df<-df[,-11]
colnames(df)[2]<-'external_gene_name'
df<-df[order(df$padj, decreasing=FALSE),]
df<-subset(df,padj< 0.05)
df<-df[,-12]
#setwd('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/bulk_tissue/csv_results')
#write.csv(df,file = 'bulk_tissues_deseq2.csv')
#--------------
library(pathview)
foldchanges = df$log2FoldChange
names(foldchanges) = df$entrezgene_id
head(foldchanges)
setwd('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/bulk_tissue/kegg')
pv_bulk_tissues <- pathview(gene.data = foldchanges, pathway.id = "04060",
                   species = "mmu", out.suffix = "mmu_test")
break 
