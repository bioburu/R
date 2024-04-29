library(DESeq2)
library(pheatmap)
library(dplyr)
setwd('/home/em_b/work_stuff/FCCC/T3_RNAseq/bulk_RNAseq_raw_files')
matrix<-read.csv('biomart.mm39.csv')
str(matrix)
row.names(matrix)<-make.names(matrix$Gene,unique = TRUE)
matrix<-matrix[,-1]
head(matrix)
gnp_0h_v_pbs<-matrix[,-c(1:18,25:27)]
gnp_0h_v_pbs<-round(gnp_0h_v_pbs)
#---remove all zeros
#-------------------------------------------------------------------------------
names<-colnames(gnp_0h_v_pbs)
condition<-c('A','A','A','B','B','B')
type<-c('paired','paired','paired','paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
row.names(coldata)<-coldata$names
coldata<-coldata[,-1]
head(coldata)
str(gnp_0h_v_pbs)
summary(gnp_0h_v_pbs)
#-----------------------------------------------------------------------------
deseq<-DESeqDataSetFromMatrix(countData = gnp_0h_v_pbs,colData = coldata,design = ~ condition)
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
#------------------------------------
library(biomaRt)
library(dplyr)
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useDataset("mmusculus_gene_ensembl",
                   mart = ensembl)
filters = listFilters(ensembl)
attributes=listAttributes(ensembl)
#-----------------------------------------------------------------
geneid <- df$Gene
head(geneid)
listFilters(ensembl)
listAttributes(ensembl)
genes <-getBM(attributes = c('external_gene_name','entrezgene_id'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
genes$external_gene_name
colnames(genes)[1]<- 'Gene'
colnames(df)[1]<-'Gene'
Ftable <- merge(genes, df, by="Gene")

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
