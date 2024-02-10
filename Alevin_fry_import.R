#remotes::install_github("mikelove/fishpond")
library(fishpond)
library(SummarizedExperiment)
library(biomaRt)
library(dplyr)
library(Seurat)
files<-'/home/deviancedev/Desktop/drive_jan2024/FCCC/scRNAseq_test/SRR27442073/mm10_quant'
x<-loadFry(files,outputFormat = 'scRNA',nonzero = FALSE,quiet = FALSE)
x
#----
gc()
counts <- assay(x, "counts")
counts
dim(counts)
df<-data.frame(counts)
gc()
#------
matrix<-df
matrix<-cbind(row.names(matrix),matrix)
colnames(matrix)[1]<-'Gene'
#-----
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl<-useEnsembl(biomart = 'genes',
                   dataset = 'mmusculus_gene_ensembl',
                   version=102)
filters = listFilters(ensembl)
attributes=listAttributes(ensembl)
#-----------------------------------------------------------------
dim(matrix)
geneid <- matrix$Gene
head(geneid)
dim(geneid)
listFilters(ensembl)
listAttributes(ensembl)
genes <-getBM(attributes = c('ensembl_gene_id','external_gene_name','chromosome_name','p_value','start_position','end_position','description'),
              filters = 'ensembl_gene_id',
              values = geneid,
              mart = ensembl)
head(genes)
genes<-genes[!duplicated(genes$ensembl_gene_id), ]
df <- matrix[matrix$Gene %in% unique(genes$ensembl_gene_id),]
dim(df)
dim(matrix)
colnames(genes)[1]<- 'Gene'
Ftable <- merge(genes, df, by="Gene")
gc()
colnames(Ftable)[1] <- c('ensembl_id')
rm(attributes,counts,df,ensembl,filters,genes,matrix,x)
gc()
DF<-Ftable[!duplicated(Ftable$external_gene_name), ]
break
setwd('/home/deviancedev/Desktop/drive_jan2024/FCCC/scRNAseq_test/SRR27442073')
write.csv(Ftable,file = 'GSM8004749.Sciatic.nerve.CD45p.spon.autoimmune.neuropathy.csv')
