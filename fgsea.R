library(fgsea)
library(biomaRt)
library(dplyr)
genes<-read.csv('/home/em_b/work_stuff/rnaseq/for_genes_set_neuron_diff_0hrMB.csv')#,header = FALSE)
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
ensembl=useMart("ensembl")
ensembl=useDataset("mmusculus_gene_ensembl",
                   mart = ensembl)
#geneid <- genes$V1
geneid <- genes$Gene

head(geneid)
genes <-getBM(attributes = c('external_gene_name','entrezgene_id'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
matrix<-read.csv('/home/em_b/work_stuff/rnaseq/mm39_entrez.matrix.csv',
                 row.names = 1)
matrix<-matrix[,-1]
#matrix<-matrix[,-c(4:9,13:33)]
#matrix<-matrix[,-c(1:3,7:9,13:33)]
matrix<-matrix[,-c(16:33)]

plotCoregulationProfile(pathway=genes$entrezgene_id, 
                        E=matrix,
                        scale = TRUE)
summary(matrix)
#-------------------------------------------------------------------------------
genes<-read.csv('/home/em_b/work_stuff/rnaseq/for_genes_set_neuron_diff_0hrMB.csv')#,header = FALSE)
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
ensembl=useMart("ensembl")
ensembl=useDataset("hsapiens_gene_ensembl",
                   mart = ensembl)
#geneid <- genes$V1
geneid <- genes$Gene

head(geneid)
genes <-getBM(attributes = c('external_gene_name','entrezgene_id'),
              filters = 'external_gene_name',
              values = geneid,
              mart = ensembl)
head(genes)
matrix<-read.csv('/home/em_b/work_stuff/rnaseq/hg38_entrez.matrix.csv',
                 row.names = 1)
matrix<-matrix[,-1]
#matrix<-matrix[,-c(4:9,13:33)]
#matrix<-matrix[,-c(1:3,7:9,13:33)]
#matrix<-matrix[,-c(16:33)]

plotCoregulationProfile(pathway=genes$entrezgene_id, 
                        E=matrix,
                        scale = TRUE)
summary(matrix)
