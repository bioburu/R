library(Rsubread)
library(biomaRt)
counts<-featureCounts(files = '/Users/burudpc/Desktop/geo_t3_rnaseq/SRR23386663.rmdup.sort.chr1.19.bam',
                      annot.inbuilt = 'mm39',
                      isPairedEnd = TRUE)
matrix<-data.frame(counts$counts)
matrix<-cbind(row.names(matrix),matrix)
colnames(matrix)<-c('Gene','MB.T3_1')
str(matrix)
#--------------------biomart
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl=useMart('ENSEMBL_MART_ENSEMBL')
listDatasets(ensembl)
ensembl=useDataset("mmusculus_gene_ensembl",
                   mart = ensembl)
geneid <- matrix$Gene
head(geneid)
listFilters(ensembl)
listAttributes(ensembl)
genes <-getBM(attributes = c('external_gene_name','entrezgene_id'),
              filters = 'entrezgene_id',
              values = geneid,
              mart = ensembl)
head(genes)
colnames(genes)[2]<-'Gene'
colnames(genes)
output <- merge(genes, matrix, by='Gene')
setwd('/Users/burudpc/Desktop')
write.csv(output,file='SRR23386662_output.csv')

