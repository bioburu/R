library(biomaRt)
#set your proxy
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
#list available biomarts
listMarts()
#set mart object
ensembl=useMart("ensembl")
#list the different datasets available from the biomart database
listDatasets(ensembl)
#update the mart object by addition of dataset
ensembl=useDataset("hsapiens_gene_ensembl",mart = ensembl)
#add filter argument
filters = listFilters(ensembl)
filters[1:5,]
#add attributes argument
attributes=listAttributes(ensembl)
attributes[1:5,]
#load in whatever and clean
setwd('/home/amp_prog/rstudio/CAR_T/GSE120649_kidney_rejection')
matrix<-read.delim('GSM3406955_s1.txt.gz',header = FALSE)
dim(matrix)
geneid <-matrix$V1
head(geneid)
#align identifiers
listFilters(ensembl)
listAttributes(ensembl)
genes <-getBM(attributes = c('ensembl_gene_id_version','external_gene_name'),
              filters = 'ensembl_gene_id_version',
              values = geneid,
              mart = ensembl)
head(genes)
df <- matrix[matrix$V1 %in% unique(genes$ensembl_gene_id_version),]
dim(df)
dim(matrix)
genes$ensembl_gene_id_version
colnames(genes)[1]<- 'V1'
matrix <- merge(genes, df, by="V1")
matrix<- matrix[,-c(1)]
colnames(matrix)<-c('Gene','AntibodyMrejection_1')

