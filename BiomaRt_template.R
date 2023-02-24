library(biomaRt)
#set your proxy
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
#list available biomarts
listMarts()
#set mart object
ensembl=useMart("ensembl")
#mouse=useMart("ENSEMBL_MART_MOUSE")
#snp=useMart("ENSEMBL_MART_SNP")
#funcgen=useMart("ENSEMBL_MART_FUNCGEN")
#list the different datasets available from the biomart database
listDatasets(ensembl)
#listDatasets(mouse)
#listDatasets(snp)
#listDatasets(funcgen)
#update the mart object by addition of dataset
ensembl=useDataset("hsapiens_gene_ensembl",mart = ensembl)
#add filter argument
filters = listFilters(ensembl)
filters[1:5,]
#add attributes argument
attributes=listAttributes(ensembl)
attributes[1:5,]
#load in whatever and clean
setwd("/home/amp_prog/Downloads")
matrix <- read.delim("/home/amp_prog/Downloads/GSM6947006_CXCR4-IL10-MSC1_featureCounts.txt",sep = "")
dim(matrix)
View(matrix)
colnames(matrix)<- matrix[1,]
matrix<-matrix[-1,]
geneid <- (matrix$Geneid)
head(geneid)
dim(geneid)
#align identifiers
listFilters(ensembl)
listAttributes(ensembl)
genes <-getBM(attributes = c('ensembl_gene_id','external_gene_name','description'),
      filters = 'ensembl_gene_id',
      values = geneid,
      mart = ensembl)
head(genes)
df <- matrix[matrix$Geneid %in% unique(genes$ensembl_gene_id),]
dim(df)
dim(matrix)
View(df)
genes$ensembl_gene_id
View(genes)
colnames(genes)[1]<- 'Geneid'
Ftable <- merge(genes, df, by="Geneid")
View(Ftable)
Ftable<- Ftable[,-c(4:8)]

