library(biomaRt)
setwd('/home/amp_prog/rstudio/HLH')
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl <- useMart("ensembl")
listDatasets(ensembl)
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
filters <- listFilters(ensembl)
filters[1:10,]
attributes <- listAttributes(ensembl)
attributes[1:10,]
attributes
data <- read.csv('/home/amp_prog/rstudio/HLH/ensembl.csv')
data <- data[,-1]
#
geneid <- data$ensembl_id
head(geneid)
dim(geneid)
genes <- getBM(attributes = c('ensembl_transcript_id_version','ensembl_gene_id','external_gene_name','description',
                              'transcript_length','transcript_count','ensembl_exon_id'),
               filters = 'ensembl_transcript_id_version',
               values = geneid,
               mart = ensembl)
head(genes)
dim(data)
df <- data[data$ensembl_id %in% unique(genes$ensembl_transcript_id_version),]
dim(df)
genes$ensembl_transcript_id_version
colnames(genes)[1] <- 'ensembl_id'
table <- merge(genes,df,by="ensembl_id")
matrix <- table[!(is.na(table$external_gene_name) | table$external_gene_name==""),]
dim(matrix)
row.names(matrix)<- matrix$external_gene_name

break
row.names(matrix) <- make.names(matrix$external_gene_name, unique = TRUE)
matrix<-matrix[,-c(1:4)]
setwd('/home/amp_prog/rstudio/HLH')
write.csv(matrix, file = 'external_gene_name.csv')
