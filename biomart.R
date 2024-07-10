library(biomaRt)
library(dplyr)
setwd('/home/em_b/Desktop/AST_OPC_files')
matrix <- read.csv('salmon_out.csv')
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl=useDataset("mmusculus_gene_ensembl",
                   mart = ensembl)
filters = listFilters(ensembl)
attributes=listAttributes(ensembl)
#-----------------------------------------------------------------
dim(matrix)
geneid <- matrix$Gene
head(geneid)
listFilters(ensembl)
listAttributes(ensembl)
genes <-getBM(attributes = c('ensembl_transcript_id_version','external_gene_name'),
              filters = 'ensembl_transcript_id_version',
              values = geneid,
              mart = ensembl)
head(genes)
df <- matrix[matrix$Gene %in% unique(genes$ensembl_transcript_id_version),]
dim(df)
dim(matrix)
genes$ensembl_transcript_id_version
colnames(genes)[1]<- 'Gene'
colnames(df)[1]<-'Gene'
Ftable <- merge(genes, df, by="Gene")
Ftable<-Ftable[order(Ftable$external_gene_name),]
matrix<-Ftable
colnames(matrix)
matrix<-matrix %>%
  group_by(external_gene_name) %>%
  summarize(norm.ast_1 = sum(norm.ast_1),
            norm.ast_2 = sum(norm.ast_2),
            norm.ast_3 = sum(norm.ast_3),
            norm.ast_4 = sum(norm.ast_4),
            norm.opc_1 = sum(norm.opc_1),
            norm.opc_2 = sum(norm.opc_2),
            norm.opc_3 = sum(norm.opc_3),
            norm.opc_4 = sum(norm.opc_4),
            tumor.ast_1 = sum(tumor.ast_1),
            tumor.ast_2 = sum(tumor.ast_2),
            tumor.opc_1 = sum(tumor.opc_1),
            tumor.opc_2 = sum(tumor.opc_2))
omit<-matrix$external_gene_name[!grepl('Gm',matrix$external_gene_name)]
omit<-omit[!grepl('Rik',omit)]
matrix<-data.frame(matrix)
matrix<-subset(matrix,matrix$external_gene_name%in%omit)
matrix<-matrix[-1,]
entrez_id <-getBM(attributes = c('entrezgene_id','external_gene_name'),
                  filters = 'external_gene_name',
                  values = matrix$external_gene_name,
                  mart = ensembl)
final<-merge(entrez_id,matrix,by='external_gene_name')
break 
write.csv(final,file = 'matrix_mm39.csv')
#----if using older mm10 mouse assembly
ensembl<-useEnsembl(biomart = 'genes',
                    dataset = 'mmusculus_gene_ensembl',
                    version=102)
