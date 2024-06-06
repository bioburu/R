library(biomaRt)
library(dplyr)
setwd('/home/em_b/work_stuff/rnaseq')
matrix <- read.csv('mm39_biomart.csv')
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
genes <-getBM(attributes = c('ensembl_transcript_id_version','external_gene_name','entrezgene_id'),
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
colnames(Ftable)[1] <- c('transcript_version')
colnames(Ftable)[2] <- c('Gene')
rownames(Ftable)<-make.names(Ftable$Gene,unique = TRUE)
Ftable<-cbind(rownames(Ftable),Ftable)
colnames(Ftable)[1] <- c('tracking_id')
head(Ftable)
matrix<-Ftable
matrix<-matrix[,-1]
colnames(matrix)
matrix<-matrix %>%
  group_by(entrezgene_id) %>%
  summarize(MB_0hr_1 = sum(MB_0hr_1),
            MB_0hr_2 = sum(MB_0hr_2),
            MB_0hr_3 = sum(MB_0hr_3),
            MB_6hr_1 = sum(MB_6hr_1),
            MB_6hr_2 = sum(MB_6hr_2),
            MB_6hr_3 = sum(MB_6hr_3),
            MB.PBS_1 = sum(MB.PBS_1),
            MB.PBS_2 = sum(MB.PBS_2),
            MB.PBS_3 = sum(MB.PBS_3),
            MB.T3_1 = sum(MB.T3_1),
            MB.T3_2 = sum(MB.T3_2),
            MB.T3_3 = sum(MB.T3_3),
            MB.ezh2i_1 = sum(MB.ezh2i_1),
            MB.ezh2i_2 = sum(MB.ezh2i_2),
            MB.ezh2i_3 = sum(MB.ezh2i_3),
            MB.gfp_1 = sum(MB.gfp_1),
            MB.gfp_2 = sum(MB.gfp_2),
            MB.gfp_3 = sum(MB.gfp_3),
            MB.neurod1_1 = sum(MB.neurod1_1),
            MB.neurod1_2 = sum(MB.neurod1_2),
            MB.neurod1_3 = sum(MB.neurod1_3),
            MB.sh.neurod1_1 = sum(MB.sh.neurod1_1),
            MB.sh.neurod1_2 = sum(MB.sh.neurod1_2),
            MB.sh.neurod1_3 = sum(MB.sh.neurod1_3),
            GNP.0h_1 = sum(GNP.0h_1),
            GNP.0h_2 = sum(GNP.0h_2),
            GNP.0h_3 = sum(GNP.0h_3),
            GNP.PBS_1 = sum(GNP.PBS_1),
            GNP.PBS_2 = sum(GNP.PBS_2),
            GNP.PBS_3 = sum(GNP.PBS_3),
            GNP.T3_1 = sum(GNP.T3_1),
            GNP.T3_2 = sum(GNP.T3_2),
            GNP.T3_3 = sum(GNP.T3_3),)
matrix<-na.omit(matrix)
row.names(matrix)<-matrix$entrezgene_id

break 
write.csv(matrix,file = 'mm39_entrez.matrix.csv')

