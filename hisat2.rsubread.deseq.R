library(Rsubread)
library(biomaRt)
library(dplyr)
library(DESeq2)
counts<-featureCounts(files = '/home/em_b/work_stuff/chipseq2/rna/SRR23386662/SRR23386662.rmdup.sort.chr.bam',
                      annot.inbuilt = 'mm39',
                      isPairedEnd = TRUE)
T3_1<-data.frame(counts$counts)
T3_1<-cbind(row.names(T3_1),T3_1)
colnames(T3_1)<-c('Gene','MB.T3_1')
head(T3_1)
#-------------------------------------------------------------------------------
counts<-featureCounts(files = '/home/em_b/work_stuff/chipseq2/rna/SRR23386663/SRR23386663.rmdup.sort.chr.bam',
                      annot.inbuilt = 'mm39',
                      isPairedEnd = TRUE)
T3_2<-data.frame(counts$counts)
T3_2<-cbind(row.names(T3_2),T3_2)
colnames(T3_2)<-c('Gene','MB.T3_2')
head(T3_2)
#-------------------------------------------------------------------------------
counts<-featureCounts(files = '/home/em_b/work_stuff/chipseq2/rna/SRR23386664/SRR23386664.rmdup.sort.chr.bam',
                      annot.inbuilt = 'mm39',
                      isPairedEnd = TRUE)
T3_3<-data.frame(counts$counts)
T3_3<-cbind(row.names(T3_3),T3_3)
colnames(T3_3)<-c('Gene','MB.T3_3')
head(T3_3)
#-=-----------------------------------------------------------------------------
counts<-featureCounts(files = '/home/em_b/work_stuff/chipseq2/rna/SRR23386665/SRR23386665.rmdup.sort.chr.bam',
                      annot.inbuilt = 'mm39',
                      isPairedEnd = TRUE)
PBS_1<-data.frame(counts$counts)
PBS_1<-cbind(row.names(PBS_1),PBS_1)
colnames(PBS_1)<-c('Gene','MB.PBS_1')
head(PBS_1)
#-------------------------------------------------------------------------------
counts<-featureCounts(files = '/home/em_b/work_stuff/chipseq2/rna/SRR23386666/SRR23386666.rmdup.sort.chr.bam',
                      annot.inbuilt = 'mm39',
                      isPairedEnd = TRUE)
PBS_2<-data.frame(counts$counts)
PBS_2<-cbind(row.names(PBS_2),PBS_2)
colnames(PBS_2)<-c('Gene','MB.PBS_2')
head(PBS_2)
#-------------------------------------------------------------------------------
counts<-featureCounts(files = '/home/em_b/work_stuff/chipseq2/rna/SRR23386667/SRR23386667.rmdup.sort.chr.bam',
                      annot.inbuilt = 'mm39',
                      isPairedEnd = TRUE)
PBS_3<-data.frame(counts$counts)
PBS_3<-cbind(row.names(PBS_3),PBS_3)
colnames(PBS_3)<-c('Gene','MB.PBS_3')
head(PBS_3)
#-------------------------------------------------------------------------------
matrix<-cbind(PBS_1,
              PBS_2$MB.PBS_2,
              PBS_3$MB.PBS_3,
              T3_1$MB.T3_1,
              T3_2$MB.T3_2,
              T3_3$MB.T3_3)
colnames(matrix)<-c('Gene','PBS_1','PBS_2','PBS_3','T3_1','T3_2','T3_3')
#--------------------biomart
Sys.setenv("http_proxy"="http://my.proxy.org:9999")
listMarts()
ensembl=useMart('ENSEMBL_MART_ENSEMBL')
#listDatasets(ensembl)
ensembl=useDataset("mmusculus_gene_ensembl",
                   mart = ensembl)
geneid <- matrix$Gene
head(geneid)
#listFilters(ensembl)
#listAttributes(ensembl)
genes <-getBM(attributes = c('external_gene_name','entrezgene_id'),
              filters = 'entrezgene_id',
              values = geneid,
              mart = ensembl)
head(genes)
colnames(genes)[2]<-'Gene'
colnames(genes)
output <- merge(genes, matrix, by='Gene')
#---------------------------------------------------------
output<-data.frame(output %>%
  group_by(external_gene_name) %>%
  summarize(T3_1 = sum(T3_1),
            T3_2 = sum(T3_2),
            T3_3 = sum(T3_3),
            PBS_1 = sum(PBS_1),
            PBS_2 = sum(PBS_2),
            PBS_3 = sum(PBS_3)))
row.names(output)<-output$external_gene_name
output<-output[-1,-1]
setwd('/home/em_b/work_stuff/chipseq2/rna')
#write.csv(output,file='T3_rnaseq_df.csv')
#-------------------------------------------------------------------------------
names<-colnames(output)
names
condition<-c('B','B','B','A','A','A')
type<-c('paired','paired','paired','paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
row.names(coldata)<-make.names(coldata$names,
                               unique=TRUE)
coldata<-coldata[,-1]
coldata
#-----------------------------------------------------------------------------
colnames(output)
deseq<-DESeqDataSetFromMatrix(countData = output,colData = coldata,design = ~ condition)
deseq
DE<-DESeq(deseq)
plotMA(DE,ylim=c(-5,5))
plotDispEsts(DE)
#-------------------------------------------------------------------------------
results<-results(DE)
results
x<-data.frame(results)
x<-x[order(x$pvalue, decreasing=FALSE),]
x<-subset(x,padj< 0.05)
summary(x)
upreg<-subset(x,log2FoldChange> 1.5)
downreg<-subset(x,log2FoldChange< -1.5)
#write.csv(upreg,file='upreg_rna_df.csv')
#write.csv(downreg,file='downreg_rna_df.csv')
