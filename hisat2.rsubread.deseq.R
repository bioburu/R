library(Rsubread)
library(biomaRt)
library(dplyr)
library(DESeq2)
library(EnhancedVolcano)
library(ggvenn)
library(gplots)
counts<-featureCounts(files = '/home/em_b/Desktop/dio2_il6/bam/GA1_norm_ast.rmdup.sort.chr1.19.bam',
                      annot.inbuilt = 'mm39',
                      isPairedEnd = TRUE)
norm_1<-data.frame(counts$counts)
norm_1<-cbind(row.names(norm_1),norm_1)
colnames(norm_1)<-c('Gene','norm_1')
head(norm_1)
#-------------------------------------------------------------------------------
counts<-featureCounts(files = '/home/em_b/Desktop/dio2_il6/bam/GA2_norm_ast.rmdup.sort.chr1.19.bam',
                      annot.inbuilt = 'mm39',
                      isPairedEnd = TRUE)
norm_2<-data.frame(counts$counts)
norm_2<-cbind(row.names(norm_2),norm_2)
colnames(norm_2)<-c('Gene','norm_2')
head(norm_2)
#-------------------------------------------------------------------------------
counts<-featureCounts(files = '/home/em_b/Desktop/dio2_il6/bam/GA3_norm_ast.rmdup.sort.chr1.19.bam',
                      annot.inbuilt = 'mm39',
                      isPairedEnd = TRUE)
norm_3<-data.frame(counts$counts)
norm_3<-cbind(row.names(norm_3),norm_3)
colnames(norm_3)<-c('Gene','norm_3')
head(norm_3)
#-=-----------------------------------------------------------------------------
counts<-featureCounts(files = '/home/em_b/Desktop/dio2_il6/bam/AS1_tumor_ast.rmdup.sort.chr1.19.bam',
                      annot.inbuilt = 'mm39',
                      isPairedEnd = TRUE)
tumor_1<-data.frame(counts$counts)
tumor_1<-cbind(row.names(tumor_1),tumor_1)
colnames(tumor_1)<-c('Gene','tumor_1')
head(tumor_1)
#-------------------------------------------------------------------------------
counts<-featureCounts(files = '/home/em_b/Desktop/dio2_il6/bam/AS4_tumor_ast.rmdup.sort.chr1.19.bam',
                      annot.inbuilt = 'mm39',
                      isPairedEnd = TRUE)
tumor_2<-data.frame(counts$counts)
tumor_2<-cbind(row.names(tumor_2),tumor_2)
colnames(tumor_2)<-c('Gene','tumor_2')
head(tumor_2)
#-------------------------------------------------------------------------------
counts<-featureCounts(files = '/home/em_b/Desktop/dio2_il6/bam/AS5_tumor_ast.rmdup.sort.chr1.19.bam',
                      annot.inbuilt = 'mm39',
                      isPairedEnd = TRUE)
tumor_3<-data.frame(counts$counts)
tumor_3<-cbind(row.names(tumor_3),tumor_3)
colnames(tumor_3)<-c('Gene','tumor_3')
head(tumor_3)
#-------------------------------------------------------------------------------
matrix<-cbind(norm_1,
              norm_2$norm_2,
              norm_3$norm_3,
              tumor_1$tumor_1,
              tumor_2$tumor_2,
              tumor_3$tumor_3)
colnames(matrix)<-c('Gene','norm.ast_1','norm.ast_2','norm.ast_3','tumor.ast_1','tumor.ast_2','tumor.ast_3')
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
colnames(output)
#---------------------------------------------------------
output<-data.frame(output %>%
                     group_by(external_gene_name) %>%
                     summarize(norm.ast_1 = sum(norm.ast_1),
                               norm.ast_2 = sum(norm.ast_2),
                               norm.ast_3 = sum(norm.ast_3),
                               tumor.ast_1 = sum(tumor.ast_1),
                               tumor.ast_2 = sum(tumor.ast_2),
                               tumor.ast_3 = sum(tumor.ast_3)))
row.names(output)<-output$external_gene_name
output<-output[-1,-1]
genes<-row.names(output)
genes
genes<-genes[!grepl('Gm',genes)]
genes<-genes[!grepl('Rik',genes)]
genes
output<-subset(output,row.names(output)%in%genes)
row.names(output)
setwd('/home/em_b/Desktop/dio2_il6')
#write.csv(output,file='tumor_norm_astrocytes.csv')
#-------------------------------------------------------------------------------
names<-colnames(output)
names
condition<-c('A','A','A','B','B','B')
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

EnhancedVolcano(results,
                lab = row.names(results),
                x='log2FoldChange',
                y='pvalue',
                title = 'Normal astrocytes vs Tumor-associated astrocytes',
                subtitle = 'Bulk RNAseq volcano plot. p=<0.05 log2FC=>1.5',
                legendLabels = NULL,
                legendIconSize = -1,
                legendPosition = 'bottom',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                shape = 5)
x<-data.frame(results)
x<-x[order(x$pvalue, decreasing=FALSE),]
x<-subset(x,padj< 0.05)
summary(x)
upreg<-subset(x,log2FoldChange> 1.5)
downreg<-subset(x,log2FoldChange< -1.5)
#-------------------------------------------------------------------------------
total_genes<-row.names(output)
upreg_genes<-row.names(upreg)
downreg_genes<-row.names(downreg)
venn <- list(upregulated = upreg_genes,
             downregulated = downreg_genes,
             total = total_genes)
venn
ggvenn(
  venn, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
v.table <- venn(venn)
print(v.table)
