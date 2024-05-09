cat('find all downregulated transcripts to subset with')
library(DESeq2)
library(pheatmap)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(dplyr)
library(RColorBrewer)
setwd('/home/em_b/work_stuff/FCCC/T3_RNAseq')
matrix<-read.csv('dataframe2.csv')
colnames(matrix)
#---if subsetting genes
#list<-out$SYMBOL
#list
#matrix<-subset(matrix, subset = Gene %in% list)
row.names(matrix)<-make.names(matrix$Gene,unique = TRUE)
matrix<-matrix[,-1]
head(matrix)
tester<-matrix[,-c(4:6,10:27)]
names<-colnames(tester)
names
condition<-c('A','A','A','B','B','B')
type<-c('paired','paired','paired','paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
row.names(coldata)<-make.names(coldata$names,
                               unique=TRUE)
coldata<-coldata[,-1]
coldata
#-----------------------------------------------------------------------------
tester<-round(tester)
colnames(tester)
row.names(coldata)
deseq<-DESeqDataSetFromMatrix(countData = tester,colData = coldata,design = ~ condition)
deseq
DE<-DESeq(deseq)
plotMA(DE,ylim=c(-5,5))
plotDispEsts(DE)
#-------------------------------------------------------------------------------
results<-results(DE)
results
x<-data.frame(results)
x<-x[order(x$pvalue, decreasing=FALSE),]
x<-subset(x,pvalue< 0.05)
downreg<-subset(x,log2FoldChange< 0)
downregulated_rna<-cbind(row.names(downreg),downreg)
colnames(downregulated_rna)[1]<-'Gene'
library(DESeq2)
library(csaw)
library(ChIPseeker)
library(ReactomePA)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(pheatmap)
setwd('/home/em_b/work_stuff/FCCC/chipseq')
txdb<-TxDb.Mmusculus.UCSC.mm39.refGene
bam.files<-c(H3K27me3_0hr.E1='/home/em_b/work_stuff/FCCC/chipseq/bam/E1.0h.ME3.rmdup.sort.bam',
             #H3K27me3_0hr.G1='/home/em_b/work_stuff/FCCC/chipseq/bam/F1.0h.Tra.rmdup.sort.bam',
             H3K27me3_0hr.I1='/home/em_b/work_stuff/FCCC/chipseq/bam/I1.0h.ME3.rmdup.sort.bam',
             H3K27me3_t3.E3='/home/em_b/work_stuff/FCCC/chipseq/bam/E3_t3_me3.rmdup.sort.bam',
             H3K27me3_t3.G3='/home/em_b/work_stuff/FCCC/chipseq/bam/G3.T3.Me3.rmdup.sort.bam',
             H3K27me3_t3.I3='/home/em_b/work_stuff/FCCC/chipseq/bam/I3.T3.Me3.rmdup.sort.bam')
#---get paired end sizes in bam files 
bam1<-getPESizes('/home/em_b/work_stuff/FCCC/chipseq/bam/E1.0h.ME3.rmdup.sort.bam',
                 param=readParam(pe='both'))
bam2<-getPESizes('/home/em_b/work_stuff/FCCC/chipseq/bam/I1.0h.ME3.rmdup.sort.bam',
                 param=readParam(pe='both'))
bam3<-getPESizes('/home/em_b/work_stuff/FCCC/chipseq/bam/E3_t3_me3.rmdup.sort.bam',
                 param=readParam(pe='both'))
bam4<-getPESizes('/home/em_b/work_stuff/FCCC/chipseq/bam/G3.T3.Me3.rmdup.sort.bam',
                 param=readParam(pe='both'))
bam5<-getPESizes('/home/em_b/work_stuff/FCCC/chipseq/bam/I3.T3.Me3.rmdup.sort.bam',
                 param=readParam(pe='both'))
summary(bam1$sizes)
summary(bam2$sizes)
summary(bam3$sizes)
summary(bam4$sizes)
summary(bam5$sizes)
param <- readParam(minq=30,
                   max.frag = 300,
                   pe='both')
data <- windowCounts(bam.files,
                     ext=200,
                     width=10,
                     param=param,
                     bin = TRUE)
#---choose windows to disregard
binned <- windowCounts(bam.files,
                     #ext=200,
                     width=10000,
                     param=param,
                     bin = TRUE)
#binned <- windowCounts(bam.files,
#                       bin=TRUE,
#                       width=10000,
#                       param=param)
break 
cat('keep data and binned if starting over')
gc()
cat('find all upregulated me3 peaks')
#-- the higher the filter value, the more restrictive 
log2(5)
keep <- filterWindowsGlobal(data=data,
                            background=binned)$filter > 0.99
Data <- data[keep,]
head(assay(Data))
Data<-normFactors(binned,se.out = Data)
colData(Data)
counts<-assays(Data)$counts
ranges<-data.frame(Data@rowRanges)
matrix<-cbind(ranges,counts)
gc()
counts_in<-matrix
matrix<-matrix[,-c(1:5)]
names<-colnames(matrix)
names
condition<-c('A','A','B','B','B')
type<-c('paired','paired','paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
row.names(coldata)<-coldata$names
coldata<-coldata[,-1]
head(coldata)
#-----------------------------------------------------------------------------
deseq<-DESeqDataSetFromMatrix(countData = matrix,colData = coldata,design = ~ condition)
deseq<-estimateSizeFactors(deseq)
deseq<-estimateDispersions(deseq)
test <- DESeq(deseq, test="LRT", reduced= ~ 1)
resLRT <- results(test)
resLRT<-data.frame(resLRT)
df<-cbind(counts_in,resLRT)
df<-subset(df, pvalue<0.05)
df<-subset(df,log2FoldChange>0)
list<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19')
DF<-subset(df, subset = seqnames %in% list)
summary(DF$seqnames)
DF<-GRanges(DF)
peakAnno <- annotatePeak(DF,
                         TxDb=txdb, 
                         annoDb="org.Mm.eg.db")
out<-data.frame(peakAnno)
summary(out$log2FoldChange)
upregulated_chip<-subset(out,
            log2FoldChange> 0)
summary(upregulated_chip$log2FoldChange)
setwd('/home/em_b/Downloads')
#write.csv(out,file='all_db_downregulations_b22.T3.csv')
break 
cat('find all that match in RNA and ChIP')
length(downregulated_rna$Gene)
length(upregulated_chip$SYMBOL)
list<-downregulated_rna$Gene
subset_list<-subset(upregulated_chip, SYMBOL %in% list)
length(subset_list$SYMBOL)
library(pheatmap)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(dplyr)
library(RColorBrewer)
setwd('/home/em_b/work_stuff/FCCC/T3_RNAseq/')
#----------------all conditions ------------------------------------------------
#---------------- neuronal differentiation  ------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
matrix<-read.csv('/home/em_b/work_stuff/FCCC/T3_RNAseq/dataframe2.csv')
list<-subset_list$SYMBOL
list
colnames(matrix)
matrix<-subset(matrix, subset = Gene %in% list)
row.names(matrix)<-make.names(matrix$Gene,unique = TRUE)
matrix<-matrix[,-1]
head(matrix)
subset<-matrix[,-c(4:6,10:21)]
colnames(subset)
pheatmap(subset,
         scale = 'row',
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'grey',
         fontsize = 7,
         color = colorRampPalette(rev(brewer.pal(n=7,name='RdYlBu')))(100))
setwd('/home/em_b/Downloads')
#write.csv(subset_list,file='T3_0hr_downRNA.upChIP.csv')
rm(coldata,counts,counts_in,Data,deseq,df,DF,matrix,peakAnno,ranges,resLRT,subset,test)
gc()
cat('gotcha bitch')
