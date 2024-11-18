library(DESeq2)
library(csaw)
library(ChIPseeker)
library(ReactomePA)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(clusterProfiler)
library(enrichplot)
library(pathview)
setwd('/home/em_b/work_stuff/FCCC/chipseq')
txdb<-TxDb.Mmusculus.UCSC.mm39.refGene
#-----minq= minimum quality score
param <- readParam(minq=20)
#----------generate counts for bins in all files
bam.files<-c(THRA_0hr='/home/em_b/work_stuff/FCCC/chipseq/bam/D1.0h.tra.rmdup.sort.bam',
             THRA_0hr='/home/em_b/work_stuff/FCCC/chipseq/bam/F1.0h.Tra.rmdup.sort.bam',
             THRA_0hr='/home/em_b/work_stuff/FCCC/chipseq/bam/H1.0h.Tra.rmdup.sort.bam',
             THRA_t3='/home/em_b/work_stuff/FCCC/chipseq/bam/D3.t3.tra.rmdup.sort.bam',
             THRA_t3='/home/em_b/work_stuff/FCCC/chipseq/bam/F3.T3.Tra.rmdup.sort.bam',
             THRA_t3='/home/em_b/work_stuff/FCCC/chipseq/bam/H3.T3_tra.rmdup.sort.bam')
bam.files
#-----count read overlaps in windows. ext= average length of fragments. width= width of window. 
data <- windowCounts(bam.files, ext=500, width=50, param=param,bin = TRUE)
data
#----get larger bins for normalization 
binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)
binned
data
binned$totals
data$totals
#---------filter out uninformative regions ie low counts
#-- the higher the filter value, the more restrictive 
keep <- filterWindowsGlobal(data=data,
                            background=binned)$filter > 5
summary(keep)
dim(data)
data <- data[keep,]
dim(data)

head(assay(data))
rowRanges(data)
colData(data)
counts<-assays(data)$counts
ranges<-data.frame(data@rowRanges)
matrix<-cbind(ranges,counts)
head(matrix)
summary(matrix)
gc()
counts_in<-matrix
matrix<-matrix[,-c(1:5)]
names<-colnames(matrix)
names
condition<-c('A','A','A','B','B','B')
type<-c('paired','paired','paired','paired','paired','paired')
coldata<-data.frame(cbind(names,condition,type))
row.names(coldata)<-coldata$names
coldata<-coldata[,-1]
head(coldata)
#-----------------------------------------------------------------------------
deseq<-DESeqDataSetFromMatrix(countData = matrix,colData = coldata,design = ~ condition)
deseq
deseq<-estimateSizeFactors(deseq)
deseq<-estimateDispersions(deseq)
test <- DESeq(deseq, test="LRT", reduced= ~ 1)
test<-nbinomLRT(deseq,reduced = 1)
resLRT <- results(test)
resLRT<-data.frame(resLRT)
resLRT
df<-cbind(counts_in,resLRT)
df<-subset(df, pvalue<0.05)
list<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19')
DF<-subset(df, subset = seqnames %in% list)
summary(DF$seqnames)
DF<-GRanges(DF)
peakAnno <- annotatePeak(DF,
                         TxDb=txdb, 
                         annoDb="org.Mm.eg.db")
out<-data.frame(peakAnno)
out
#write.csv(out,file = 'initial_deg_chipseq.csv')
#-------------------------
gr_out<-GRanges(out)
gr_out
GO_result <- enrichGO(gene = gr_out$geneId,
                      #universe = allGeneIDs,
                      OrgDb = org.Mm.eg.db,
                      ont = "BP")
GO_result_df <- data.frame(GO_result)
GO_result_df[1:2, ]
GO_result_df<-GO_result_df[order(GO_result_df$Count, decreasing=TRUE),]
GO_result_plot <- pairwise_termsim(GO_result)
emapplot(GO_result_plot, showCategory = 40)
#--------------------------------------
reactome <- enrichPathway(as.vector(gr_out$geneId),
                          organism = 'mouse')
reactome
dotplot(reactome)
foldchanges = out$log2FoldChange
names(foldchanges) = out$geneId
head(foldchanges)
#-------------------------------------------------
#setwd('/home/deviancedev/Desktop/drive_jan2024/FCCC/exosome_exp/bulk_tissue/kegg')
pv_bulk_tissues <- pathview(gene.data = foldchanges, pathway.id = "04722",
                            species = "mmu", out.suffix = "mmu_test")
