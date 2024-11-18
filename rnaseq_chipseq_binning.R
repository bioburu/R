library(DESeq2)
library(pheatmap)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(biomaRt)
library(RColorBrewer)
library(csaw)
library(ChIPseeker)
library(ReactomePA)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(enrichplot)
library(pathview)
library(igvR)
library(GenomicRanges)
library(enrichplot)
library(ReactomePA)
library(plyranges)
library(BRGenomics)
library(ggplot2)
library(UpSetR)
library(ChIPpeakAnno)
setwd('/home/T3_RNAseq')
matrix<-read.csv('dataframe2.csv')
colnames(matrix)
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
tester<-round(tester)
colnames(tester)
row.names(coldata)
deseq<-DESeqDataSetFromMatrix(countData = tester,colData = coldata,design = ~ condition)
deseq
DE<-DESeq(deseq)
plotMA(DE,ylim=c(-5,5))
plotDispEsts(DE)
results<-results(DE)
results
x<-data.frame(results)
x<-x[order(x$pvalue, decreasing=FALSE),]
x<-subset(x,pvalue< 0.05)
downreg<-subset(x,log2FoldChange< 0)
downregulated_rna<-cbind(row.names(downreg),downreg)
colnames(downregulated_rna)[1]<-'Gene'
txdb<-TxDb.Mmusculus.UCSC.mm39.refGene
bam.files<-c(H3K27me3_0hr.E1='/home/E1.0h.ME3.rmdup.sort.bam',
             H3K27me3_0hr.G1='/home/F1.0h.Tra.rmdup.sort.bam',
             H3K27me3_0hr.I1='/home/I1.0h.ME3.rmdup.sort.bam',
             H3K27me3_t3.E3='/home/E3_t3_me3.rmdup.sort.bam',
             H3K27me3_t3.G3='/home/G3.T3.Me3.rmdup.sort.bam',
             H3K27me3_t3.I3='/home/I3.T3.Me3.rmdup.sort.bam')
#---get paired end sizes in bam files ONLY IF NEEDED
bam1<-getPESizes('E1.0h.ME3.rmdup.sort.bam',
                 param=readParam(pe='both'))
bam2<-getPESizes('I1.0h.ME3.rmdup.sort.bam',
                 param=readParam(pe='both'))
bam3<-getPESizes('E3_t3_me3.rmdup.sort.bam',
                 param=readParam(pe='both'))
bam4<-getPESizes('G3.T3.Me3.rmdup.sort.bam',
                 param=readParam(pe='both'))
bam5<-getPESizes('I3.T3.Me3.rmdup.sort.bam',
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
gc()
cat('find all upregulated h3k27me3 peaks')
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
#write.csv(out,file='all_db_downregulations_b22.T3.csv')
cat('find all that match in RNA and ChIP')
length(downregulated_rna$Gene)
length(upregulated_chip$SYMBOL)
list<-downregulated_rna$Gene
subset_list<-subset(upregulated_chip, SYMBOL %in% list)
length(subset_list$SYMBOL)
setwd('/home/T3_RNAseq/')
#----------------ChIPseq heatmap  ------------------------------------------------
chip_htmp<-subset_list
row.names(chip_htmp)<-make.names(chip_htmp$SYMBOL,
                                 unique = TRUE)
chip_htmp<-chip_htmp[,-c(1:5,11:28)]
pheatmap(chip_htmp,
         scale = 'row',
         #labels_row = '',
         labels_col = c('H3K27me3_0hr','H3K27me3_0hr',
                        'H3K27me3_T3','H3K27me3_T3','H3K27me3_T3'),
         fontsize_row = 5,
         fontsize_col = 15,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'grey',
         fontsize = 15,
         color = colorRampPalette(rev(brewer.pal(n=7,name='RdYlBu')))(100),
         angle_col = 315,
         main = 'H3K27me3_0hr_vs_T3_ChIP_peaks')
#---------------- RNAseq heatmap  ------------------------------------
matrix<-read.csv('/home/T3_RNAseq/dataframe2.csv')
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
         labels_col = c('MB_0hr','MB_0hr','MB_0hr',
                        'MB_T3','MB_T3','MB_T3'),
         fontsize_row = 5,
         fontsize_col = 15,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'grey',
         fontsize = 15,
         color = colorRampPalette(rev(brewer.pal(n=7,name='RdYlBu')))(100),
         angle_col = 315,
         main = 'MB_tumors_0hr_vs_T3_RNAseq_transcripts')
setwd('/home/ChIP_RNA_comparison_templates')
#write.csv(subset_list,file='H3K27me3_T3_0hr_RNAChIP.csv')
rm(coldata,counts,counts_in,Data,deseq,df,DF,matrix,peakAnno,ranges,resLRT,subset,test)
gc()
#------DAVID pathway heatmaps RNAseq
negative_regulation_of_neuron_differentiation<-read.delim('/home/em_b/work_stuff/ChIP_RNA_comparison_templates/davidBI_GOTerm_negative regulation of neuron differentiation.txt')
list<-row.names(negative_regulation_of_neuron_differentiation)
pheatmap(subset[list,],
         scale = 'row',
         labels_col = c('MB_0hr','MB_0hr','MB_0hr',
                        'MB_T3','MB_T3','MB_T3'),
         fontsize_row = 18,
         fontsize_col = 18,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'grey',
         fontsize = 10,
         color = colorRampPalette(rev(brewer.pal(n=7,name='RdYlBu')))(100),
         angle_col = 315,
         main = 'GO Term: Negative_regulation_of_neuron_differentiation')
negative_regulation_of_transcription_from_RNA_polymerase_II_promoter<-read.delim('/home/em_b/work_stuff/ChIP_RNA_comparison_templates/davidBI_GOTerm_negative_regulation_of_transcription_from_RNA_polymerase_II_promoter.txt')
list<-row.names(negative_regulation_of_transcription_from_RNA_polymerase_II_promoter)
pheatmap(subset[list,],
         scale = 'row',
         labels_col = c('MB_0hr','MB_0hr','MB_0hr',
                        'MB_T3','MB_T3','MB_T3'),
         fontsize_row = 18,
         fontsize_col = 18,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'grey',
         fontsize = 10,
         color = colorRampPalette(rev(brewer.pal(n=7,name='RdYlBu')))(100),
         angle_col = 315,
         main = 'GO Term:Negative_regulation_of_transcription_from_RNA_polyII_promoter')
chromatin_binding<-read.delim('/home/em_b/work_stuff/ChIP_RNA_comparison_templates/davidBI_GOTerm_chromatin_binding.txt')
list<-row.names(chromatin_binding)
pheatmap(subset[list,],
         scale = 'row',
         labels_col = c('MB_0hr','MB_0hr','MB_0hr',
                        'MB_T3','MB_T3','MB_T3'),
         fontsize_row = 18,
         fontsize_col = 18,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'grey',
         fontsize = 15,
         color = colorRampPalette(rev(brewer.pal(n=7,name='RdYlBu')))(100),
         angle_col = 315,
         main = 'GO Term: Chromatin_binding')
#------DAVID pathway heatmaps RNAseq
chip<-read.csv('/home/H3K27me3_T3_0hr_RNAChIP.csv',
               row.names = 1)
row.names(chip)<-make.names(chip$SYMBOL,unique = TRUE)
chip<-chip[,-c(1:5,11:28)]
negative_regulation_of_neuron_differentiation<-read.delim('/home/davidBI_GOTerm_negative regulation of neuron differentiation.txt')
list<-row.names(negative_regulation_of_neuron_differentiation)
list
pheatmap(chip[list,],
         scale = 'row',
         labels_col = c('MB_0hr','MB_0hr','MB_0hr',
                        'MB_T3','MB_T3','MB_T3'),
         fontsize_row = 18,
         fontsize_col = 18,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'grey',
         fontsize = 10,
         color = colorRampPalette(rev(brewer.pal(n=7,name='RdYlBu')))(100),
         angle_col = 315,
         main = 'GO Term: Negative_regulation_of_neuron_differentiation')
negative_regulation_of_transcription_from_RNA_polymerase_II_promoter<-read.delim('/home/davidBI_GOTerm_negative_regulation_of_transcription_from_RNA_polymerase_II_promoter.txt')
list<-row.names(negative_regulation_of_transcription_from_RNA_polymerase_II_promoter)
pheatmap(chip[list,],
         scale = 'row',
         labels_col = c('MB_0hr','MB_0hr','MB_0hr',
                        'MB_T3','MB_T3','MB_T3'),
         fontsize_row = 18,
         fontsize_col = 18,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'grey',
         fontsize = 8,
         color = colorRampPalette(rev(brewer.pal(n=7,name='RdYlBu')))(100),
         angle_col = 315,
         main = 'GO Term:Negative_regulation_of_transcription_from_RNA_polyII_promoter')
chromatin_binding<-read.delim('/home/davidBI_GOTerm_chromatin_binding.txt')
list<-row.names(chromatin_binding)
pheatmap(chip[list,],
         scale = 'row',
         labels_col = c('MB_0hr','MB_0hr','MB_0hr',
                        'MB_T3','MB_T3','MB_T3'),
         fontsize_row = 18,
         fontsize_col = 18,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         border_color = 'grey',
         fontsize = 15,
         color = colorRampPalette(rev(brewer.pal(n=7,name='RdYlBu')))(100),
         angle_col = 315,
         main = 'GO Term: Chromatin_binding')
#------------------------------------------------------------------------------
igv<-igvR()
setBrowserWindowTitle(igv, "ChIP-seq")
setGenome(igv, 'mm39')
showGenomicRegion(igv, 'ISL1')
#showGenomicRegion(igv, 'EPHA2')
#for(i in 1:4) zoomOut(igv)
roi <- getGenomicRegion(igv)
gr.roi <- with(roi, GRanges(seqnames=chrom, ranges = IRanges(start, end)))
param <- ScanBamParam(which=gr.roi, what = scanBamWhat())
#---enter in bam file
bamFile <- '/home/I1.0h.ME3.rmdup.sort.bam'
alignments <- readGAlignments(bamFile, use.names=TRUE, param=param)
track <- GenomicAlignmentTrack(trackName='h3k27me3_0hr', alignments, visibilityWindow=10000000, trackHeight=200) 
track
displayTrack(igv, track)
bamFile <- '/home/I3.T3.Me3.rmdup.sort.bam'
alignments <- readGAlignments(bamFile, use.names=TRUE, param=param)
track <- GenomicAlignmentTrack(trackName='h3k27me3_T3', alignments, visibilityWindow=10000000, trackHeight=200) 
track
displayTrack(igv, track)
tbl.pk <- read.delim('/home/I3_C3_peaks.narrowPeak',
                     header=FALSE)
colnames(tbl.pk)<-c('chrom','start','end','name','score')
head(tbl.pk)
tbl.pk<-tbl.pk[,-c(4,6:10)]
unlist(lapply(tbl.pk, class))
track <- DataFrameQuantitativeTrack("h3k27me3_T3_peaks", tbl.pk, color="red", autoscale=TRUE)
displayTrack(igv, track)
