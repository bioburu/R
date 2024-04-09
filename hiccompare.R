library(HiCcompare)
library(EnsDb.Mmusculus.v79)
library(ChIPpeakAnno)
library(tibble)
setwd('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output')
bl6<-read.delim('bl6.chr7.103564001.104036000.tsv',sep='')
bl6<-as.matrix(bl6)
str(bl6)
balbc<-read.delim('balbc.chr7.103564001.104036000.tsv',sep='')
balbc<-as.matrix(balbc)
str(balbc)
gc()
hic_table<-create.hic.table(bl6,balbc,chr = 'chr7')
head(hic_table)
#--joint normalization 
hic_table<-hic_loess(hic_table,Plot = TRUE,Plot.smooth = FALSE)
#---filter data
filter_params(hic_table)
#---set filter and compare
hic_table <- hic_compare(hic_table, A.min=NA, adjust.dist = TRUE, p.method = 'fdr', Plot = TRUE)
hic_table<-hic_table[order(hic_table$p.adj,decreasing=FALSE),]
df<-subset(hic_table,p.adj< 0.1)
IntSet <- make_InteractionSet(df)
IntSet
#--plots
MD.plot1(M = hic_table$M,
         D = hic_table$D,
         mc = hic_table$mc,
         smooth = TRUE)
MD.plot2(M = hic_table$adj.M,
         D = hic_table$D,
         hic_table$p.value,
         smooth = FALSE)
break 
colnames(df)[1:3]<-c('seqnames','start','end')
Loops1 <- toGRanges(df)
Loops1
annoData<-toGRanges(EnsDb.Mmusculus.v79)
anno<-annotatePeakInBatch(Loops1,AnnotationData = annoData)
anno <- addGeneIDs(anno, orgAnn="org.Mm.eg.db", 
                   feature_id_type="ensembl_gene_id",
                   IDs2Add=c("symbol"))
anno
Loops1<-data.frame(Loops1)
Loops1<-add_column(Loops1, frag1_genes=anno$symbol, .after = 'end')
Loops1
df<-Loops1
break
colnames(df)[c(1:3,7:9)]<-c('chr1','start1','end1','seqnames','start','end')
Loops1 <- toGRanges(df)
Loops1
anno<-annotatePeakInBatch(Loops1,AnnotationData = annoData)
anno <- addGeneIDs(anno, orgAnn="org.Mm.eg.db", 
                   feature_id_type="ensembl_gene_id",
                   IDs2Add=c("symbol"))
anno
Loops1<-data.frame(Loops1)
Loops1<-add_column(Loops1, frag2_genes=anno$symbol, .after = 'end')
colnames(Loops1)[1:3]<-c('chr2','start2','end2')
df<-Loops1[,-c(5:6,18)]
cat('we need more replicates')
