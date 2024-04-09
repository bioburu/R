library(HiCool)
library(HiContacts)
library(GenomicRanges)
library(HiCExperiment)
library(BiocParallel)
library(HiCcompare)
library(patchwork)
library(GenomicRanges)
library(Matrix)
resolution<-8000
setwd('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output')
#-------------------------------------------------------------------------------
mcool_path<-file.path('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/matrices/SRR26855486_BL6.mcool')
pairs_path<-file.path('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/pairs/SRR26855486_BL6.pairs')
log_path<-file.path('/home/em_b/Desktop/hicool/GSE248029_il4_macro_DpnII_HinfI/mm10_output/logs/SRR26855486_BL6.log')
bl6_cf<-CoolFile(mcool_path,
                 pairs_path,
                 metadata = list(log=log_path),
                 resolution=resolution)
bl6_cf
bl6_hic<-import(bl6_cf)
bl6_hic
gc()
interactions(bl6_hic)
df<-data.frame(interactions(bl6_hic))
summary(df$seqnames1)
df <- df[df$seqnames1 %in% c('chr1','chr2','chr3','chr4','chr5',
                                 'chr6','chr7','chr8','chr9','chr10',
                                 'chr11','chr12','chr13','chr14','chr15',
                                 'chr16','chr17','chr18','chr19'),] 
summary(df$seqnames1)
gc()
df<-df[,-c(1,3:10,12:20,22)]
colnames(df)<-c('region1','region2','IF')
df<-as.matrix(df)
str(df)
gc()
sparseMtx<-Matrix(df, sparse = TRUE)
head(sparseMtx)
break 
writeMM(sparseMtx,file='bl6_spMtx.mtx')
#write.csv(df,file = 'interaction_counts.csv')
