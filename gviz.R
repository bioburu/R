library(Gviz)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm39)
library(ChIPseeker)
library(BRGenomics)
#----GenomeAxisTrack----------------------------------------------------------
gen<-'mm39'
chr<-'chr9'
bm <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                 dataset = "mmusculus_gene_ensembl")
biomTrack <- BiomartGeneRegionTrack(genome = "mm39",
                                    #chromosome = chr, 
                                    #start = 79282981,
                                    #end = 79286980,
                                    name = "ENSEMBL",
                                    biomart = bm,
                                    symbol = 'Gnb5')
plotTracks(biomTrack,
           transcriptAnnotation='symbol')
gtrack <- GenomeAxisTrack()
plotTracks(list(gtrack,biomTrack),
           transcriptAnnotation='symbol')
#------show chromosome ideogram
itrack <- IdeogramTrack(genome = gen,chromosome=chr)
plotTracks(list(gtrack,itrack,biomTrack),
           transcriptAnnotation='symbol')
cat('add bases to genomic tracks') 
strack <- SequenceTrack(Mmusculus, chromosome = chr)
plotTracks(list(gtrack,itrack,strack,biomTrack),
           transcriptAnnotation='symbol')

#-----enter in chip seq files here 
chip1<-read.delim('/home/em_b/work_stuff/chipseq2/thra_0hr_comb/peaks/thra_0hr_IgG_peaks.narrowPeak',
                 header = FALSE)
head(chip1)
colnames(chip1)<-c('chr','start','end','peak','score','NA','fold_change','-log10pval','-log10qval','summit')
head(chip1)
list<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19')
chip1<-subset(chip1, subset = chr %in% list)
table(chip1$chr)
chip1<-DataTrack(data = chip1$fold_change, start = chip1$start,
                end = chip1$end, chromosome = chr, genome = gen, 
                name = 'THRA_0hr_peaks')
#----------------------
chip2<-read.delim('/home/em_b/work_stuff/chipseq2/thra_b22_comb/peaks/thra_b22_IgG_peaks.narrowPeak',
                  header = FALSE)
head(chip2)
colnames(chip2)<-c('chr','start','end','peak','score','NA','fold_change','-log10pval','-log10qval','summit')
head(chip2)
list<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19')
chip2<-subset(chip2, subset = chr %in% list)
table(chip2$chr)
chip2<-DataTrack(data = chip2$fold_change, start = chip2$start,
                 end = chip2$end, chromosome = chr, genome = gen, 
                 name = 'THRA-B22_peaks')
#----------------------
chip3<-read.delim('/home/em_b/work_stuff/chipseq2/thra_t3_comb/peaks/thra_t3_IgG_peaks.narrowPeak',
                  header = FALSE)
head(chip3)
colnames(chip3)<-c('chr','start','end','peak','score','NA','fold_change','-log10pval','-log10qval','summit')
head(chip3)
list<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19')
chip3<-subset(chip3, subset = chr %in% list)
table(chip3$chr)
chip3<-DataTrack(data = chip3$fold_change, start = chip3$start,
                 end = chip3$end, chromosome = chr, genome = gen, 
                 name = 'THRA-T3_peaks')
plotTracks(list(gtrack,itrack,strack,biomTrack,
                chip1,chip2,chip3),
           chromosome = chr,
           transcriptAnnotation='symbol',
           extend.left = 0.5,
           extend.right = 100000,
           col=c('grey'),
           from = 75013570,
           to=75545923,
           type=c('histogram'),
           legend=TRUE)
gc()

#------enter in rnaseq files here 
rna1 <- DataTrack(range ='/home/em_b/work_stuff/chipseq2/rna/SRR23386662/SRR23386662.rmdup.sort.chr.bam',
                     genome = "mm39",
                     type = "horizon", 
                     name = 'Rnaseq_T3',
                     #window = -1, 
                     chromosome = 'chr9')
rna2 <- DataTrack(range ='/home/em_b/work_stuff/chipseq2/rna/SRR23386665/SRR23386665.rmdup.sort.chr.bam',
                  genome = "mm39",
                  type = "horizon", 
                  name = 'Rnaseq_pbs',
                  #window = -1, 
                  chromosome = 'chr9')
cat('options for types are: histogram,horizon,heatmap,polygon,mountain,b and you can layer')
plotTracks(list(gtrack,itrack,strack,biomTrack,
                chip1,chip2,chip3,rna1,rna2),
           chromosome = chr,
           transcriptAnnotation='symbol',
           #extend.left = 0.5,
           #extend.right = 100000,
           col=c('black'),
           from = 75013570,
           to=75545923,
           type=c('b','histogram'),
           legend=TRUE)


