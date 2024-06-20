library(Gviz)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm39)
library(ChIPseeker)
library(BRGenomics)
#----Neurod2 and thra 
gen<-'mm39'
chr<-'chr11'
start<-98116241
end<-98759832
chip_yaxis<-c(0,50)
#bm <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
#                 dataset = "mmusculus_gene_ensembl")
biomTrack <- BiomartGeneRegionTrack(genome = "mm39",
                                    name = "ENSEMBL",
                                    biomart = bm,
                                    symbol = 'Thra',
                                    background.panel = "#FFFEDB",
                                    background.title = '#999883',
                                    stacking = 'squish')
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = gen,chromosome=chr)
strack <- SequenceTrack(Mmusculus, chromosome = chr)
bam1 <- DataTrack(range ='/home/em_b/work_stuff/chip_thra/bam/D1.0h.tra.rmdup.sort.bam',
                  genome = "mm39",
                  name = '0hr_thra',
                  #chromosome = chr,
                  ylim=chip_yaxis,
                  background.title = "#999883",
                  col='red')
bam1
bam2 <- DataTrack(range ='/home/em_b/work_stuff/chip_thra/bam/D3.t3.tra.rmdup.sort.bam',
                  genome = "mm39",
                  name = 'T3_thra',
                  #chromosome = chr,
                  ylim=chip_yaxis,
                  background.title = "#999883",
                  col='green')
bam2
bam3 <- DataTrack(range ='/home/em_b/work_stuff/chip_thra/bam/E1.0h.ME3.rmdup.sort.bam',
                  genome = "mm39",
                  name = '0hr_Me3',
                  #chromosome = chr,
                  ylim=chip_yaxis,
                  background.title = "#999883",
                  col='blue')
bam3
bam4 <- DataTrack(range ='/home/em_b/work_stuff/chip_thra/bam/E3_t3_me3.rmdup.sort.bam',
                  genome = "mm39",
                  name = 'T3_Me3',
                  #chromosome = chr,
                  ylim=chip_yaxis,
                  background.title = "#999883",
                  col='skyblue')
bam4
bam5 <- DataTrack(range ='/home/em_b/work_stuff/chip_thra/bam/J1.0h.IgG.sort.rmdup.bam',
                  genome = "mm39",
                  name = 'IgG_ctrl',
                  #chromosome = chr,
                  ylim=chip_yaxis,
                  background.title = "#999883",
                  col='grey')
bam5
plotTracks(list(itrack,gtrack,biomTrack,strack,bam1,bam2,bam3,bam4,bam5),
           transcriptAnnotation='symbol',
           #col=c('black'),
           from = start,
           to=end,
           type=c('h','','',''),
           legend=TRUE,
           cex=1,
           #col.histogram='red',
           reverseStrand = FALSE,
           pch=21,
           showBandId=FALSE,
           cex.bands=1)
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

