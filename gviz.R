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
                                    #chromosome = 9, 
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
#-----DataTrack constructor is important
set.seed(255)
lim <- c(75213570, 75345923)
lim
lim[1]
lim[2]
cat('These are randomly sample coordinates')
coords <- sort(c(lim[1], 
                 sample(seq(from = lim[1], to = lim[2]), 99), 
                 lim[2]))
coords
dat <- runif(100, min = -10, max = 10)
dat
cat('look hard at DataTrack constructor')
coords[-length(coords)]
coords[-1]
dtrack <- DataTrack(data = dat, start = coords[-length(coords)],
                    end = coords[-1], chromosome = chr, genome = gen, 
                    name = "Uniform")
Test<-read.delim('/home/em_b/work_stuff/chipseq2/thra_0hr_comb/peaks/thra_0hr_IgG_peaks.narrowPeak',
                 header = FALSE)
head(Test)
colnames(Test)<-c('chr','start','end','peak','score','na','FC','pval','qval','summit')
head(Test)
list<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19')
Test<-subset(Test, subset = chr %in% list)
table(Test$chr)
#thra_0hr<-GRanges(Test)
#thra_0hr<-DataTrack(thra_0hr,
#                    name = 'thra_0hr')

test<-DataTrack(data = Test$FC, start = Test$start,
                end = Test$end, chromosome = chr, genome = gen, 
                name = "Test")
cat('options for types are: histogram,horizon,heatmap,polygon,mountain,b and you can layer')
plotTracks(list(gtrack,itrack,strack,biomTrack,dtrack,test),
           chromosome = chr,
           transcriptAnnotation='symbol',
           extend.left = 0.5,
           extend.right = 100000,
           col=c('red'),
           from = 75013570,
           to=75545923,
           type=c('histogram'),
           legend=TRUE)
summary(Test$FC)
break 
