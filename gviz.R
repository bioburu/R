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
chip<-read.delim('/home/em_b/work_stuff/chipseq2/thra_0hr_comb/peaks/thra_0hr_IgG_peaks.narrowPeak',
                 header = FALSE)
head(chip)
colnames(chip)<-c('chr','start','end','peak','score','na','FC','pval','qval','summit')
head(chip)
list<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19')
chip<-subset(chip, subset = chr %in% list)
table(chip$chr)
chip<-DataTrack(data = chip$FC, start = chip$start,
                end = chip$end, chromosome = chr, genome = gen, 
                name = 'Chip-seq')
plotTracks(list(gtrack,itrack,strack,biomTrack,chip),
           chromosome = chr,
           transcriptAnnotation='symbol',
           extend.left = 0.5,
           extend.right = 100000,
           col=c('red'),
           from = 75013570,
           to=75545923,
           type=c('histogram'),
           legend=TRUE)
gc()
rna <- DataTrack(range = '/home/em_b/work_stuff/chipseq2/rna/SRR23386662/SRR23386662.sorted.bam',
                     genome = "mm39",
                     type = "h", 
                     name = 'RNAseq',
                     #window = -1, 
                     chromosome = 'chr9')
class(rna)
rna
cat('options for types are: histogram,horizon,heatmap,polygon,mountain,b and you can layer')

plotTracks(list(gtrack,itrack,strack,biomTrack,chip,rna),
           chromosome = chr,
           transcriptAnnotation='symbol',
           extend.left = 0.5,
           extend.right = 100000,
           col=c('skyblue'),
           from = 75013570,
           to=75545923,
           type=c('b'),
           legend=TRUE)

