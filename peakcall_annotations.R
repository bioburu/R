D3_A3_peaks.xls <- '/home/em_b/work_stuff/FCCC/chipseq/macs/D3/D3_A3_peaks.xls'
#------------read in peaks file
peaks.xls <- read.delim(D3_A3_peaks.xls,comment.char="#")
table(peaks.xls$chr)
#----------remove chromosomes not needed
list<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19')
peaks.xls<-subset(peaks.xls, subset = chr %in% list)
table(peaks.xls$chr)
reactome<-GRanges(peaks.xls)
View(data.frame(reactome))
cat('Reactome is for all peaks general output')
#-----This is the PEAK CALL ALIGNER
gene <- seq2gene(reactome,
                 tssRegion = c(-1000, 1000),
                 flankDistance = 3000,
                 TxDb=txdb)
############################################3
or
peakAnno <- annotatePeak(macsPeaks_GR, tssRegion=c(-1000, 1000), 
                         TxDb=txdb, 
                         annoDb="org.Mm.eg.db")
