library(igvR)
igv<-igvR()
setBrowserWindowTitle(igv, "ChIP-seq")
setGenome(igv, 'mm10')
#showGenomicRegion(igv, "chr3:128,079,020-128,331,275")
# or
showGenomicRegion(igv, 'NEUROD2')
for(i in 1:4) zoomOut(igv)
roi <- getGenomicRegion(igv)
gr.roi <- with(roi, GRanges(seqnames=chrom, ranges = IRanges(start, end)))
param <- ScanBamParam(which=gr.roi, what = scanBamWhat())
#---enter in bam file
bamFile <- '/home/em_b/work_stuff/FCCC/chipseq/comb/bam/E2.B22.Me3.rmdup.sort.bam'
alignments <- readGAlignments(bamFile, use.names=TRUE, param=param)
track <- GenomicAlignmentTrack(trackName='b22.me3', alignments, visibilityWindow=10000000, trackHeight=200) 
track
displayTrack(igv, track)
tbl.pk <- read.delim('/home/em_b/work_stuff/FCCC/chipseq/comb/macs/E2/E2_A2_peaks.narrowPeak',
                     header=FALSE)
dim(tbl.pk) 
colnames(tbl.pk)<-c('chrom','start','end','name','score')
head(tbl.pk)
tbl.pk<-tbl.pk[,-c(4,6:10)]
unlist(lapply(tbl.pk, class))
track <- DataFrameQuantitativeTrack("macs3_peaks", tbl.pk, color="red", autoscale=TRUE)
displayTrack(igv, track)
#------------------------------------------------
bamFile <- '/home/em_b/work_stuff/FCCC/chipseq/comb/bam/E3.T3.me3.rmdup.sort.bam'
alignments <- readGAlignments(bamFile, use.names=TRUE, param=param)
track <- GenomicAlignmentTrack(trackName='t3.me3', alignments, visibilityWindow=10000000, trackHeight=200) 
track
displayTrack(igv, track)
tbl.pk <- read.delim('/home/em_b/work_stuff/FCCC/chipseq/comb/macs/E3/E3_A3_peaks.narrowPeak',
                     header=FALSE)
dim(tbl.pk) 
colnames(tbl.pk)<-c('chrom','start','end','name','score')
head(tbl.pk)
tbl.pk<-tbl.pk[,-c(4,6:10)]
unlist(lapply(tbl.pk, class))
track <- DataFrameQuantitativeTrack("macs3_peaks", tbl.pk, color="red", autoscale=TRUE)
displayTrack(igv, track)
