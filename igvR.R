library(igvR)
igv<-igvR()
setBrowserWindowTitle(igv, "ChIP-seq")
setGenome(igv, 'mm39')
#showGenomicRegion(igv, "chr3:128,079,020-128,331,275")
# or
showGenomicRegion(igv, 'NEUROD1')
for(i in 1:4) zoomOut(igv)
roi <- getGenomicRegion(igv)
gr.roi <- with(roi, GRanges(seqnames=chrom, ranges = IRanges(start, end)))
param <- ScanBamParam(which=gr.roi, what = scanBamWhat())
#---enter in bam file
bamFile <- '/home/em_b/work_stuff/chipseq2/thra_0hr_comb/thra_0hr_comb.bam'
alignments <- readGAlignments(bamFile, use.names=TRUE, param=param)
track <- GenomicAlignmentTrack(trackName='thra_0hr', alignments, visibilityWindow=10000000, trackHeight=200) 
track
displayTrack(igv, track)
#----enter in peaks-------------------------
tbl.pk <- read.delim('/home/em_b/work_stuff/chipseq2/thra_0hr_comb/peaks/thra_0hr_3rep_peaks.narrowPeak',
                     header=FALSE)
dim(tbl.pk) 
colnames(tbl.pk)<-c('chrom','start','end','name','score')
head(tbl.pk)
tbl.pk<-tbl.pk[,-c(4,6:10)]
unlist(lapply(tbl.pk, class))
track <- DataFrameQuantitativeTrack("macs3_peaks", tbl.pk, color="red", autoscale=TRUE)
displayTrack(igv, track)
#-------enter in bam file 
bamFile <- '/home/em_b/work_stuff/chipseq2/thra_0hr_comb/thra_0hr_comb_input.bam'
alignments <- readGAlignments(bamFile, use.names=TRUE, param=param)
track <- GenomicAlignmentTrack(trackName='thra_0hr_input', alignments, visibilityWindow=10000000, trackHeight=200) 
track
displayTrack(igv, track)
break 
