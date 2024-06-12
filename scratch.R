library(igvR)
#BiocManager::install('Biostrings')
igv<-igvR()
setBrowserWindowTitle(igv, "ChIP-seq")
setGenome(igv, 'mm39')
break 
showGenomicRegion(igv, 'Thra')
for(i in 1:2) zoomOut(igv)
roi <- getGenomicRegion(igv)
gr.roi <- with(roi, GRanges(seqnames=chrom, ranges = IRanges(start, end)))
param <- ScanBamParam(which=gr.roi, what = scanBamWhat())
bamFile <- '/Users/burudpc/Desktop/bam/thra_0hr_comb.bam'
alignments <- readGAlignments(bamFile, use.names=TRUE, param=param)
track <- GenomicAlignmentTrack(trackName='crtl', alignments, visibilityWindow=10000000, trackHeight=200) 
displayTrack(igv, track)

bamFile <- '/Users/burudpc/Desktop/bam/thra_t3_comb.bam'
alignments <- readGAlignments(bamFile, use.names=TRUE, param=param)
track2 <- GenomicAlignmentTrack(trackName='T3', alignments, visibilityWindow=10000000, trackHeight=200) 
displayTrack(igv,track2)

bamFile <- '/Users/burudpc/Desktop/bam/thra_0hr_comb_input.bam'
alignments <- readGAlignments(bamFile, use.names=TRUE, param=param)
track3 <- GenomicAlignmentTrack(trackName='input', alignments, visibilityWindow=10000000, trackHeight=200) 
displayTrack(igv,track3)

break 
