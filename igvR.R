library(igvR)
igv<-igvR()
setBrowserWindowTitle(igv, "ChIP-seq")
setGenome(igv, 'mm39')
showGenomicRegion(igv, 'Neurod2')
for(i in 1:4) zoomOut(igv)
roi <- getGenomicRegion(igv)
gr.roi <- with(roi, GRanges(seqnames=chrom, ranges = IRanges(start, end)))
param <- ScanBamParam(which=gr.roi, what = scanBamWhat())
bamFile <- '/home/em_b/Desktop/chip_validations_final/bam/D1.0h.tra.rmdup.sort.bam'
alignments <- readGAlignments(bamFile, use.names=TRUE, param=param)
track <- GenomicAlignmentTrack(trackName='D1', alignments, visibilityWindow=10000000, trackHeight=200) 
displayTrack(igv, track)
bamFile <- '/home/em_b/Desktop/chip_validations_final/bam/A1.0hr.input.rmdup.sort.bam'
alignments <- readGAlignments(bamFile, use.names=TRUE, param=param)
track2 <- GenomicAlignmentTrack(trackName='input', alignments, visibilityWindow=10000000, trackHeight=200) 
displayTrack(igv,track2)
peak <- read.delim('/home/em_b/Desktop/chip_validations_final/in_house/peaks_q0.01/A1_D1_peaks.xls',
                     comment.char="#")
head(peak)
peak_df<-data.frame(chrom=peak$chr,
                   start=peak$start,
                   end=peak$end,
                   score=peak$pileup)
head(peak_df)
track <- DataFrameQuantitativeTrack("macs3_peaks", peak_df, color="red", autoscale=TRUE)
displayTrack(igv, track)

#----wig files
bw.slice<-import('/home/em_b/work_stuff/chip_thra/GSM940399_chapseq/ChAPSeq_TRa.wig',
                 which=gr.roi)
track<-GRangesQuantitativeTrack('test',bw.slice,autoscale=TRUE)
displayTrack(igv,track)
