library(Rfastp)
setwd('/home/deviancedev/Desktop/drive_nov2023/FCCC/RNAseq_exosomes/RNAseq_exosomes_reads1')
pe_json_report <- rfastp(read1 = 'A_Cereb1_1.fq.gz',
                         read2 = 'A_Cereb1_2.fq.gz',
                         outputFastq = paste0('A_Cereb1_filtered'),
                         adapterTrimming = TRUE, qualityFiltering = TRUE,
                         qualityFilterPhred = 30)
summary <- qcSummary(pe_json_report)
summary
cat('Base quality plot')
plot <- curvePlot(pe_json_report)
plot
cat('GC content plot')
plot2 <- curvePlot(pe_json_report,curves = 'content_curves')
plot2
trim <- trimSummary(pe_json_report)
trim

