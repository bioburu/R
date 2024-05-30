library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(org.Mm.eg.db)
#----------single files
files<-list(Thra_0hr_1='/home/em_b/work_stuff/chipseq/final/peaks_q_0.05/D1_J1_peaks.narrowPeak',
            Thra_0hr_2='/home/em_b/work_stuff/chipseq/final/peaks_q_0.05/F1_K1_peaks.narrowPeak',
            #Thra_0hr_3='/home/em_b/work_stuff/chipseq/final/peaks_q_0.05/H1_L1_peaks.narrowPeak',
            Thra_B22_1='/home/em_b/work_stuff/chipseq/final/peaks_q_0.05/D2_J2_peaks.narrowPeak',
            #Thra_B22_2='/home/em_b/work_stuff/chipseq/final/peaks_q_0.05/F2_K2_peaks.narrowPeak',
            Thra_B22_3='/home/em_b/work_stuff/chipseq/final/peaks_q_0.05/H2_L2_peaks.narrowPeak',
            Thra_T3_1='/home/em_b/work_stuff/chipseq/final/peaks_q_0.05/D3_J3_peaks.narrowPeak')#,
            #Thra_T3_2='/home/em_b/work_stuff/chipseq/final/peaks_q_0.05/F3_K3_peaks.narrowPeak',
            #Thra_T3_3='/home/em_b/work_stuff/chipseq/final/peaks_q_0.05/H3_L3_peaks.narrowPeak')
#-----------merged files 
files<-list(Thra_0hr='/home/em_b/work_stuff/chipseq/thra_0hr_comb/peaks/Thra_IgG_peaks.narrowPeak',
            Thra_B22='/home/em_b/work_stuff/chipseq/thra_b22_comb/peaks/Thra_b22_IgG_peaks.narrowPeak',
            Thra_T3='/home/em_b/work_stuff/chipseq/thra_t3_comb/peaks/Thra_t3_IgG_peaks.narrowPeak')
#peakAnnoList <- lapply(files, annotatePeak, TxDb=TxDb.Mmusculus.UCSC.mm39.knownGene,)
peakAnnoList <- lapply(files, annotatePeak, TxDb=TxDb.Mmusculus.UCSC.mm39.knownGene,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
peakAnnoList
plotAnnoBar(peakAnnoList)
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)
enrichPeakOverlap(queryPeak     = files[[1]],
                  targetPeak    = unlist(files[3]),
                  TxDb          = TxDb.Mmusculus.UCSC.mm39.knownGene,
                  pAdjustMethod = "BH",
                  nShuffle      = 50,
                  chainFile     = NULL,
                  verbose       = FALSE,
                  mc.cores = 3)
break 
#---------------enrichGO biological processes
test <- compareCluster(geneCluster   = genes,
                           fun           = 'enrichGO',
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH",
                           OrgDb=org.Mm.eg.db,
                       ont='BP')
#--subset
test@compareClusterResult$ID
dotplot(test, showCategory = 30, title = 'Biologial Process',size='Count')
#--------------enrichGO molecular function
test <- compareCluster(geneCluster   = genes,
                       fun           = 'enrichGO',
                       pvalueCutoff  = 0.05,
                       pAdjustMethod = "BH",
                       OrgDb=org.Mm.eg.db,
                       ont='MF')
dotplot(test, showCategory = 50, title = 'Molecular Function',size='Count')
#--------------enrichGO cellular component
test <- compareCluster(geneCluster   = genes,
                       fun           = 'enrichGO',
                       pvalueCutoff  = 0.05,
                       pAdjustMethod = "BH",
                       OrgDb=org.Mm.eg.db,
                       ont='CC')
dotplot(test, showCategory = 50, title = 'cellular stuff',size='Count')
