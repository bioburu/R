library(DiffBind)
samples<-read.csv('/home/em_b/work_stuff/ChIP_RNA_comparison_templates/diffbind_test_samplesheet.csv')
samples
#dbObj<-dba(sampleSheet = samples)
#dbObj
#plot(dbObj)
test<-dba.analyze(samples)#,
                  #method = DBA_ALL_METHODS,
                  #bBlacklist = TRUE,
                  #bGreylist = TRUE,
                  #bRetrieveAnalysis = TRUE,
                  #bParallel = TRUE
                  #)
cat('ERROR-----Unable to complete analysis.
Warning messages:
1: No contrasts added. There must be at least two sample groups with at least three replicates. 
2: No contrasts added. There must be at least two sample groups with at least three replicates. ')
