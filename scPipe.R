library(scPipe)
output_folder <- '/home/em_b/Desktop/SRR22324456/output_folder'
r1<-file.path('/home/em_b/Desktop/SRR22324456/SRR22324456_S1_L001_R1_001.fastq.gz') 
r2<-file.path('/home/em_b/Desktop/SRR22324456/SRR22324456_S1_L001_R2_001.fastq.gz') 
r3<-file.path('/home/em_b/Desktop/SRR22324456/SRR22324456_S1_L001_R3_001.fastq.gz')
sc_atac_trim_barcode (r1            = r1, 
                      r2            = r2, 
                      bc_file       = r3, 
                      rmN           = TRUE,
                      rmlow         = TRUE,
                      output_folder = output_folder)
demux_r1<-file.path(output_folder, "demux_completematch_SRR22324456_S1_L001_R1_001.fastq.gz")
demux_r2<-file.path(output_folder, "demux_completematch_SRR22324456_S1_L001_R2_001.fastq.gz")
aligned_bam<-sc_aligning(ref = '/home/em_b/cellranger/scatacseq/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa', 
                R1 = demux_r1, 
                R2 = demux_r2, 
                nthreads  = 20,
                output_folder = output_folder)
sorted_tagged_bam <- sc_atac_bam_tagging (inbam = aligned_bam, 
                       output_folder =  output_folder, 
                       bam_tags      = list(bc="CB", mb="OX"), 
                       nthreads      =  20)
removed <- sc_atac_remove_duplicates(sorted_tagged_bam,output_folder = output_folder)
if (!isFALSE(removed))
  sorted_tagged_bam <- removed
sc_atac_create_fragments(inbam = sorted_tagged_bam,
                         output_folder = output_folder)
cat('macs3 callpeak -f BAMPE -t demux_completematch_aligned_tagged_sorted_markdup.bam -g mm -n name_me -B -q 0.01 --call-summits')
features<-file.path('/home/em_b/Desktop/SRR22324456/output_folder/SRR22324456_peaks.narrowPeak')
fragments<-file.path('/home/em_b/Desktop/SRR22324456/output_folder/fragments.bed')
sc_atac_feature_counting (fragment_file = fragments,
                          feature_input = features, 
                          bam_tags      = list(bc="CB", mb="OX"), 
                          feature_type  = "peak",
                          organism      = "mm10",
                          cell_calling  = "filter",
                          min_uniq_frags = 0,
                          min_frac_peak = 0,
                          min_frac_promoter = 0,
                          yieldsize     = 1000000,
                          exclude_regions = TRUE,
                          output_folder = output_folder,
                          fix_chr       = "none"
                          )
feature_matrix <- readRDS(file.path("/home/em_b/Desktop/SRR22324456/output_folder/unfiltered_feature_matrix.rds"))
dplyr::glimpse(feature_matrix)
sparseM <- readRDS(file.path("/home/em_b/Desktop/SRR22324456/output_folder/sparse_matrix.rds"))
dplyr::glimpse(sparseM)
sce <- sc_atac_create_sce(input_folder = output_folder,
                   organism     = "mm10",
                   feature_type = "peak",
                   pheno_data   = NULL,
                   report       = FALSE)
data<-readRDS(file.path('/home/em_b/Desktop/SRR22324456/output_folder/scPipe_atac_SCEobject.rds'))
df<-data.frame(as.matrix(counts(data)))
#----make fragments.tsv file from fragments.bed
break 
