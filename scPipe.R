library(scPipe)
output_folder <- '/home/em_b/Desktop/FCCC/GSE218223_scATACseq_mouseCD8/my_results'
r1<-file.path('/home/em_b/Desktop/FCCC/GSE218223_scATACseq_mouseCD8/my_results/combine_R1.fastq.gz') 
r2<-file.path('/home/em_b/Desktop/FCCC/GSE218223_scATACseq_mouseCD8/my_results/combine_R2.fastq.gz') 
r3<-file.path('/home/em_b/Desktop/FCCC/GSE218223_scATACseq_mouseCD8/my_results/combine_R3.fastq.gz')
cat('valid_barcode_file is optional')
sc_atac_trim_barcode (r1            = r1, 
                      r2            = r2, 
                      bc_file       = r3, 
                      rmN           = TRUE,
                      rmlow         = TRUE,
                      output_folder = output_folder)
demux_r1<-file.path(output_folder, "demux_completematch_combine_R1.fastq.gz")
demux_r2<-file.path(output_folder, "demux_completematch_combine_R2.fastq.gz")
aligned_bam<-sc_aligning(ref = '/home/em_b/cellranger/scatacseq/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa', 
                R1 = demux_r1, 
                R2 = demux_r2, 
                nthreads  = 24,
                output_folder = output_folder,
                index_path = '/home/em_b/subread/mm10/genome_index',
                tech = 'atac',
                input_format = 'FASTQ',
                type = 'dna')
sorted_tagged_bam <- sc_atac_bam_tagging (inbam = aligned_bam, 
                       output_folder =  output_folder, 
                       bam_tags      = list(bc="CB", mb="OX"), 
                       nthreads      =  24)

break 
removed <- sc_atac_remove_duplicates(sorted_tagged_bam,output_folder = output_folder)
removed
sc_atac_create_fragments(inbam = sorted_tagged_bam,
                         output_folder = output_folder)
#sort -k1,1 -k2,2n -k3,3n fragments.bed >fragments.sorted.bed
cat('convert fragment.bed to .tsv.gz and .tsv.gz.tbi files and place into final_output')
cat('sort -k1,1 -k2,2n -k3,3n fragments.bed >fragments.sorted.bed')
cat('bgzip fragments.sorted.bed')
cat('tabix -p bed fragments.sorted.bed.gz')
cat('convert all bed files to tsv')
break 
cat('macs3 callpeak -f BAMPE -t tagged_sorted_markdup.bam -g mm -n combined -B -q 0.01 --call-summits')
#-------------------------------------------------------------------------------
features<-file.path(output_folder,'combined_peaks.narrowPeak')
fragments<-file.path(output_folder,'fragments.bed')
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
sce <- sc_atac_create_sce(input_folder = output_folder,
                   organism     = "mm10",
                   feature_type = "peak",
                   pheno_data   = NULL,
                   report       = FALSE)

