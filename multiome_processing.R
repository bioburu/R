#--download scRNAseq fastqs
fastq-dump SRR26312478.sra --split-files
#--cell ID and UMI data in Read 1 (~28 bp)
#--Read 2 is the RNA insert (~91 bp)
#--get 2 index files and 2 read files and rename and format for cellranger arc:
#--rnaseq_S1_L001_I1_001.fastq rnaseq_S1_L001_I2_001.fastq rnaseq_S1_L001_R1_001.fastq rnaseq_S1_L001_R2_001.fastq
