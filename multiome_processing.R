#--download scRNAseq fastqs
fastq-dump SRR26312478.sra --split-files
#--cell ID and UMI data in Read 1 (~28 bp)
#--Read 2 is the RNA insert (~91 bp)
