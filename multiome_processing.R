#--download scRNAseq fastqs
fastq-dump SRR26312478.sra --split-files
#--cell ID and UMI data in Read 1 (~28 bp)
#--Read 2 is the RNA insert (~91 bp)
#--get 2 index files and 2 read files and rename and format for cellranger arc:
#--rnaseq_S1_L001_I1_001.fastq rnaseq_S1_L001_I2_001.fastq rnaseq_S1_L001_R1_001.fastq rnaseq_S1_L001_R2_001.fastq
#--place all files into folder names 'rna'
pigz -p 20 <all rna files>

#---download scATACseq fastqs
#-- after cellranger-atac mkfastq, there four fastq.gz files will be generated. 
#-- I1, R1, R2 and R3. I1 is the 8 bp sample barcode, R1 is the forward read. 
#-- R2 is the 16 bp 10x feature barcode and R3 is the reverse read. 
pigz -p 20 <all atac files>
