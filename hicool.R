library(ShortRead)
library(HiCool)
library(HiContacts)
library(HiCcompare)
library(HiCDOC)
r1<-file.path('/home/em_b/Desktop/hicool/SRR21038697/SRR21038697split/SRR21038697filtered_R1.part_001.fastq.gz')
x<-readFastq(r1)
x
head(x@sread)
r2<-file.path('/home/em_b/Desktop/hicool/SRR21038697/SRR21038697split/SRR21038697filtered_R2.part_001.fastq.gz')
x2<-readFastq(r2)
x2
head(x2@sread)
setwd('/home/em_b/Desktop/hicool/SRR21038697/SRR21038697split')
HiCool(r1,
       r2,
       restriction = 'DpnII',
       genome = 'mm10',
       threads = 22,
       exclude_chr = 'Mito|chrM|MT',
       output = './HiCool_output/',
       keep_bam = TRUE,
       build_report = TRUE,
       resolutions = c(1000,4000,8000,16000))
