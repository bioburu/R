library(ShortRead)
library(HiCool)
library(HiContacts)
library(HiCcompare)
library(HiCDOC)
setwd('/home/em_b/Desktop/hicool/SRR22161720.musC8.didnt work')
r1<-file.path('SRR22161720_1.fastq')
#--files are too big to load into environment
#x<-readFastq(r1)
#x
#head(x@sread)
r2<-file.path('SRR22161720_2.fastq')
#x2<-readFastq(r2)
#x2
#head(x2@sread)
#HiCool(r1,
#       r2,
#       restriction = 'DpnII',
#       genome = 'mm10',
#       threads = 24,
#       exclude_chr = 'Mito|chrM|MT',
#       output = './HiCool_output/',
#       keep_bam = TRUE,
#       build_report = TRUE)
HiCool(
  r1 = r1, 
  r2 = r2, 
  restriction = 'MboI', 
  #resolutions = c(1000,4000,8000,16000), 
  genome = 'mm10', 
  output = './HiCool/',
  threads = 24,
  keep_bam = TRUE,
  build_report = TRUE)
break 
