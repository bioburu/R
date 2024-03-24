library(ShortRead)
library(HiCool)
library(HiContacts)
library(HiCcompare)
library(HiCDOC)
setwd('/home/em_b/Desktop/hicool/SRR22161720')
r1<-file.path('SRR22161720_1.part_001.fastq.gz')
r2<-file.path('SRR22161720_2.part_001.fastq.gz')
hcf <- HiCool(
  r1 = r1, 
  r2 = r2, 
  restriction = 'MboI', 
  resolutions = c(1000,4000, 8000, 16000), 
  genome = 'mm10', 
  output = './HiCool/',
  keep_bam = TRUE,
  build_report = TRUE)
#-----------------------------------
#-----Arima are DpnII,HinfI.
