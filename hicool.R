library(ShortRead)
library(HiCool)
library(HiContacts)
library(HiCcompare)
library(HiCDOC)
setwd('/home/em_b/Desktop/hicool/SRR22161720_muCD8_MboI')
r1<-file.path('SRR22161720_muCD8_MboI_filtered_R1.fastq.gz')
r2<-file.path('SRR22161720_muCD8_MboI_filtered_R2.fastq.gz')
hcf <- HiCool(
  r1 = r1, 
  r2 = r2, 
  restriction = 'MboI', 
  resolutions = c(1000,4000, 8000, 16000), 
  genome = 'mm10', 
  output = './HiCool/',
  keep_bam = TRUE,
  build_report = TRUE,
  threads = 24)

#-----------------------------------
#-----Arima are DpnII,HinfI.
