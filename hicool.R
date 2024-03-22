library(ShortRead)
library(HiCool)
library(HiContacts)
library(HiCcompare)
library(HiCDOC)
setwd('/home/em_b/Desktop/hicool/SRR21038697/550mil')
r1<-file.path('SRR21038697filtered_R1.part_001.fastq.gz')
#--files are too big to load into environment
#x<-readFastq(r1)
#x
#head(x@sread)
r2<-file.path('SRR21038697filtered_R2.part_001.fastq.gz')
#x2<-readFastq(r2)
#x2
#head(x2@sread)
HiCool(r1,
       r2,
       restriction = 'DpnII,HinfI',
       genome = 'mm10',
       threads = 22,
       exclude_chr = 'Mito|chrM|MT',
       output = './HiCool_output/',
       keep_bam = TRUE,
       build_report = TRUE,
       resolutions = c(1000,4000,8000,16000))
break 

try
python3 -m pip install --user pipx
python3 -m pipx ensurepath 
pipx install chromosight 
chromosight -h 

chromosight detect c241d99d5e2_7833^mapped-R64-1-1^LYRTWZ.mcool::resolutions/1000 results
