setwd('/home/deviancedev/Downloads')
matrix<-read.delim('GSE244778_Read_Counts.txt.gz')
matrix$GENEID<-gsub("\\_.*","",matrix$GENEID)
