# https://bioconductor.org/books/devel/OHCA/pages/visualization.html
library(ggplot2)
library(GenomicRanges)
library(InteractionSet)
library(HiCExperiment)
library(HiContactsData)
library(HiContacts)
library(rtracklayer)
r1 <- HiContactsData(sample = 'yeast_wt', format = 'fastq_R1')
View(r1)
r2 <- HiContactsData(sample = 'yeast_wt', format = 'fastq_R2')
r1
r2
library(HiCool)
HiCool(
  r1, 
  r2, 
  restriction = 'DpnII,HinfI', 
  resolutions = c(4000, 8000, 16000), 
  genome = 'R64-1-1', 
  output = './HiCool/'
)
hic<-import('/Users/burudpc/Desktop/HiCool/matrices/test.mcool',
            format = 'cool',
            resolution=4000,
            focus='XII')
hic
plotMatrix(hic)
plotMatrix(hic,maxDistance=2000000)
getLoops(hic)

hic@topologicalFeatures$loops

