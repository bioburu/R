#https://jserizay.com/OHCA/docs/devel/pages/topological-features.html
#library(ShortRead)
library(HiCool)
library(HiContacts)
library(GenomicRanges)
library(HiCExperiment)
library(BiocParallel)
#library(HiCcompare)
#library(HiCDOC)
library(WGCNA)
library(ggplot2)
library(patchwork)
setwd('/Users/burudpc/Desktop')
cf<-CoolFile('twofiles^mapped-mm10^CCP6TK.mcool')
pairs_file<- PairsFile('/Users/burudpc/Desktop/twofiles^mapped-mm10^CCP6TK.pairs')
availableResolutions(cf)
hic<-import(cf,
            resolution=16000,
            pairs=pairs_file,
            format = 'mcool')
hic
#---annotate A/B compartments
phasing_track<-BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
compts<-getCompartments(hic,chromosomes = 'chr1',genome = phasing_track)
compts
compartdf<-data.frame(topologicalFeatures(compts,'compartments'))
plotSaddle(compts)
metadata(compts)$eigens
compts<-autocorrelate(compts)
