library(HiContacts)
mcool_file <- CoolFile('/home/dev/Desktop/HiCUP_test/fuckingwork.cool')
range <- 'chr1:20000-80000' # range of interest
availableResolutions(mcool_file)
hic <- HiCExperiment::import('/home/dev/Desktop/HiCUP_test/fuckingwork.cool',
                             format = 'cool',
                             focus=range)#,
                             #resolution = 1000)
hic
interaction_table<-data.frame(interactions(hic))
summary(interaction_table$count)
plotMatrix(hic)
#https://jserizay.com/HiContacts/
