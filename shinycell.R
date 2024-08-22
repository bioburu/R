#----make shiny app with object of choice. 
library(Seurat)
library(ShinyCell)
setwd('/home/deviancedev01/Desktop')
data<-readRDS('/home/deviancedev01/work_stuff/multiome_MB_tumorigenesis/rna_p7.rds')
scConf<-createConfig(data)
makeShinyApp(data,
             scConf,
             gene.mapping = TRUE,
             shiny.title = 'GSE243609:Medulloblastoma_tumorigenesis:P7:Ptch1-Atoh1EGFPmice')
#-----deploy app on shinyapps.io
library(rsconnect)
rsconnect::setAccountInfo(name='9a17md-chakkapong-burudpakdee',
                          token='EB963E0629DB74B4F40C8CE2621114B3',
                          secret='deDqnsWWrRZfMUby+x+4ppKNwS2BJ0Ns+IigjjLV')
rsconnect::deployApp('/home/deviancedev01/Desktop/shinyApp')
