library(Gviz)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm39)
library(ChIPseeker)
library(BRGenomics)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ReactomePA)
library(plyranges)
library(BRGenomics)
library(ggplot2)
library(UpSetR)
library(ChIPpeakAnno)
library(pathview)
setwd('/home/em_b/work_stuff/chipseq/final')
df<-read.csv('thra_peaks.csv')
list<-read.csv('upreg_rna.csv',
                row.names = 1)
list
#write.csv(downreg,file='downreg_rna.csv')
list<-row.names(list)
rna_chip_up<-subset(df, subset = SYMBOL %in% list)
cat(rna_chip_up$SYMBOL)
upreg_reactome <- enrichPathway(rna_chip_up$geneId,
                                organism = 'mouse')
upreg_reactome
dotplot(upreg_reactome,
        x='Count',
        color= 'qvalue',
        size= NULL,
        split=NULL,
        font.size=15,
        label_format=50,
        title='30 total Thra peaks')
View(data.frame(upreg_reactome))
edox <- setReadable(upreg_reactome, 'org.Mm.eg.db', 'ENTREZID')
cnetplot(edox,
         #foldChange=gene)#,
         node_label='all',
         colorEdge=TRUE)
cnetplot(edox,
         #foldChange=gene,
         node_label='all',
         colorEdge=TRUE,
         circular=TRUE)
cat('Do go next')
go <- enrichGO(gene = rna_chip_up$geneId,OrgDb = org.Mm.eg.db,ont = "BP")
View(data.frame(go))

pathways<-c('GO:0007416','GO:0048167','GO:0035249','GO:0021953','GO:1905606',
            'GO:0099174','GO:0099643','GO:0097106','GO:0046928','GO:0021549',
            'GO:0045666','GO:0034329') 
go@result = go@result[go@result$ID %in% pathways,]
go@result


GO_plot <- pairwise_termsim(go)
View(data.frame(GO_plot))
emapplot(GO_plot,
         color='p.adjust',
         layout.params = list(layout = 'dh'),
         edge.params = list(min = 0.5),
         cex.params = list(line = 0.1),
         pie='equal',
         shadowtext=TRUE,
         repel=TRUE,
         cex_label_category=1,
         showCategory = 50)
goplot(GO_plot,
       showCategory = 17,
       color = "p.adjust",
       layout = "sugiyama",
       geom = "text")
break 
#----GenomeAxisTrack----------------------------------------------------------
#bm <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
#                 dataset = "mmusculus_gene_ensembl")

#--------------Rora 
gen<-'mm39'
chr<-'chr9'
#---far far away peaks
#start<-79279000
#end<-79800000
rna_yaxis<-c(0,100)
chip_yaxis<-c(0,100)

#--far distal peaks
#start<-79279000
#end<-79400000

#--closer peaks
#start<-79282000
#end<-79288000

#--tester
start<-68450627
end<-69298000

biomTrack <- BiomartGeneRegionTrack(genome = "mm39",
                                    #chromosome = chr, 
                                    #start = 79282981,
                                    #end = 79286980,
                                    name = "ENSEMBL",
                                    biomart = bm,
                                    symbol = 'Rora',
                                    background.panel = "#FFFEDB",
                                    background.title = '#999883',
                                    stacking = 'squish')
break 
#---------------Neurod2
gen<-'mm39'
chr<-'chr11'

rna_yaxis<-c(0,800)

#--close peaks
chip_yaxis<-c(0,50)


#--far peaks
start<-98203000
end<-98223471
chip_yaxis<-c(0,300)

biomTrack <- BiomartGeneRegionTrack(genome = "mm39",
                                    #chromosome = chr, 
                                    #start = 79282981,
                                    #end = 79286980,
                                    name = "ENSEMBL",
                                    biomart = bm,
                                    symbol = 'Neurod2',
                                    background.panel = "#FFFEDB",
                                    background.title = '#999883',
                                    stacking = 'squish')
#------------------------------------------------------------------------------
#-----------Neurod1
gen<-'mm39'
chr<-'chr2'
#---far far away peaks
start<-79279000
end<-79300000
rna_yaxis<-c(0,11500)
chip_yaxis<-c(0,30)

#--far distal peaks
#start<-78279000
#end<-80400000

#--closer peaks
#start<-79282000
#end<-79288000

#--tester
#start<-79382981
#end<-79386980
biomTrack <- BiomartGeneRegionTrack(genome = "mm39",
                                    #chromosome = chr, 
                                    #start = 79282981,
                                    #end = 79286980,
                                    name = "ENSEMBL",
                                    biomart = bm,
                                    symbol = 'Neurod1',
                                    background.panel = "#FFFEDB",
                                    background.title = '#999883',
                                    stacking = 'squish')

#------Rbfox1
gen<-'mm39'
chr<-'chr16'
#---far far away peaks
rna_yaxis<-c(0,500)
chip_yaxis<-c(0,40)
start<-5600000
end<-7260344

#--close peaks
rna_yaxis<-c(0,500)
chip_yaxis<-c(0,20)
start<-7190000
end<-7230344



biomTrack <- BiomartGeneRegionTrack(genome = "mm39",
                                    #chromosome = chr, 
                                    #start = 79282981,
                                    #end = 79286980,
                                    name = "ENSEMBL",
                                    biomart = bm,
                                    symbol = 'Rbfox1',
                                    background.panel = "#FFFEDB",
                                    background.title = '#999883',
                                    stacking = 'squish')
#----Adcy1 ----------------------------------------------------------
gen<-'mm39'
chr<-'chr11'
start<-7003433
end<-7139506
rna_yaxis<-c(0,850)
chip_yaxis<-c(0,200)

#-------------------------------------------------------------------------------
biomTrack <- BiomartGeneRegionTrack(genome = "mm39",
                                    #chromosome = chr, 
                                    #start = 79282981,
                                    #end = 79286980,
                                    name = "ENSEMBL",
                                    biomart = bm,
                                    symbol = 'Adcy1',
                                    background.panel = "#FFFEDB",
                                    background.title = '#999883')
#----Nrxn3
gen<-'mm39'
chr<-'chr12'
start<-89089454
end<-90301709
rna_yaxis<-c(0,80)
chip_yaxis<-c(0,50)

#bm <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
#                 dataset = "mmusculus_gene_ensembl")
#-------------------------------------------------------------------------------
biomTrack <- BiomartGeneRegionTrack(genome = "mm39",
                                    #chromosome = chr, 
                                    #start = 79282981,
                                    #end = 79286980,
                                    name = "ENSEMBL",
                                    biomart = bm,
                                    symbol = 'Nrxn3',
                                    background.panel = "#FFFEDB",
                                    background.title = '#999883',
                                    stacking = 'squish')
#----
gen<-'mm39'
chr<-'chr5'
start<-28651838
end<-28772099
rna_yaxis<-c(0,80)
chip_yaxis<-c(0,50)

#bm <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
#                 dataset = "mmusculus_gene_ensembl")
#-------------------------------------------------------------------------------
biomTrack <- BiomartGeneRegionTrack(genome = "mm39",
                                    #chromosome = chr, 
                                    #start = 79282981,
                                    #end = 79286980,
                                    name = "ENSEMBL",
                                    biomart = bm,
                                    symbol = 'Shh',
                                    background.panel = "#FFFEDB",
                                    background.title = '#999883',
                                    stacking = 'squish')

#--------------------------
#-------------------------
#------------------------

gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = gen,chromosome=chr)
strack <- SequenceTrack(Mmusculus, chromosome = chr)

rna1 <- DataTrack(range ='/home/em_b/work_stuff/chipseq/rna/replicates_combined/RNAseq_T3_3reps.bam',
                  genome = "mm39",
                  type = "horizon", 
                  name = 'RNAseq_T3',
                  #window = -1, 
                  #chromosome = chr,
                  ylim=rna_yaxis,
                  #col='red',
                  #background.panel = "white",
                  background.title = "#999883")
rna1
rna2 <- DataTrack(range ='/home/em_b/work_stuff/chipseq/rna/replicates_combined/RNAseq_PBS_3reps.bam',
                  genome = "mm39",
                  type = "horizon", 
                  name = 'RNAseq_PBS',
                  #window = -1, 
                  #chromosome = chr,
                  ylim=rna_yaxis,
                  #col='red',
                  #background.panel = "white",
                  background.title = "#999883")
rna2
chip1 <- DataTrack(range ='/home/em_b/work_stuff/chipseq/thra_0hr_comb/thra_0hr_comb.bam',
                   genome = "mm39",
                   type = "horizon", 
                   name = 'THRA_0hr',
                   #window = -1, 
                   #chromosome = chr,
                   ylim=chip_yaxis,
                   #col='red',
                   #background.panel = "white",
                   background.title = "#999883")
chip2 <- DataTrack(range ='/home/em_b/work_stuff/chipseq/thra_b22_comb/thra_b22_comb.bam',
                   genome = "mm39",
                   type = "horizon", 
                   name = 'THRA_b22',
                   #window = -1, 
                   #chromosome = chr,
                   ylim=chip_yaxis,
                   #col='red',
                   #background.panel = "white",
                   background.title = "#999883")
chip3 <- DataTrack(range ='/home/em_b/work_stuff/chipseq/thra_t3_comb/thra_t3_comb.bam',
                   genome = "mm39",
                   type = "horizon", 
                   name = 'THRA_b22',
                   #window = -1, 
                   #chromosome = chr,
                   ylim=chip_yaxis,
                   #col='red',
                   #background.panel = "white",
                   background.title = "#999883")
chip4 <- DataTrack(range ='/home/em_b/work_stuff/chipseq/h3k27me3/E1.0h.ME3.rmdup.sort.bam',
                   genome = "mm39",
                   type = "horizon", 
                   name = 'H3k27me33_0hr',
                   #window = -1, 
                   #chromosome = chr,
                   ylim=chip_yaxis,
                   #col='red',
                   #background.panel = "white",
                   background.title = "#999883")
chip5 <- DataTrack(range ='/home/em_b/work_stuff/chipseq/h3k27me3/E2.B22.Me3.rmdup.sort.bam',
                   genome = "mm39",
                   type = "horizon", 
                   name = 'H3k27me33_b22',
                   #window = -1, 
                   #chromosome = chr,
                   ylim=chip_yaxis,
                   #col='red',
                   #background.panel = "white",
                   background.title = "#999883")
chip6 <- DataTrack(range ='/home/em_b/work_stuff/chipseq/h3k27me3/E3_t3_me3.rmdup.sort.bam',
                   genome = "mm39",
                   type = "horizon", 
                   name = 'H3k27me3_T3',
                   #window = -1, 
                   #chromosome = chr,
                   ylim=chip_yaxis,
                   #col='red',
                   #background.panel = "white",
                   background.title = "#999883")
chip7 <- DataTrack(range ='/home/em_b/work_stuff/chipseq/h3k27me3/G1.0h.ME3.rmdup.sort.bam',
                   genome = "mm39",
                   type = "horizon", 
                   name = 'H3k27me33_0hr',
                   #window = -1, 
                   #chromosome = chr,
                   ylim=chip_yaxis,
                   #col='red',
                   #background.panel = "white",
                   background.title = "#999883")
chip8 <- DataTrack(range ='/home/em_b/work_stuff/chipseq/h3k27me3/I1.0h.ME3.rmdup.sort.bam',
                   genome = "mm39",
                   type = "horizon", 
                   name = 'H3k27me33_0hr',
                   #window = -1, 
                   #chromosome = chr,
                   ylim=chip_yaxis,
                   #col='red',
                   #background.panel = "white",
                   background.title = "#999883")
chip9 <- DataTrack(range ='/home/em_b/work_stuff/chipseq/h3k27me3/G2.B22.Me3.rmdup.sort.bam',
                   genome = "mm39",
                   type = "horizon", 
                   name = 'H3k27me33_b22',
                   #window = -1, 
                   #chromosome = chr,
                   ylim=chip_yaxis,
                   #col='red',
                   #background.panel = "white",
                   background.title = "#999883")
chip10 <- DataTrack(range ='/home/em_b/work_stuff/chipseq/h3k27me3/I2.B22.Me3.rmdup.sort.bam',
                   genome = "mm39",
                   type = "horizon", 
                   name = 'H3k27me33_b22',
                   #window = -1, 
                   #chromosome = chr,
                   ylim=chip_yaxis,
                   #col='red',
                   #background.panel = "white",
                   background.title = "#999883")
chip11 <- DataTrack(range ='/home/em_b/work_stuff/chipseq/h3k27me3/G3.T3.Me3.rmdup.sort.bam',
                    genome = "mm39",
                    type = "horizon", 
                    name = 'H3k27me3_T3',
                    #window = -1, 
                    #chromosome = chr,
                    ylim=chip_yaxis,
                    #col='red',
                    #background.panel = "white",
                    background.title = "#999883")
chip12 <- DataTrack(range ='/home/em_b/work_stuff/chipseq/h3k27me3/I3.T3.Me3.rmdup.sort.bam',
                    genome = "mm39",
                    type = "horizon", 
                    name = 'H3k27me3_T3',
                    #window = -1, 
                    #chromosome = chr,
                    ylim=chip_yaxis,
                    #col='red',
                    #background.panel = "white",
                    background.title = "#999883")
h3k27me3_0h_ot<-OverlayTrack(trackList=list(chip4,chip7,chip8),
                             background.title = "#999883")
h3k27me3_b22_ot<-OverlayTrack(trackList=list(chip5,chip9,chip10),
                              background.title = "#999883")
h3k27me3_t3_ot<-OverlayTrack(trackList=list(chip6,chip11,chip12),
                             background.title = "#999883")
#---tester
plotTracks(list(itrack,gtrack,biomTrack,strack,
                chip1,chip2,chip3,h3k27me3_0h_ot,h3k27me3_b22_ot,h3k27me3_t3_ot,
                rna1,rna2),
           transcriptAnnotation='symbol',
           col=c('black'),
           from = start,
           to=end,
           type=c('h','','gradient','g'),
           legend=TRUE,
           cex=1,
           col.histogram='red',
           reverseStrand = FALSE,
           pch=21,
           showBandId=FALSE,
           cex.bands=1)


