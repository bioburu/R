#--------------TME CD14+ genes--------------------------------------------
#-------------Cytoprotective autophagy phenotype---------------------- 
VlnPlot(data, features = c('VAMP8','CD59','HLA-A','HMGB1','HSPB1','CD68','DDX58','IFIH1','MILR1'),pt.size = 2)
VlnPlot(data,features = c('CD55','CD59','CD68','NR1H3','IL7','HLA-A'), pt.size = 2)
VlnPlot(data, features = c('HLA-B','HLA-C','MKI67','NFKB1','IL18'),pt.size = 2)
#-------------monocyte markers---------------------------------------------------
VlnPlot(data, features = c("CD14","FCGR3A","ITGAM","ITGAX","CD38","CD4",
                           "CD86","FCGR1A","CD68","CD163","MILR1","CD80"), pt.size = 2)
#-------Memory/transcription/stim Markers-------------------------------------------------
VlnPlot(data, features = c('CCR7','CD27','SELL','CD69','ITGAE','FOXP3','CD40LG','CD2','CD55'),pt.size = 2)
#-------Activation/exhaustion markers-------------------------------------------
VlnPlot(data, features = c('LAG3','CD69','FYN','LTB','FAS','FASLG','CTLA4','LSP1','HCST'),pt.size = 2)
VlnPlot(data, features = c("JAK2","STAT3","IRAK1","IRAK3","MAPK1","POGLUT1","CCL22","CXCL4","SH2D1A"),pt.size = 2)
VlnPlot(data, features = c("SLAMF7","SLAMF8","CD44",'TIGIT','PDCD1','CD38','MKI67','TGFB1'),pt.size = 2)
#-----metabolic activity--------------------------------------------------------------------------------------------
VlnPlot(data, features = c("GAPDH","ACACA","IDH2","HK1","G6PD","PRDX2",
                           "ATP1A1","ATP1B1","CPT1A"), pt.size = 2)
#activation/exhaustion markers
VlnPlot(data, features = c("JAK2","STAT3","IRAK1","IRAK3","MAPK1",
                           "POGLUT1","LAG3","CCL22","CCL8","SH2D1A",
                           "SLAMF7","SLAMF8"), pt.size = 2)
#effector markers
VlnPlot(data, features = c("GZMA","GZMB","PRF1","PSAP","LAMP1","LAMP2",
                           "ATG7","MRC1","MSR1","NR1H3"), pt.size = 2)
#metabolic activity
VlnPlot(data, features = c("GAPDH","ACACA","IDH2","HK1","G6PD","PRDX2","PRDX2","ATP1A1","ATP1B1","CPT1A"),pt.size = 2)
#effector markers
VlnPlot(data, features = c("GZMA","GZMB","PRF1","PSAP","LAMP1","LAMP2","LAG3","GNLY"),pt.size = 2)
####cYTOKEINSSSSSS-------------------------
####cYTOKEINSSSSSS-------------------------
VlnPlot(data, features = c("IFNG","IFNAR1","IFNLR1","IFNAR2","IL1B","IL1RAP","IL1RN","IL2","IL2RA"),pt.size = 2)
VlnPlot(data, features = c("IL2RB","IL2RG","IL4I1","IL4R","IL6","IL6R","IL6ST","IL7","IL7R","CXCL8"),pt.size = 2)
VlnPlot(data, features = c("IL10","IL10RA","IL10RB","IL10RB-AS1","IL11","IL12RA","IL12RB1","IL12RB2","IL13"),pt.size = 2)
VlnPlot(data, features = c("IL13RA1","IL15","IL15RA","IL16","IL17RA","IL17RB","IL17RC","IL18","IL18BP"),pt.size = 2)
VlnPlot(data, features = c("IL21R","IL23A","IL24","IL27","IL27RA","IL32","TNF","TNFAIP1","TNFAIP3"),pt.size = 2)
VlnPlot(data, features =c("TNFSF10","CSF1","CSF1R","CSF2","CSFRA","CSFRB"),pt.size = 2)
#---------MHC--------------------------------------------------------------------
VlnPlot(data, features = c("HLA-A","HLA-B","HLA-C","HLA-G","HLA-H","HLA-J","HLA-L","HLA-DRA","HLA-DRB1"),pt.size = 2)
VlnPlot(data, features = c("HLA.A","HLA.B","HLA.C","HLA.G","HLA.H","HLA.J","HLA.L","HLA.DRA","HLA.DRB1"),pt.size = 2)

VlnPlot(data, features = c("LILRB1","RFX5","CIITA"),pt.size = 2)
