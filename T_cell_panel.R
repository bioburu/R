#-------General T cell markers----------------------------------------------
VlnPlot(data, features = c('CD3D','CD3G','CD3E','CD4','CD8A','CD8B'))
#-------MAIT----------------------------------------------------------------------------
VlnPlot(data, features = c('CD44','KLRB1','IL18R1','CCR5','CXCR6','CCR6'))
VlnPlot(data, features = c())

#-------General  markers-------------------------------------------------
VlnPlot(data, features = c('PTPRC','CCR7','FOXP3','LCK','CD28','LAP3'))
VlnPlot(data, features = c('ZAP70','TRAC','TRBC1','TRBC2','PLCG1','ITK'))
VlnPlot(data, features = c('FYN','ADA','','CD69','CD2','CD40LG'))
VlnPlot(data, features = c('FAS','FASLG','VAV1','CD6','FLT3LG','CTLA4'))

#---------T cell receptor alpha variables-------------------------------------
VlnPlot(data, features = c('TRAV1-1','TRAV1-2','TRAV10','TRAV12-1','TRAV12-2','TRAV12-3'))
VlnPlot(data, features = c('TRAV13-1','TRAV13-2','TRAV14DV4','TRAV16','TRAV17','TRAV18'))
VlnPlot(data, features = c('TRAV19','TRAV2','TRAV20','TRAV21','TRAV22','TRAV23DV6'), pt.size = 1)
VlnPlot(data, features = c('TRAV24','TRAV25','TRAV26.1','TRAV26.2','TRAV27','TRAV29DV5'), pt.size = 1)
VlnPlot(data, features = c('TRAV3','TRAV30','TRAV34','TRAV36DV7','TRAV38.1','TRAV38.2DV8'), pt.size = 1)
VlnPlot(data, features = c('TRAV39','TRAV4','TRAV40','TRAV41','TRAV5','TRAV6'), pt.size = 1)
VlnPlot(data, features = c('TRAV8.1','TRAV8.2','TRAV8.3','TRAV8.4','TRAV8.5','TRAV9.2'), pt.size = 1)
#--------T cell receptor beta variables-------------------------------------------
VlnPlot(data, features = c('TRBV10.1','TRBV10.2','TRBV10.3','TRBV11.1','TRBV11.2','TRBV11.3'), pt.size = 1)
VlnPlot(data, features = c('TRBV12.3','TRBV12.4','TRBV12.5','TRBV13','TRBV14','TRBV15'), pt.size = 1)
#metabolic activity
VlnPlot(data, features = c("GAPDH","ACACA","IDH2","HK1","G6PD","PRDX2"))
VlnPlot(data, features = c("PRDX2","ATP1A1","ATP1B1","CPT1A"))
#activation/exhaustion markers
VlnPlot(data, features = c("JAK2","STAT3","IRAK1","IRAK3","MAPK1","POGLUT1"))
VlnPlot(data, features = c("LAG3","CCL22","CXCL4","SH2D1A","SLAMF7","SLAMF8"))
#effector markers
VlnPlot(data, features = c("GZMA","GZMB","PRF1","PSAP","LAMP2","LAG3"))
####cYTOKEINSSSSSS-------------------------
VlnPlot(data, features = c("IL2","IL2RA","IL3","IL4","IL4I1","IL4R"))
VlnPlot(data, features = c("IFNG","IL1RAP","IL6","IL7R","CXCL8","IL10"))
VlnPlot(data, features = c("IL10RA","IL10RB","IL10RB.AS1","IL11","IL12RB1","IL12RB2"))
VlnPlot(data, features = c("IL13","IL13RA1","IL15","IL15RA","IL16","IL17RA"))
VlnPlot(data, features = c("IL17RB","IL17RC","IL18","IL18BP","IL21R","IL23A",""))
VlnPlot(data, features = c("","IL27RA","IL32","TNF","TNFAIP1","TNFAIP3","TNFSF10"))
VlnPlot(data, features = c("IL31","CSF1","CSF2","LDHA","CD46","GIMAP1"))
####PRR and MHC-----------------------------------------------------------
RidgePlot(data, features = c("TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7",
                             "TLR8","TLR9","TLR10"), cols  = c('grey','red','blue'))
RidgePlot(data, features = c("HLA.A","HLA.B","HLA.C","HLA.G","HLA.H","HLA.J","HLA.L",
                             "HLA.DRA","HLA.DRB1",
                             "LILRB1","RFX5","CIITA"), cols  = c('grey','red','blue'))
