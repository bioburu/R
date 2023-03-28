#-------General T cell markers----------------------------------------------
VlnPlot(data, features = c('CD3D','CD3G','CD4','CD8A','TRAC','TRBC1','CD28','CD7','PTPRC'))
#-------Memory/transcription/stim Markers-------------------------------------------------
VlnPlot(data, features = c('CCR7','CD27','SELL','CD69','ITGAE','FOXP3','CD40LG','CD2','CD55'))
#-------Activation/exhaustion markers-------------------------------------------
VlnPlot(data, features = c('LAG3','CD69','FYN','LTB','FAS','FASLG','CTLA4','LSP1','HCST'))
#---------T cell receptor alpha variables-------------------------------------
VlnPlot(data, features = c('TRAV1-1','TRAV1-2','TRAV10','TRAV12-1','TRAV12-2','TRAV12-3','TRAV13-1','TRAV13-2','TRAV14DV4'))
VlnPlot(data, features = c('TRAV16','TRAV17','TRAV18','TRAV19','TRAV2','TRAV20','TRAV21','TRAV22','TRAV23DV6'))
VlnPlot(data, features = c('TRAV24','TRAV25','TRAV26-1','TRAV26-2','TRAV27','TRAV29DV5','TRAV3','TRAV30','TRAV34'))
VlnPlot(data, features = c('TRAV8-1','TRAV8-2','TRAV8-3','TRAV8-4','TRAV8-5','TRAV9-2'))
#--------T cell receptor beta variables-------------------------------------------
VlnPlot(data, features = c('TRBV10-1','TRBV10-2','TRBV10-3','TRBV11-1','TRBV11-2','TRBV11-3'))
VlnPlot(data, features = c('TRBV12-3','TRBV12-4','TRBV12-5','TRBV13','TRBV14','TRBV15'))
#metabolic activity
VlnPlot(data, features = c("GAPDH","ACACA","IDH2","HK1","G6PD","PRDX2","PRDX2","ATP1A1","ATP1B1","CPT1A"))
#activation markers
VlnPlot(data, features = c("JAK2","STAT3","IRAK1","IRAK3","MAPK1","POGLUT1","CCL22","CXCL4","SH2D1A"))
VlnPlot(data, features = c("SLAMF7","SLAMF8"))
#effector markers
VlnPlot(data, features = c("GZMA","GZMB","PRF1","PSAP","LAMP1","LAMP2","LAG3","GNLY"))
####cYTOKEINSSSSSS-------------------------
####cYTOKEINSSSSSS-------------------------
VlnPlot(data, features = c("IFNG","IFNAR1","IFNLR1","IFNAR2","IL1B","IL1RAP","IL1RN","IL2","IL2RA"))
VlnPlot(data, features = c("IL2RB","IL2RG","IL4I1","IL4R","IL6","IL6R","IL6ST","IL7","IL7R","CXCL8"))
VlnPlot(data, features = c("IL10","IL10RA","IL10RB","IL10RB-AS1","IL11","IL12RA","IL12RB1","IL12RB2","IL13"))
VlnPlot(data, features = c("IL13RA1","IL15","IL15RA","IL16","IL17RA","IL17RB","IL17RC","IL18","IL18BP"))
VlnPlot(data, features = c("IL21R","IL23A","IL24","IL27","IL27RA","IL32","TNF","TNFAIP1","TNFAIP3"))
VlnPlot(data, features =c("TNFSF10","CSF1","CSF1R","CSF2","CSFRA","CSFRB"))
#---------MHC--------------------------------------------------------------------
VlnPlot(data, features = c("HLA-A","HLA-B","HLA-C","HLA-G","HLA-H","HLA-J","HLA-L","HLA-DRA","HLA-DRB1"))
VlnPlot(data, features = c("LILRB1","RFX5","CIITA"))
#---damage associated molecular pattern (alarmins)-----------------------
VlnPlot(data, features = c('HMGB1','HMGB2','HMGB3','VIM','TP53','HSP90AA1','HSPB1','HSPA1A','HSF1')) 
VlnPlot(data, features = c('HSF2','HSPA4','HSPA5','HSP90AB1','HSPA1B','HSPD1','HSPD1','HSPH1','S100A1')) 
VlnPlot(data, features = c('S100A10','S100A11','S100A13','S100A2','S100A3','S100A4','S100A5','S100A6','S100A7A')) 
VlnPlot(data, features = c('S100A8','S100A9','S100B','S100G','S100P','S100PBP','S100Z','MRPL1','SAAL1')) 
VlnPlot(data, features = c('AGER','SYK','CARD6','CARD8','CARD9','CARD11','CARD14','LGALS1','SAP130')) 
VlnPlot(data, features = c('TREM1','TREML2','IFIH1','DDX58','DDX12P','LY96','PTAFR','SCARB1','SCARB2')) 
VlnPlot(data, features = c('MAVS','RNF135','SEC14L1','TRIM25','TFAM')) 
#-------------phagosome genes----------------------------------
VlnPlot(data, features = c('VAMP3',"VAMP5","VAMP8",'MR1','CORO1C','TLR8')) 
VlnPlot(data, features = c('SEC22B','FCAR','FCGRT','STX18','DCSTAMP','MARCO')) 
VlnPlot(data, features = c('CD36','RAB10','RAB5A','LAMP3','LAMP5','M6PR')) 
VlnPlot(data, features=c('PIKFYVE','TAP1','TAP2','SEC61A1','SEC61B','SEC61G'))
#----------------RIG1--------------------------------------------
VlnPlot(data, features = c('DDX58','IFIH1','IRF7','TRIM25','CCL8','TAP1')) 
#-----------------il6---------------------------
VlnPlot(data, features = c('IL6','CD55','LAMP3','VAMP8','MRC1')) 
VlnPlot(data, features = c('CD14','MR1'),pt.size = 2) 

