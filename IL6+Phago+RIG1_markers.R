#-------------phagosome genes----------------------------------
VlnPlot(data, features = c('VAMP3',"VAMP5","VAMP8",'MR1','CORO1C','TLR8'),pt.size = 2) 
VlnPlot(data, features = c('SEC22B','FCAR','FCGRT','STX18','DCSTAMP','MARCO'),pt.size = 2) 
VlnPlot(data, features = c('CD36','RAB10','RAB5A','LAMP3','LAMP5','M6PR'),pt.size = 2) 
VlnPlot(data, features=c('PIKFYVE','TAP1','TAP2','SEC61A1','SEC61B','SEC61G'),pt.size = 2)
#----------------RIG1--------------------------------------------
VlnPlot(data, features = c('DDX58','IFIH1','IRF7','TRIM25','CCL8','TAP1'),pt.size = 2) 
#-----------------il6---------------------------
VlnPlot(data, features = c('IL6','CD55','LAMP3','VAMP8','MRC1'),pt.size = 2) 
VlnPlot(data, features = c('CD14','MR1'),pt.size = 2) 
