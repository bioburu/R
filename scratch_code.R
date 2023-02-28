hlh3yR1 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr1_3.4yr.txt')
dim(hlh3yR1)
head(hlh3yR1)
hlh3yR1 <- hlh3yR1[,-c(2:4)] 
head(hlh3yR1)
#---------------------------------------------------------------
hlh3yR2 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr2_3.4yr.txt')
dim(hlh3yR2)
head(hlh3yR2)
hlh3yR2 <- hlh3yR2[,-c(2:4)] 
head(hlh3yR2)
data <- merge(hlh3yR1, hlh3yR2, by='Name')
colnames(data) <- c('ensembl_id','3.5y1','3.5y2')
head(data)
#-----------------------------------------------------------------
hlh3yR3 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr3_3.4yr.txt')
dim(hlh3yR3)
head(hlh3yR3)
hlh3yR3 <- hlh3yR3[,-c(2:4)] 
head(hlh3yR3)
colnames(hlh3yR3)[1] <- 'ensembl_id' 
head(data)
data <- merge(data, hlh3yR3, by='ensembl_id')
colnames(data) <- c('ensembl_id','3.5y1','3.5y2','3.5y3')
head(data)
#------------------------------------------------------------------
hlh3yR4 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr4_3.4yr.txt')
dim(hlh3yR4)
head(hlh3yR4)
hlh3yR4 <- hlh3yR4[,-c(2:4)] 
head(hlh3yR4)
colnames(hlh3yR4)[1] <- 'ensembl_id' 
head(data)
data <- merge(data, hlh3yR3, by='ensembl_id')
colnames(data) <- c('ensembl_id','3.5y1','3.5y2','3.5y3','3.5y4')
head(data)
#----------------------------------------------------------------
hlh4y1 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr1_4.3yr.txt')
dim(hlh4y1)
head(hlh4y1)
hlh4y1 <- hlh4y1[,-c(2:4)] 
head(hlh4y1)
colnames(hlh4y1)[1] <- 'ensembl_id' 
head(hlh4y1)
head(data)
data <- merge(data, hlh4y1, by='ensembl_id')
colnames(data) <- c('ensembl_id','3.5_y1','3.5_y2','3.5_y3','3.5_y4','4.3_y1')
head(data)
#---------------------------------------------------------------------
hlh4y2 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr2_4.3yr.txt')
dim(hlh4y2)
head(hlh4y2)
hlh4y2 <- hlh4y2[,-c(2:4)] 
head(hlh4y2)
colnames(hlh4y2)[1] <- 'ensembl_id' 
head(hlh4y2)
head(data)
data <- merge(data, hlh4y2, by='ensembl_id')
colnames(data) <- c('ensembl_id','3.5_y1','3.5_y2','3.5_y3','3.5_y4','4.3_y1',
                    '4.2_y2')
head(data)
#-------------------------------------------------------------------
hlh4y3 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr3_4.3yr.txt')
dim(hlh4y3)
head(hlh4y3)
hlh4y3 <- hlh4y3[,-c(2:4)] 
head(hlh4y3)
colnames(hlh4y3)[1] <- 'ensembl_id' 
head(hlh4y3)
head(data)
data <- merge(data, hlh4y3, by='ensembl_id')
colnames(data) <- c('ensembl_id','3.5_y1','3.5_y2','3.5_y3','3.5_y4','4.3_y1',
                    '4.3_y2','4.3_y3')
head(data)
#--------------------------------------------------------------------
hlh4y4 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr4_4.3yr.txt')
dim(hlh4y4)
head(hlh4y4)
hlh4y4 <- hlh4y4[,-c(2:4)] 
head(hlh4y4)
colnames(hlh4y4)[1] <- 'ensembl_id' 
head(hlh4y4)
head(data)
data <- merge(data, hlh4y4, by='ensembl_id')
colnames(data) <- c('ensembl_id','3.5_y1','3.5_y2','3.5_y3','3.5_y4','4.3_y1',
                    '4.3_y2','4.3_y3','4.3_y4')
head(data)
#----------------------------------------------------------------------
hlh12y1 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr1_12yr.txt')
dim(hlh12y1)
head(hlh12y1)
hlh12y1 <- hlh12y1[,-c(2:4)] 
head(hlh12y1)
colnames(hlh12y1)[1] <- 'ensembl_id' 
head(hlh12y1)
head(data)
data <- merge(data, hlh12y1, by='ensembl_id')
colnames(data) <- c('ensembl_id','3.5_y1','3.5_y2','3.5_y3','3.5_y4','4.3_y1',
                    '4.3_y2','4.3_y3','4.3_y4','12_y1')
head(data)
#----------------------------------------------------------------------
hlh12y2 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr2_12yr.txt')
dim(hlh12y2)
head(hlh12y2)
hlh12y2 <- hlh12y2[,-c(2:4)] 
head(hlh12y2)
colnames(hlh12y2)[1] <- 'ensembl_id' 
head(hlh12y2)
head(data)
data <- merge(data, hlh12y2, by='ensembl_id')
colnames(data) <- c('ensembl_id','3.5_y1','3.5_y2','3.5_y3','3.5_y4','4.3_y1',
                    '4.3_y2','4.3_y3','4.3_y4','12_y1','12_y2')
head(data)
#----------------------------------------------------------------------
hlh12y3 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr3_12yr.txt')
dim(hlh12y3)
head(hlh12y3)
hlh12y3 <- hlh12y3[,-c(2:4)] 
head(hlh12y3)
colnames(hlh12y3)[1] <- 'ensembl_id' 
head(hlh12y3)
head(data)
data <- merge(data, hlh12y3, by='ensembl_id')
colnames(data) <- c('ensembl_id','3.5_y1','3.5_y2','3.5_y3','3.5_y4','4.3_y1',
                    '4.3_y2','4.3_y3','4.3_y4','12_y1','12_y2','12_y3')
head(data)
#---------------------------------------------------------------------
norm1 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr1_norm.txt')
dim(norm1)
head(norm1)
norm1 <- norm1[,-c(2:4)] 
head(norm1)
colnames(norm1)[1] <- 'ensembl_id' 
head(norm1)
head(data)
data <- merge(data, norm1, by='ensembl_id')
colnames(data) <- c('ensembl_id','3.5_y1','3.5_y2','3.5_y3','3.5_y4','4.3_y1',
                    '4.3_y2','4.3_y3','4.3_y4','12_y1','12_y2','12_y3','norm_1')
head(data)
#---------------------------------------------------------------------------------
norm2 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr2_norm.txt')
dim(norm2)
head(norm2)
norm2 <- norm2[,-c(2:4)] 
head(norm2)
colnames(norm2)[1] <- 'ensembl_id' 
head(norm2)
head(data)
data <- merge(data, norm2, by='ensembl_id')
colnames(data) <- c('ensembl_id','3.5_y1','3.5_y2','3.5_y3','3.5_y4','4.3_y1',
                    '4.3_y2','4.3_y3','4.3_y4','12_y1','12_y2','12_y3','norm_1',
                    'norm_2')
head(data)
#-----------------------------------------------------------------------------------
norm3 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr3_norm.txt')
dim(norm3)
head(norm3)
norm3 <- norm3[,-c(2:4)] 
head(norm3)
colnames(norm3)[1] <- 'ensembl_id' 
head(norm3)
head(data)
data <- merge(data, norm3, by='ensembl_id')
colnames(data) <- c('ensembl_id','3.5_y1','3.5_y2','3.5_y3','3.5_y4','4.3_y1',
                    '4.3_y2','4.3_y3','4.3_y4','12_y1','12_y2','12_y3','norm_1',
                    'norm_2','norm_3')
head(data)
#-----------------------------------------------------------------------------------
norm4 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr4_norm.txt')
dim(norm4)
head(norm4)
norm4 <- norm4[,-c(2:4)] 
head(norm4)
colnames(norm4)[1] <- 'ensembl_id' 
head(norm4)
head(data)
data <- merge(data, norm4, by='ensembl_id')
colnames(data) <- c('ensembl_id','3.5_y1','3.5_y2','3.5_y3','3.5_y4','4.3_y1',
                    '4.3_y2','4.3_y3','4.3_y4','12_y1','12_y2','12_y3','norm_1',
                    'norm_2','norm_3','norm_4')
head(data)
#-----------------------------------------------------------------------------------
norm5 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr5_norm.txt')
dim(norm5)
head(norm5)
norm5 <- norm5[,-c(2:4)] 
head(norm5)
colnames(norm5)[1] <- 'ensembl_id' 
head(norm5)
head(data)
data <- merge(data, norm5, by='ensembl_id')
colnames(data) <- c('ensembl_id','3.5_y1','3.5_y2','3.5_y3','3.5_y4','4.3_y1',
                    '4.3_y2','4.3_y3','4.3_y4','12_y1','12_y2','12_y3','norm_1',
                    'norm_2','norm_3','norm_4','norm_5')
head(data)
#-----------------------------------------------------------------------------------
norm6 <- read.delim('/home/amp_prog/rstudio/HLH/famHLHr6_norm.txt')
dim(norm6)
head(norm6)
norm6 <- norm6[,-c(2:4)] 
head(norm6)
colnames(norm6)[1] <- 'ensembl_id' 
head(norm6)
head(data)
data <- merge(data, norm6, by='ensembl_id')
colnames(data) <- c('ensembl_id','3.5_y1','3.5_y2','3.5_y3','3.5_y4','4.3_y1',
                    '4.3_y2','4.3_y3','4.3_y4','12_y1','12_y2','12_y3','norm_1',
                    'norm_2','norm_3','norm_4','norm_5','norm_6')
head(data)
str(data)
break 
setwd('/home/amp_prog/rstudio/HLH')
write.csv(data, file = 'ensembl.csv')
