ERROR: dependencies ‘impute’, ‘shinyFiles’, ‘openxlsx’, ‘GSA’ are not available for package ‘samr’
* removing ‘/home/deviancedev/R/x86_64-pc-linux-gnu-library/4.1/samr’

BiocManager::install('impute')
library(impute)
install.packages('shinyFiles')
library(shinyFiles)
install.packages('openxlsx')
library(openxlsx)
install.packages('GSA')
library('GSA')
#---now install samr from .tar.gz
