#--example 1
biopsy<-read.csv('/home/deviancedev01/Downloads/WA_Fn-UseC_-Telco-Customer-Churn.csv',
                 stringsAsFactors = TRUE)%>%mutate_if(is.numeric, round)%>%mutate_if(is.numeric,as.integer)
biopsy$customerID <- NULL
biopsy <- na.omit(biopsy)
