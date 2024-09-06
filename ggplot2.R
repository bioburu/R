#------------simple
library(ggplot2)
library(reshape2)
library(graphics)
df<-read.csv('/home/deviancedev01/Desktop/interesting_exps/ggplot_worldecon/Table-4.-Govt-revenue.csv')
df$year<-sub(":.*", "", df$year)
df$year<-as.integer(df$year)
str(df)
df<-na.omit(df)
#melt data frame into long format
df <- melt(df ,  id.vars = 'year', variable.name = 'series')

#create line plot for each column in data frame
ggplot(df, aes(year, value)) +
  geom_line(aes(colour = series))+geom_text(
    label=df$series,
    nudge_x=0.45, nudge_y=0.1,
    check_overlap=T,
    size=2)

#-----complex
#---kaggle dataset
library(ggplot2)
library(reshape2)
library(graphics)
library(tidyverse) 
library(dplyr)
library(MASS)
library(lime)
library(corrplot)

df<-read.csv('/home/deviancedev01/Desktop/interesting_exps/ggplot_worldecon/Global Dataset of Inflation.csv',
             check.names = FALSE)
column.names<-df$Country
df<-na.omit(df)
df<-df %>% 
  group_by(Country) %>% 
  reframe(across(where(is.character)),
            across(where(is.numeric), sum, na.rm = T))
df<-data.frame(unique(df),check.names = FALSE)
row.names(df)<-df$Country
df<-df[,-c(1:2)]
df<-data.frame(t(df))
str(df)
df<-log10(df)
df<-cbind(row.names(df),df)
colnames(df)[1]<-'Years'
str(df)
df <- melt(df ,  id.vars = 'Years', variable.name = 'series')
ggplot(df, aes(Years, value)) +
  geom_line(aes(colour = series))+geom_text(
    label=df$series,
    nudge_x=0.45, nudge_y=0.1,
    check_overlap=T,
    size=3)+ theme(legend.position="none")+theme(axis.text.x=element_text(size=8))+labs(y= "Inflation_Rates", x = "Years")

usa<-subset(df,df$series=='United.States')
ggplot(usa, aes(Years, value)) +
  geom_line(aes(colour = series))+geom_text(
    label=usa$series,
    nudge_x=0.45, nudge_y=0.1,
    check_overlap=T,
    size=3)+ theme(legend.position="none")+geom_point()+theme(axis.text.x=element_text(size=8))+labs(y= "Inflation_Rates", x = "Years")

uk<-subset(df,df$series=='United.Kingdom')
ggplot(uk, aes(Years, value)) +
  geom_line(aes(colour = series))+geom_text(
    label=uk$series,
    nudge_x=0.45, nudge_y=0.1,
    check_overlap=T,
    size=3)+ theme(legend.position="none")+geom_point()+theme(axis.text.x=element_text(size=8))+labs(y= "Inflation_Rates", x = "Years")

germ<-subset(df,df$series=='Germany')
ggplot(germ, aes(Years, value)) +
  geom_line(aes(colour = series))+geom_text(
    label=germ$series,
    nudge_x=0.45, nudge_y=0.1,
    check_overlap=T,
    size=3)+ theme(legend.position="none")+geom_point()+theme(axis.text.x=element_text(size=8))+labs(y= "Inflation_Rates", x = "Years")

jap<-subset(df,df$series=='Japan')
ggplot(jap, aes(Years, value)) +
  geom_line(aes(colour = series))+geom_text(
    label=jap$series,
    nudge_x=0.45, nudge_y=0.1,
    check_overlap=T,
    size=3)+ theme(legend.position="none")+geom_point()+theme(axis.text.x=element_text(size=8))+labs(y= "Inflation_Rates", x = "Years")

thai<-subset(df,df$series=='Thailand')
ggplot(thai, aes(Years, value)) +
  geom_line(aes(colour = series))+geom_text(
    label=thai$series,
    nudge_x=0.45, nudge_y=0.1,
    check_overlap=T,
    size=3)+ theme(legend.position="none")+geom_point()+theme(axis.text.x=element_text(size=8))+labs(y= "Inflation_Rates", x = "Years")

china<-subset(df,df$series=='China')
ggplot(china, aes(Years, value)) +
  geom_line(aes(colour = series))+geom_text(
    label=china$series,
    nudge_x=0.45, nudge_y=0.1,
    check_overlap=T,
    size=3)+ theme(legend.position="none")+geom_point()+theme(axis.text.x=element_text(size=8))+labs(y= "Inflation_Rates", x = "Years")

russ<-subset(df,df$series=='Russian.Federation')
ggplot(russ, aes(Years, value)) +
  geom_line(aes(colour = series))+geom_text(
    label=russ$series,
    nudge_x=0.45, nudge_y=0.1,
    check_overlap=T,
    size=3)+ theme(legend.position="none")+geom_point()+theme(axis.text.x=element_text(size=8))+labs(y= "Inflation_Rates", x = "Years")

arg<-subset(df,df$series=='Argentina')
ggplot(arg, aes(Years, value)) +
  geom_line(aes(colour = series))+geom_text(
    label=arg$series,
    nudge_x=0.45, nudge_y=0.1,
    check_overlap=T,
    size=3)+ theme(legend.position="none")+geom_point()+theme(axis.text.x=element_text(size=8))+labs(y= "Inflation_Rates", x = "Years")

vene<-subset(df,df$series=='Venezuela..RB')
ggplot(vene, aes(Years, value)) +
  geom_line(aes(colour = series))+geom_text(
    label=vene$series,
    nudge_x=0.45, nudge_y=0.1,
    check_overlap=T,
    size=3)+ theme(legend.position="none")+geom_point()+theme(axis.text.x=element_text(size=8))+labs(y= "Inflation_Rates", x = "Years")

#------------
# create data sets
df<-read.csv('/home/deviancedev01/Desktop/interesting_exps/ggplot_worldecon/Global Dataset of Inflation.csv',
             check.names = FALSE)
#column.names<-df$Country
df<-na.omit(df)
df<-df %>% 
  group_by(Country) %>% 
  reframe(across(where(is.character)),
          across(where(is.numeric), sum, na.rm = T))
df<-data.frame(unique(df),check.names = FALSE)
row.names(df)<-df$Country
df<-df[,-c(1:2)]
df<-data.frame(t(df))
set.seed(4)
corr.matrix <- cor(df)
corrplot(corr.matrix,
         method = 'circle',
         tl.cex = 0.35,
         bg='white')
