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
