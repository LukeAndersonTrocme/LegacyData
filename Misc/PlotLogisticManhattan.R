##Plot Logistic P values
library(ggplot2)
library(data.table)
#/Users/luke/Documents/Regression/CHR22.Logistic.txt
Logistic<-fread('~/Documents/Regression/CHR22.Regression_JPT.csv', header=T)

ggplot(Logistic,aes(x=p), color='grey70')+geom_density()+theme_classic()

ggplot(Logistic,aes(x=Pos,y=Plog10, color=Pop))+geom_density(stat='density')+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())