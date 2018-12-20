##Manhattan
library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(multtest)
library(cowplot)

dir='~/Documents/GenomeWideRegression/'
out='~/Documents/GenomeWideRegression/'

Reg = fread('~/Documents/GenomeWideRegression/ALLPOP_GenomeWide.Regression.csv', header=F, sep=' ', col.names=c('Chr','Pos','Pop','dev'))
Reg$Chr<-as.numeric(as.character(Reg$Chr))
Reg$Pos<-as.numeric(as.character(Reg$Pos))
Reg$dev<-as.numeric(as.character(Reg$dev))


for(i in 1:22){
print(i)	
Sub<- Reg[which(Reg$Chr == i),]
#chisquared random variable
#taking the sum of the deviances (mean is to take into account the DF)
Sub <- Sub %>% group_by(Chr,Pos) %>% summarize(SumDev = sum(dev), MeanDev = mean(dev), Count = n())
#write results to file
if(i == 1){
	write.table(Sub,file='~/Documents/GenomeWideRegression/GenomeWide.MeanDev_ALLPOPS.csv', quote=F, row.names=F, append=FALSE)}
if(i > 1){
    write.table(Sub,file='~/Documents/GenomeWideRegression/GenomeWide.MeanDev_ALLPOPS.csv', quote=F, row.names=F, append=TRUE)}
}


	

