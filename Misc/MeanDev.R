##Manhattan
library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(multtest)
library(cowplot)

dir='~/Documents/Regression/'
out='~/Documents/Regression/'

for(chrom in seq(1,5)){
	print(chrom)
	#read in the regressions for each pop
	fileNames = list.files(path=dir,pattern=paste("CHR",chrom,"_NOT.Regression",sep=''), full.names = T)
	if(length(fileNames) != 0){
	Reg = do.call(rbind, lapply(fileNames, function(x) read.table(x, header=T, sep=' ')))
	
	print(head(Reg))
	#had duplicate headers, so had to remove them and reset columns to numeric
	Reg<-Reg[which(Reg$Pos != 'Pos'),]
	Reg$Chr<-as.numeric(as.character(Reg$Chr))
	Reg$Pos<-as.numeric(as.character(Reg$Pos))
	Reg$dev<-as.numeric(as.character(Reg$dev))
	#chisquared random variable
	#taking the sum of the deviances (mean is to take into account the DF)
	Reg <- Reg %>% group_by(Chr,Pos) %>% summarize(SumDev = sum(dev), MeanDev = mean(dev), Count = n())
	
	write.table(Reg,paste(out,chrom,'NOT_MeanDev.csv',sep=''))
	}
	}
