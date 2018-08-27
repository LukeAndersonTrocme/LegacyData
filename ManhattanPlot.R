##Manhattan
library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(multtest)

dir='/Volumes/gravel/luke_projects/1000Genomes/Regression/'
out='/Volumes/gravel/luke_projects/1000Genomes/MeanDev/'

if(F){ #this is the loop that took the mean of the deviances
for(chrom in seq(1,22)){
	print(chrom)
	#read in the regressions for each pop
	fileNames = list.files(path=dir,pattern=paste("CHR",chrom,".Regression",sep=''), full.names = T)
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
	
	write.table(Reg,paste(out,chrom,'_MeanDev.csv',sep=''))
	}
}
}

#read in rsIDs (used for GWAS catalogue search)
rsID<-fread('/Users/luke/genomes/genomes/1kGP_4bedFiltered/1kGP_GenomeWide_Chr.Pos.rsID.txt')
colnames(rsID)<-c('Chr','Pos','rsID')

#read in list of global singletons
singles<-fread(
	'~/Documents/QualityPaper/Misc/Singletons_Chr.Pos.rsID.txt',
	col.names=c('Chr','Pos','rsID'))
singles$drop<-'Yes' #dummy variable for merging
	
doubles<-fread(
	'/Users/luke/Documents/QualityPaper/Misc/Doubletons.txt',
	col.names=c('Chr','Pos','rsID'))
doubles$drop<-'Yes' #dummy variable for merging

OnesTwos<-rbind(singles, doubles)
rm(singles)
rm(doubles)

#read in mean deviances for all pops
Reg <- fread('/Volumes/gravel/luke_projects/1000Genomes/MeanDev/AllPops_Dev.csv', fill=T)
colnames(Reg)<-c('V1','Chr','Pos','SumDev','MeanDev','Count')

#had duplicate headers, so had to remove them and reset columns to numeric
Reg$V1<-NULL
Reg<-Reg[which(Reg$Chr != 'Pos'),]
Reg$Chr<-as.numeric(as.character(Reg$Chr))
Reg$Pos <-as.numeric(as.character(Reg$Pos))
Reg$SumDev <-as.numeric(as.character(Reg$SumDev))
Reg$MeanDev <-as.numeric(as.character(Reg$MeanDev))

#join the singles to the Reg to get info of what sites to remove
Reg<-left_join(Reg, OnesTwos, by=c('Chr','Pos'))

rm(singles)
#get p values from sum of deviances
Reg$plog10 <- - pchisq(Reg$SumDev, Reg$Count, lower.tail=F, log.p=T)/log(10)
Reg$p <- pchisq(Reg$SumDev, Reg$Count, lower.tail=F)

#write the table to make it easier to load later on
write.table(Reg, file='~/Documents/QualityPaper/Misc/Reg.txt', quote=F, row.names=F)
#Reg <- fread('~/Documents/QualityPaper/Misc/Reg.txt')

#this removes the singletons
NoSingles<-Reg[which(is.na(Reg$drop)),]
NoSingles<-left_join(NoSingles, rsID, by=c('Chr','Pos')) #get rsID for these sites

##RICK Adjustment
procedures <- c( "TSBH") # two-stage Benjamini & Hochberg (2006) step-up FDR-controlling
#0.05
adjusted <- mt.rawp2adjp(NoSingles$p, procedures, alpha = 0.05)
adj <- as.data.frame(adjusted$adj[order(adjusted$index), ])
NoSingles$adjusted.05 <- -log10(adj$TSBH_0.05)
NoSingles$adjustedP.05 <- adj$TSBH_0.05
#0.01
adjusted <- mt.rawp2adjp(NoSingles$p, procedures, alpha = 0.01)
adj <- as.data.frame(adjusted$adj[order(adjusted$index), ])
NoSingles$adjusted.01 <- -log10(adj$TSBH_0.01)
NoSingles$adjustedP.01 <- adj$TSBH_0.01

adj.sig.05<-NoSingles[which((NoSingles$adjusted.05 >= -log10(0.05))&(NoSingles$Count > 1)),]
adj.sig.01<-NoSingles[which((NoSingles$adjusted.01 >= -log10(0.01))&(NoSingles$Count > 1)),]

write.table(adj.sig.05, 
	'~/Documents/QualityPaper/Significant0.5SNPs_adjusted.txt', 
	quote=F, row.names=F)
write.table(adj.sig.01, 
	'~/Documents/QualityPaper/Significant0.1SNPs_adjusted.txt', 
	quote=F, row.names=F)
write.table(adj.sig.05$rsID.y, 
	'~/Documents/QualityPaper/Significant0.5SNPs_adjusted_rsID.txt', 
	quote=F, row.names=F, col.names=F)
write.table(adj.sig.01$rsID.y, 
	'~/Documents/QualityPaper/Significant0.1SNPs_adjusted_rsID.txt', 
	quote=F, row.names=F, col.names=F)
write.table(adj.sig.05[,c('Chr','Pos')],
	'~/Documents/QualityPaper/Significant0.5SNPs_adjusted_POS.txt', 
	quote=F, row.names=F, col.names=F,  sep='\t')
write.table(adj.sig.01[,c('Chr','Pos')], 
	'~/Documents/QualityPaper/Significant0.1SNPs_adjusted_POS.txt', 
	quote=F, row.names=F, col.names=F, sep='\t')
