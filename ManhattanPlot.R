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
singles<-fread('~/Documents/QualityPaper/Misc/Singletons_Chr.Pos.rsID.txt', col.names=c('Chr','Pos','rsID'))
singles$single<-'Yes' #dummy variable for merging

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
Reg<-left_join(Reg, singles, by=c('Chr','Pos'))

rm(singles)
#get p values from sum of deviances
Reg$plog10 <- - pchisq(Reg$SumDev, Reg$Count, lower.tail=F, log.p=T)/log(10)
Reg$p <- pchisq(Reg$SumDev, Reg$Count, lower.tail=F)

#write the table to make it easier to load later on
write.table(Reg, file='~/Documents/QualityPaper/Misc/Reg.txt', quote=F, row.names=F)
#Reg <- fread('~/Documents/QualityPaper/Misc/Reg.txt')
#this removes the singletons
NoSingles<-Reg[which(is.na(Reg$single)),]
NoSingles<-left_join(NoSingles, rsID, by=c('Chr','Pos')) #get rsID for these sites
##RICK Adjustment
procedures <- c( "TSBH") # two-stage Benjamini & Hochberg (2006) step-up FDR-controlling
adjusted <- mt.rawp2adjp(NoSingles$p, procedures, alpha = 0.05)
adj<-as.data.frame(adjusted$adj[order(adjusted$index), ])
NoSingles$adjusted<--log10(adj$TSBH_0.05)
NoSingles$adjustedP<-adj$TSBH_0.05

adj.less<-NoSingles[which((NoSingles$adjusted < 20)&(NoSingles$Count > 1)),]

adj.20<-NoSingles[which((NoSingles$adjusted >= 20)&(NoSingles$Count > 1)),]
adj.20$adjusted=20
adj.sig<-NoSingles[which((NoSingles$adjusted >= -log10(0.05))&(NoSingles$Count > 1)),]

write.table(adj.sig, '~/Documents/QualityPaper/Significant6SNPs_adjusted.txt', quote=F, row.names=F)

ggplot(adj.less, aes(x=Pos, y = adjusted, color=as.factor(Chr)))+
geom_point(data=adj.20, aes(x=Pos, y = plog10,),shape=1)+
facet_grid(~Chr, scales='free_x', space='free_x', switch='x')+
scale_fill_manual (values=getPalette(colourCount))+
scale_y_continuous(expand=c(0,0))+
geom_point(alpha=0.3, size=1)+theme_classic()+
geom_hline(yintercept=-log10(0.05), color='blue')+
labs(y='adjusted -log10(p)', x='Chromosome')+
theme(plot.title = element_text(hjust = 0.5), 
plot.subtitle = element_text(hjust = 0.5), 
axis.text.x=element_blank(), 
axis.ticks.x=element_blank(),
axis.line.x=element_blank(), 
strip.background = element_blank(),
strip.text.x = element_text(size = 6))+
guides(color=F)+expand_limits(y=c(0,21))
ggsave('~/Documents/QualityPaper/ManhattanPlot_adjusted.jpg', height=5, width=10)
