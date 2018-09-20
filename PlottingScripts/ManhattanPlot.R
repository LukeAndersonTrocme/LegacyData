##Manhattan
library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(multtest)

#dir='/Volumes/gravel/luke_projects/1000Genomes/Regression/'
#out='/Volumes/gravel/luke_projects/1000Genomes/MeanDev/'
dir='~/Documents/Regression/'
out='~/Documents/Regression/'

if(F){ #this is the loop that took the mean of the deviances
for(chrom in seq(4,22)){
	print(chrom)
	#read in the regressions for each pop
	fileNames = list.files(path=dir,pattern=paste("CHR",chrom,"_indels.Regression",sep=''), full.names = T)
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
	
	write.table(Reg,paste(out,chrom,'indels_MeanDev.csv',sep=''))
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
Reg <- fread('/Users/luke/Documents/Regression/ALL_indels_MeanDev.csv', fill=T)
#Reg <- fread('/Volumes/gravel/luke_projects/1000Genomes/MeanDev/AllPops_Dev.csv', fill=T)
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
write.table(Reg, file='~/Documents/QualityPaper/Misc/Reg_indels.txt', quote=F, row.names=F)
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

sig.20<-NoSingles[which((NoSingles $adjusted.01 >= 20)&(NoSingles $Count > 1)),]
sig.20$adjusted.01=20
sig.19<-NoSingles[which((NoSingles $adjusted.01 < 20)&(NoSingles $Count > 1)),]

ggplot(sig.19, aes(x=Pos, y = adjusted.01, color=as.factor(Chr)))+
facet_grid(~Chr, scales='free_x', space='free_x', switch='x')+
scale_fill_manual (values=getPalette(colourCount))+
scale_y_continuous(expand=c(0,0))+
geom_point(alpha=0.3, size=1)+theme_classic()+
geom_point(data=sig.20, aes(x=Pos, y = adjusted.01),shape=1)+
geom_hline(yintercept=-log10(0.01), color='blue')+
labs(y='-log10(p)', x='Chromosome')+
theme(plot.title = element_text(hjust = 0.5), 
plot.subtitle = element_text(hjust = 0.5), 
axis.text.x=element_blank(), 
axis.ticks.x=element_blank(),
axis.line.x=element_blank(), 
strip.background = element_blank(),
strip.text.x = element_text(size = 6))+
guides(color=F)+expand_limits(y=c(0,21))

ggsave('~/Documents/QualityPaper/Figures/ManhattanPlot_indels.jpg', height=5, width=10)
ggsave('~/Documents/QualityPaper/Figures/ManhattanPlot_indels.tiff', height=5, width=10)


#make min for dev
#geom_vline(xintercept=
#min(adj.sig.01[which(adj.sig.01$Count==df),]$SumDev),
#color='grey60',linetype=2)+
#labs(x=paste('df =',df), y='density')+

MakePlot<-function(file){
	ggplot(file, aes(x=SumDev))+
				stat_function(aes(x=seq(0,10,
				length=nrow(file))),fun=dchisq, 
				args=list(df=df), color='blue')+
				geom_density()+xlim(c(0,80))+
				geom_vline(xintercept=
				qchisq(0.000001, df=df, lower.tail=F),
				color='grey60',linetype=2)+
				labs(x='deviance', y='density',subtitle=paste('df =',df))
}

LP<-list()
for(df in seq(1,26)){
	file <- NoSingles[which(NoSingles$Count==df),]
	LP[[df]] <- MakePlot(file)	
	}

pf <- grid.arrange(LP[[1]],LP[[2]], LP[[3]],LP[[4]],LP[[5]],LP[[6]],LP[[7]],LP[[8]],LP[[9]],LP[[10]], LP[[11]],LP[[12]],LP[[13]],LP[[14]],LP[[15]],LP[[16]],LP[[17]],LP[[18]],LP[[19]],LP[[20]], LP[[21]],LP[[22]],LP[[23]], LP[[24]],LP[[25]],LP[[26]])

ggsave('~/Documents/QualityPaper/Figures/AllDeviances.jpg',pf, height=16,width=22)
ggsave('~/Documents/QualityPaper/Figures/AllDeviances.tiff',pf, height=16,width=22)
	
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

##indels
write.table(adj.sig.01[which(is.na(adj.sig.01$rsID.y)==F),]$rsID.y, 
	'~/Documents/QualityPaper/Significant0.1INDELs_adjusted_rsID.txt', 
	quote=F, row.names=F, col.names=F)