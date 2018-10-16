##Total Regression
library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(multtest)
library(cowplot)

##read Regression
Regression<-fread('~/Documents/QualityPaper/Misc/Reg_ALL_SITES_5col.txt', fill=T)
Regression <-Regression[which(Regression $Pos != 'Pos'),]
Regression $Chr<-as.numeric(as.character(Regression $Chr))
Regression $Pos<-as.numeric(as.character(Regression $Pos))
Regression $SumDev<-as.numeric(as.character(Regression $SumDev))
Regression $MeanDev<-as.numeric(as.character(Regression $MeanDev))
Regression $Count<-as.numeric(as.character(Regression $Count))

#read in list of global singletons
singles<-fread(
	'/Users/luke/genomes/genomes/1kGP_NotFiltered/GenomeWide_doubles.frq',
	col.names=c('Chr','Pos'))
singles$drop<-'Yes' #dummy variable for merging

#join the singles to the Reg to get info of what sites to remove
Regression <-left_join(Regression, singles, by=c('Chr','Pos'))
rm(singles)

#this removes the singletons
Regression <-Regression[which(is.na(Regression$drop)),]
Regression$drop<-NULL

##read Repeats
Repeats <- fread('~/Documents/QualityPaper/Misc/Reg_strictMask_nestedRepeats.POS',fill=T,col.names=c('Chr','Pos'))
Repeats$Type='Repeat'

##read Non-Repeats
NoRepeats <- fread('~/Documents/QualityPaper/Misc/Reg_strictMask_NOnestedRepeats.POS',fill=T,col.names=c('Chr','Pos'))
NoRepeats$Type='NOT_Repeat'

##read Indels
Indels <- fread('/Users/luke/Documents/QualityPaper/Misc/Reg_indels.POS',header=T)
Indels$Indel = 'yes'

#read in rsIDs (used for GWAS catalogue search)
rsID<-fread('/Users/luke/genomes/genomes/hg19/phase3/rsID.txt')
colnames(rsID)<-c('Chr','Pos','rsID')

Regression <- left_join(Regression, Repeats, by=c('Chr','Pos'))

Regression <- left_join(Regression, NoRepeats, by=c('Chr','Pos'))

Regression <- left_join(Regression, Indels, by=c('Chr','Pos'))

Repeat_Indels <- Regression[which(Regression$Indel=='yes' 
								& Regression$Type.x=='Repeat'),]

NORepeat_Indels <- Regression[which(Regression$Indel=='yes' 
								& Regression$Type.y=='NOT_Repeat'),]
								
Repeat_SNPs <- Regression[which(is.na(Regression$Indel) 
								& Regression$Type.x=='Repeat'),]
								
NORepeat_SNPs <- Regression[which(is.na(Regression$Indel) 
								& Regression$Type.y=='NOT_Repeat'),]			

##RICK Adjustment
# two-stage Benjamini & Hochberg (2006) step-up FDR-controlling
TSBH <-function(df){	
procedures <- c( "TSBH") 
#0.01
df$p <- pchisq(df$SumDev, df$Count, lower.tail=F)
adjusted <- mt.rawp2adjp(df$p, procedures, alpha = 0.01)
adj <- as.data.frame(adjusted$adj[order(adjusted$index), ])
df$log10P_0.01 <- -log10(adj$TSBH_0.01)
df$P_0.01 <- adj$TSBH_0.01

return(df)
}

Repeat_Indels <-TSBH(Repeat_Indels)
NORepeat_Indels <-TSBH(NORepeat_Indels)
Repeat_SNPs <-TSBH(Repeat_SNPs)
NORepeat_SNPs <-TSBH(NORepeat_SNPs)

write.table(Repeat_Indels, file='/Users/luke/Documents/QualityPaper/Misc/Repeat_Indels.csv',row.names=F)
write.table(NORepeat_Indels, file='/Users/luke/Documents/QualityPaper/Misc/NORepeat_Indels.csv',row.names=F)	
write.table(Repeat_SNPs, file='/Users/luke/Documents/QualityPaper/Misc/Repeat_SNPs.csv',row.names=F)	
write.table(NORepeat_SNPs, file='/Users/luke/Documents/QualityPaper/Misc/NORepeat_SNPs.csv',row.names=F)

makePlot <- function(df,Name){
sig.20<-df[which((df $log10P_0.01 >= 20)&(df $Count > 1)),]
sig.20$log10P_0.01=20
sig.19<-df[which((df $log10P_0.01 < 20)&(df $Count > 1)),]

ggplot(sig.19, aes(x=Pos, y = log10P_0.01, color=as.factor(Chr)))+
facet_grid(~Chr, scales='free_x', space='free_x', switch='x')+
scale_fill_manual (values=getPalette(colourCount))+
scale_y_continuous(expand=c(0,0))+
geom_point(alpha=0.3, size=1)+theme_classic()+
geom_point(data=sig.20, aes(x=Pos, y = log10P_0.01),shape=1)+
geom_hline(yintercept=-log10(0.01), color='blue')+
labs(y='-log10(p)', x='Chromosome', title=Name, subtitle=paste('Number of positions = ',nrow(df)))+
theme(plot.title = element_text(hjust = 0.5), 
plot.subtitle = element_text(hjust = 0.5), 
axis.text.x=element_blank(), 
axis.ticks.x=element_blank(),
axis.line.x=element_blank(), 
strip.background = element_blank(),
strip.text.x = element_text(size = 6))+
guides(color=F)+expand_limits(y=c(0,21))

dir='~/Documents/QualityPaper/Figures/ManhattanPlot_'
ggsave(paste(dir,Name,'.jpg',sep=''), height=5, width=10)
ggsave(paste(dir,Name,'.tiff',sep=''), height=5, width=10)
}

makePlot(Repeat_Indels, 'Repeat_Indels')

makePlot(NORepeat_Indels, 'NORepeat_Indels')

makePlot(Repeat_SNPs, 'Repeat_SNPs')

makePlot(NORepeat_SNPs, 'NORepeat_SNPs')
