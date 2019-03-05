## Deal with Regression output

library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(multtest)
library(cowplot)
library(reshape2)

Regression <-fread('~/Documents/QualityPaper/Misc/GenomeWide.Regression.txt', colClasses='numeric',col.names=c('CHROM','POS','ID','dev.p'))

Regression $ID<-as.numeric(as.character(Regression $ID))
Regression $CHROM<-as.numeric(as.character(Regression $CHROM))
Regression $POS<-as.numeric(as.character(Regression $POS))
Regression $dev.p<-as.numeric(as.character(Regression $dev.p))
Regression <- Regression[complete.cases(Regression),]


#read in list of global singletons
singles<-fread(
	'/Users/luke/genomes/genomes/1kGP_NotFiltered/GenomeWide_doubles.frq',
	col.names=c('CHROM','POS'))
singles$drop<-'Yes' #dummy variable for merging

#join the singles to the Reg to get info of what sites to remove
Regression <-left_join(Regression, singles, by=c('CHROM','POS'))
rm(singles)

#this removes the singletons
Regression <-Regression[which(is.na(Regression$drop)),]
Regression$drop<-NULL

##read Repeats
Repeats <- fread('~/Documents/QualityPaper/Misc/Reg_strictMask_nestedRepeats.POS',fill=T,col.names=c('CHROM','POS'))
Repeats$Repeat='yes'

##read Indels
Indels <- fread('/Users/luke/Documents/QualityPaper/Misc/Reg_indels.POS', skip=1 ,col.names=c('CHROM','POS'))
Indels$Indel = 'yes'

#read in rsIDs (used for GWAS catalogue search)
rsID<-fread('/Users/luke/genomes/genomes/hg19/phase3/rsID.txt')
colnames(rsID)<-c('CHROM','POS','rsID')

Regression <- left_join(Regression, Repeats, by=c('CHROM','POS'))

Regression <- left_join(Regression, Indels, by=c('CHROM','POS'))

Regression <- left_join(Regression, rsID, by=c('CHROM','POS'))

#free up memory
rm(Repeats)
rm(Indels)
rm(rsID)

write.table(Regression, file='/Users/luke/Documents/QualityPaper/Misc/TotalRegressionJustPop.csv',row.names=F)

Regression<-fread('/Users/luke/Documents/QualityPaper/Misc/TotalRegressionJustPop.csv')

Repeat_Indels <- unique(Regression[which(Regression$Indel=='yes' 
								& Regression$Repeat=='yes'),])

NORepeat_Indels <- unique(Regression[which(Regression$Indel=='yes' 
								& is.na(Regression$Repeat)),])
								
Repeat_SNPs <- unique(Regression[which(is.na(Regression$Indel) 
								& Regression$Repeat=='yes'),])
								
NORepeat_SNPs <- unique(Regression[which(is.na(Regression$Indel) 
								& is.na(Regression$Repeat)),])			

##RICK Adjustment
# two-stage Benjamini & Hochberg (2006) step-up FDR-controlling
TSBH <-function(df){	
procedures <- c( "TSBH") 
#0.01
df$p <- 10^-(df$dev.p)
adjusted <- mt.rawp2adjp(df$p, procedures, alpha = 0.01)
adj <- as.data.frame(adjusted$adj[order(adjusted$index), ])
df$log10P_0.01 <- -log10(adj$TSBH_0.01)
df$P_0.01 <- adj$TSBH_0.01

return(df)
}
#perform FDR separately for each category
Repeat_Indels <-TSBH(Repeat_Indels)
NORepeat_Indels <-TSBH(NORepeat_Indels)
Repeat_SNPs <-TSBH(Repeat_SNPs)
NORepeat_SNPs <-TSBH(NORepeat_SNPs)

#write significant positions
writeSig <- function(df, Name){
write.table(df[which(df$log10P_0.01>-log10(0.01)),]$rsID, file=paste('~/Documents/QualityPaper/sig/',Name,'_rsID.txt',sep=''), quote=F, col.names=F, row.names=F)

write.table(df[which(df$log10P_0.01>-log10(0.01)),c('CHROM','POS')], file=paste('~/Documents/QualityPaper/sig/',Name,'_POS.txt',sep=''), quote=F, col.names=F, row.names=F)
}

writeSig(Repeat_Indels, 'Repeat_Indels')
writeSig(NORepeat_Indels, 'NORepeat_Indels')
writeSig(Repeat_SNPs, 'Repeat_SNPs')
writeSig(NORepeat_SNPs, 'NORepeat_SNPs')

#combine them back together and write output file
Regression <- bind_rows(list(Repeat_Indels, NORepeat_Indels, Repeat_SNPs, NORepeat_SNPs))
write.table(Regression[which(Regression$log10P_0.01>-log10(0.01)),c('CHROM','POS','rsID','ID','dev.p','log10P_0.01')], file = '~/Documents/QualityPaper/sig/SignificantVariants.txt',row.names=F, quote=F, sep=',')
write.table(Regression[which(Regression$log10P_0.01>-log10(0.01)),c('CHROM','POS')], file = '~/Documents/QualityPaper/sig/SignificantVariantsCHROM_POS.txt',row.names=F, quote=F,col.names=F,sep='\t')

#Make the manhattan plots
makePlot <- function(df,Name, fileName){
#subset the data to deal with sites -log10(p) > 20
sig.20<-df[which((df$log10P_0.01 >= 20)),]
sig.20$log10P_0.01=20
sig.19<-df[which((df$log10P_0.01 < 20)),]
nSig <- nrow(df[which(df$log10P_0.01 > -log10(0.01)),])

#make the plot
ggplot(sig.19, aes(x=POS, y = log10P_0.01, color=as.factor(CHROM)))+
facet_grid(~CHROM, scales='free_x', space='free_x', switch='x')+
scale_color_manual(values = rep(c("grey", "black"), 22 )) +
scale_y_continuous(expand=c(0,0))+
geom_point(alpha=0.3, size=1)+theme_classic()+
geom_point(data=sig.20, aes(x=POS, y = log10P_0.01),shape=1)+
geom_hline(yintercept=-log10(0.01), color='blue')+
labs(y='-log10(p)', x='Chromosome', title=Name, subtitle=paste('Number of loci tested =',nrow(df),'\n','Number of significant loci =',nSig))+
theme(plot.title = element_text(hjust = 0.5), 
plot.subtitle = element_text(hjust = 0.5), 
axis.text.x=element_blank(), 
axis.ticks.x=element_blank(),
axis.line.x=element_blank(), 
strip.background = element_blank(),
strip.text.x = element_text(size = 6))+
guides(color=F)+expand_limits(y=c(0,21))
#save it
dir='~/Documents/QualityPaper/Figures/ManhattanPlot_'
ggsave(paste(dir, fileName,'.jpg',sep=''), height=5, width=10)
ggsave(paste(dir, fileName,'.tiff',sep=''), height=5, width=10)

}

makePlot(Repeat_Indels, 'Repeat Indels', 'Repeat_Indels')
makePlot(NORepeat_Indels, 'Non Repeat Indels', 'NORepeat_Indels')
makePlot(Repeat_SNPs, 'Repeat SNPs', 'Repeat_SNPs')
makePlot(NORepeat_SNPs, 'Non Repeat SNPs', 'NORepeat_SNPs')


#### GWAS publications,
catalog <- fread('~/Documents/GWAS_Qual/gwas_catalog_cut.txt',fill=TRUE, sep='\t')
catalog$rsID <- catalog$SNPS



catalog <- merge(catalog, Regression[,c('CHROM','POS', 'rsID','log10P_0.01')], by='rsID')

write.table(catalog[which(catalog$log10P_0.01 > -log10(0.01)),], file = '~/Documents/QualityPaper/sig/SignificantPublications.txt',row.names=F, quote=F, sep=',')

table <- catalog[which(catalog$log10P_0.01 > -log10(0.01)),c('PUBMEDID','JOURNAL','SNPS','PVALUE_MLOG','log10P_0.01')]
print(table[order(-log10P_0.01),])

###OMNI CHIP
omni <- fread('~/Downloads/HumanOmni2-5-8-v1-2-A-b138-rsIDs.txt',header=T)
omni <- merge(omni, Regression, by.x='RsID', by.y='rsID')
omni <- omni[which(omni$RsID != '.'),]

ggplot(omni[which(omni$log10P_0.01 > 2),], aes(x=ID, y= log10P_0.01))+geom_point(shape=1)+labs(x='Allele Frequency',y='Adjusted -log10(p)')+ggtitle('219 SNPs present in the Omni 2.5 Microarray')
ggsave('~/Documents/QualityPaper/Figures/Omni_AF.jpg', height=5, width=10)
####

