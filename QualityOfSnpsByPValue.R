library(ggplot2)
library(cowplot)
library(reshape2)

NamePop<-read.table('~/genomes/genomes/PopNames/Name.BigPop.Pop.txt',col.names=c('Name', 'BigPop', 'Pop'))
NamePop<-NamePop[which(NamePop$Pop != 'NAG'),]

colors<-read.table("~/Dropbox/LukeTemp/NamePopBigPop_color.html.txt", comment.char='>',col.names=c('Pop','hex','BigPop'))
MyColour <- as.character(colors $hex)
names(MyColour) <- colors $Pop

Meta<-read.table("~/Dropbox/LukeTemp/1000GenomesMetaData.txt")
Meta$SUBMISSION.DATE<-as.Date(Meta$SUBMISSION.DATE, "%Y-%m-%d")  
SubMeta<-unique(Meta[c('Name','Pop', 'average_quality_of_mapped_bases')])

qualPop<-aggregate(.~Pop,SubMeta[c('Pop', 'average_quality_of_mapped_bases')],mean)
meanQual<-as.character(qualPop$average_quality_of_mapped_bases)
names(meanQual)<-qualPop$Pop

allPops<-data.frame()
for(pop in unique(NamePop$Pop)){
print(pop)	
fname=paste('~/Documents/GWAS_Qual/Final_GWAS/1kGP_GenomeWide_',pop,'.frq', sep='')
frq <-read.table(fname)
frq <-frq[which(frq$MAF > 0),]
frq$Pop<-pop

GWAS<-read.table(paste('~/Documents/GWAS_Qual/Final_GWAS/1kGP_GenomeWide_',pop,'_INT_ci_sig6.assoc.linear',sep=''),col.names=c('CHR','SNP','BP','A1','TEST','NMISS','BETA','SE','L95','U95','STAT','P'))
GWAS $Plog10<--log10(GWAS $P)
GWAS=merge(GWAS,frq, by='SNP')
#GWAS=merge(GWAS, melt.GT, by.x=c('CHR', 'BP'), by.y=c('Chr','Pos'))
GWAS=merge(GWAS,SubMeta, by.x='Pop', by.y='Name')
GWAS=GWAS[,c('CHR','BP','SNP','BETA','Plog10','Freq','Pop','average_quality_of_mapped_bases')]

G<-unique(GWAS[,c('CHR','BP','average_quality_of_mapped_bases','Freq','Plog10','Pop')])
Mean<-aggregate(.~CHR+BP+Freq+Plog10+Pop, G,mean)
allPops<-rbind(allPops,Mean)
ggplot(Mean,aes(x= average_quality_of_mapped_bases, y=Plog10, color=Freq))+geom_point(shape=1)+ggtitle(pop)
ggsave(paste('/Users/luke/Documents/GWAS_Qual/PvalueByQuality_',pop,'_Jun28.tiff',sep=''),height=5,width=7)
}

pl=list()
for(pop in unique(NamePop$Pop)){
if(pop != 'BEB'){
	suba<-allPops[which(allPops$Pop==pop),]
	m=qualPop[which(qualPop$Pop==pop),]$average_quality_of_mapped_bases
	pl[[pop]]<-ggplot(suba,aes(x= average_quality_of_mapped_bases, y=Plog10, color=Freq))+geom_vline(xintercept=m, color='blue',linetype=2)+geom_hline(yintercept=8, color='blue',linetype=2)+guides(color=F)+scale_color_gradient(low='green',high='red')+geom_point(shape=1)+xlab(pop)
	}}
pll<-plot_grid(pl[['JPT']],pl[['CHB']],pl[['PUR']],pl[['YRI']],pl[['ASW']],pl[['LWK']],pl[['MXL']],pl[['GWD']],pl[['CHS']],pl[['CEU']],pl[['GBR']],pl[['FIN']],pl[['MSL']],pl[['ESN']],pl[['STU']],pl[['CLM']],pl[['TSI']],pl[['ITU']],pl[['PJL']],pl[['ACB']],pl[['KHV']],pl[['PEL']],pl[['CDX']],pl[['IBS']])	
	
ggsave('/Users/luke/Documents/GWAS_Qual/PvalueByQualityAllPops_mean_1.tiff',pll, height=8,width=12)

aggAll<-aggregate(.~CHR+BP+Freq+Plog10+Pop, allPops,mean)
aggAll<-merge(aggAll,qualPop, by='Pop')
ggplot(aggAll,aes(x= average_quality_of_mapped_bases, y=Plog10, color=Freq))+geom_point(shape=1)