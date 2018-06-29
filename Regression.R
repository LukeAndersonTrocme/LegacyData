#Quality Genotype Regression
#java -jar ~/bin/snpEff_latest_core/snpEff/SnpSift.jar extractFields /Users/luke/Documents/QualityPaper/PublishedSNPs_header.vcf CHROM POS "GEN[*].GT"   | sed 's/\//|/g ; s/0|0/0/g ; s/0|1/1/g ; s/1|0/1/g ; s/1|1/2/g ; /|/d'| tail -n +2 > /Users/luke/Documents/QualityPaper/PublishedSNPs.Genotypes.txt

library(reshape2)
library(ggplot2)

NamePop<-read.table('~/genomes/genomes/PopNames/Name.BigPop.Pop.txt',col.names=c('Name', 'BigPop', 'Pop'))
NamePop<-NamePop[which(NamePop$Pop != 'NAG'),]

colors<-read.table("~/Dropbox/LukeTemp/NamePopBigPop_color.html.txt", comment.char='>',col.names=c('Pop','hex','BigPop'))
MyColour <- as.character(colors $hex)
names(MyColour) <- colors $Pop

Meta<-read.table("~/Dropbox/LukeTemp/1000GenomesMetaData.txt")
Meta$SUBMISSION.DATE<-as.Date(Meta$SUBMISSION.DATE, "%Y-%m-%d")  
SubMeta<-unique(Meta[c('Name','Pop', 'average_quality_of_mapped_bases')])

GT<-read.table('~/Documents/QualityPaper/PublishedSNPs.Genotypes.txt',header=F)
names(GT)=c('Chr','Pos', as.character(NamePop$Name))
#cut=c('Chr','Pos',as.character(NamePop[which(NamePop$Pop == pop),]$Name))
#GT<-GT[cut]

melt.GT<-melt(GT,id=c('Chr','Pos'))
melt.GT$value<-gsub('2','0',melt.GT$value)
GT.Qual<-merge(melt.GT, SubMeta, by.x='variable',by.y='Name')

#ggplot(GT.Qual, aes(x= average_quality_of_mapped_bases, y=value, color=Pop))+geom_point()+facet_grid(Chr~Pop)+scale_color_manual(values=MyColour)+guides(color=F)+theme_classic()

GT.Qual$p=''
for(pop in unique(GT.Qual$Pop)){
gtQual<-GT.Qual[which(GT.Qual$Pop==pop),]
gtQual $t=''
print(pop)
for(snp in unique(jpt$Pos)){
	sub= gtQual[which(gtQual $Pos==snp),]
	if(length(unique(sub$value))!=1){
	r1=glm(as.numeric(value) ~ average_quality_of_mapped_bases,family=binomial, data = sub)
	t1=paste("Adj R2 = ",signif(summary(r1)$adj.r.squared, 3)," P =",signif(summary(r1)$coef[1,4], 3))
	gtQual[which(gtQual $Pos==snp),]$t=t1
	GT.Qual[which((GT.Qual $Pos==snp)& GT.Qual $Pop==pop),]$p=-log10(summary(r1)$coef[1,4])
	
}
}
ggplot(gtQual, aes(x= average_quality_of_mapped_bases, y=value))+geom_point(shape=3)+facet_wrap(Pos~., ncol=3)+theme_classic()

#ggplot(gtQual, aes(x= average_quality_of_mapped_bases, y=value))+geom_point()+facet_wrap(t~., ncol=2)+theme_classic()+geom_smooth(method = "lm", se = FALSE,size=0.5)+labs(title=paste("Genotype vs Quality for",pop), y='Genotype',x='Average quality of mapped bases')+theme(plot.title = element_text(hjust= 0.5))+scale_y_continuous(breaks=c(0,1,2))

#ggsave(paste('~/Documents/Regression/PublishedSNPs_',pop,'.jpg', sep=''))
}

ye<-gtQual[which(gtQual$Pos == 55566507),]
r1=glm(as.numeric(value) ~ average_quality_of_mapped_bases,family=binomial, data = ye)
ggplot(ye, aes(x= average_quality_of_mapped_bases, y=value))+geom_point(shape=3)

pops<-list.files(path="/Users/luke/Documents/SigSnps",pattern=".frq$", full.names=T)

Freq<-data.frame()
for(f in seq(1,length(pops))){
pop<-unique(read.table(pops[f], row.names=NULL,header=T))
names(pop)=c('CHROM', 'POS', 'N_ALLELES','N_CHR', 'Y.ALLELE.FREQ','X.ALLELE.FREQ')
pop$Pop=substr(pops[f],31,33)
Freq<-rbind(Freq, pop)
}

Freq$X.ALLELE.FREQ<-gsub('T:','',Freq$X.ALLELE.FREQ)
Freq$X.ALLELE.FREQ<-gsub('C:','',Freq$X.ALLELE.FREQ)
Freq$X.ALLELE.FREQ<-gsub('G:','',Freq$X.ALLELE.FREQ)
Freq$X.ALLELE.FREQ<-gsub('A:','',Freq$X.ALLELE.FREQ)

perPop<-read.table('/Users/luke/Documents/SigSnps/PerPop.txt', col.names=c('CHROM','POS','Pop'))
perPop$Sig<-'Yes'

Freq<-merge(Freq,perPop, by=c('CHROM', 'POS','Pop'),all=T)

ordr<-c('FIN', 'GBR', 'CEU', 'IBS', 'TSI', 'CHS', 'CDX', 'CHB', 'JPT', 'KHV', 'GIH', 'STU', 'PJL', 'ITU', 'BEB', 'PEL', 'MXL', 'CLM', 'PUR', 'ASW', 'ACB', 'GWD', 'YRI', 'LWK', 'ESN', 'MSL')

Freq $ordered <- factor(Freq$Pop, levels =ordr)

plt<-merge(Freq, unique(GT.Qual[,c('Pos','Pop', 'p')]), by.x=c('POS','Pop'),by.y=c('Pos','Pop'))

ggplot(plt, aes(x=ordered, y= as.numeric(p), color=ordered, group=1))+facet_wrap(.~ POS,ncol=1)+geom_line()+geom_point(shape=15)+scale_color_manual(values=MyColour)+theme_classic()+theme(strip.background=element_blank(), strip.text=element_blank(), axis.title=element_blank(), axis.ticks.y=element_blank(),axis.line=element_blank(), axis.text.x=element_text(angle=90, size=5, vjust=-0.05))+guides(color=F)+geom_point(data=plt[complete.cases(plt),], shape=1, size=4, color='black')
