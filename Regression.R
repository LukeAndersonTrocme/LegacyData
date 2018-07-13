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
#get sites, make vcf
#while read -a line; do gzcat /Users/luke/genomes/genomes/1kGP_4bedFiltered/1kGP_Chr${line[0]}.4bed_filtered.vcf.gz | awk -v rsID=${line[1]} '$3 == rsID {print $0}'; done < /Users/luke/Documents/GWAS_Qual/Final_GWAS/1kGP_GenomeWide_JPT_INT_ci_sig6.POS > /Users/luke/Documents/GWAS_Qual/Final_GWAS/JPT_sig6.vcf
#get genotypes from vcf
#java -jar ~/bin/snpEff_latest_core/snpEff/SnpSift.jar extractFields /Users/luke/Documents/GWAS_Qual/Final_GWAS/JPT_sig6_h.vcf CHROM POS "GEN[*].GT"   | sed 's/\//|/g ; s/0|0/0/g ; s/0|1/1/g ; s/1|0/1/g ; s/1|1/2/g ; /|/d' | tail -n +2 > /Users/luke/Documents/GWAS_Qual/Final_GWAS/JPT_sig6.Genotypes.txt


GT<-read.table('/Users/luke/Documents/GWAS_Qual/Final_GWAS/JPT_sig6.Genotypes.txt',header=F)
#GT<-read.table('~/Documents/QualityPaper/PublishedSNPs.Genotypes.txt',header=F)
names(GT)=c('Chr','Pos', as.character(NamePop$Name))
#cut=c('Chr','Pos',as.character(NamePop[which(NamePop$Pop == pop),]$Name))
#GT<-GT[cut]

melt.GT<-melt(GT,id=c('Chr','Pos'))
melt.GT$value<-gsub('2','0',melt.GT$value)
melt.GT$value<-as.factor(melt.GT$value)
GT.Qual<-merge(melt.GT, SubMeta, by.x='variable',by.y='Name')

#ggplot(GT.Qual, aes(x= average_quality_of_mapped_bases, y=value, color=Pop))+geom_point()+facet_grid(Chr~Pop)+scale_color_manual(values=MyColour)+guides(color=F)+theme_classic()

GT.Qual$p=''
for(pop in unique(GT.Qual$Pop)){
gtQual<-GT.Qual[which(GT.Qual$Pop==pop),]
gtQual $t=''
print(pop)
for(snp in unique(gtQual$Pos)){
	#print(snp)
	sub= gtQual[which(gtQual $Pos==snp),]
	if(length(unique(sub$value))!=1){
	sub.glm=glm(formula=value ~ average_quality_of_mapped_bases, data=sub, family=binomial)
	t1=paste(snp,"\nP =",signif(summary(sub.glm)$coefficients[2,4], 3))
	gtQual[which(gtQual $Pos==snp),]$t=t1
	#GT.Qual[which((GT.Qual $Pos==snp)& GT.Qual $Pop==pop),]$p=-log10(summary(sub.glm)$coefficients[2,4])
}
}
}
#ggplot(gtQual, aes(x= average_quality_of_mapped_bases, y=value))+geom_point(shape=3)+facet_wrap(Pos~., ncol=3)+theme_classic()
#ggplot(gtQual, aes(x= average_quality_of_mapped_bases, y=as.numeric(value)))+geom_point()+facet_wrap(t~., ncol=2)+theme_classic()+geom_smooth(method = "glm",     method.args = list(family = "binomial"), se = FALSE,size=0.5) +labs(title=paste("Genotype vs Quality for",pop), y='Genotype',x='Average quality of mapped bases')+theme(plot.title = element_text(hjust= 0.5))+scale_y_continuous(breaks=c(0,1,2))
#ggsave(paste('~/Documents/Regression/PublishedSNPs_',pop,'.jpg', sep=''),height=10, width=4)
#}



#write.table(GT.Qual, file='~/Documents/Regression/JPT_sig6_Regression.txt')
#GT.Qual <-read.table('~/Documents/Regression/JPT_sig6_Regression.txt')

##Read JPT sig 6 P vals
GwasP<-read.table('/Users/luke/Documents/GWAS_Qual/Final_GWAS/1kGP_GenomeWide_JPT_INT_ci_sig6.assoc.linear',col.names=c('Chr','SNP','Pos','A1','TEST','NMISS','BETA','SE','L95','U95','STAT','P'))
GwasP$Plog10<--log10(GwasP$P)
GwasP<-GwasP[,c('Chr','Pos','Plog10')]

G.vs.R<-merge(GwasP, GT.Qual[which(GT.Qual$Pop=='JPT'),], by=c('Chr','Pos'))
ggplot(G.vs.R, aes(x=Plog10, y=as.numeric(p)))+geom_point()+theme_classic()+labs(x='GWAS p-value',y='Regression p-value')

ggplot(GT.Qual, aes(x=reorder(Pos, Chr),y=p,  color=Pop))+geom_point()+facet_wrap(.~Pop,ncol=1)+theme(axis.text=element_blank(), axis.ticks=element_blank())

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

ggplot(plt, aes(x=ordered, y= as.numeric(p), color=ordered, group=1))+facet_wrap(.~ POS,ncol=1)+geom_hline(yintercept=6, color='grey60',linetype=3)+geom_line()+geom_point(shape=15)+scale_color_manual(values=MyColour)+theme_classic()+theme(strip.background=element_blank(), strip.text=element_blank(), axis.title=element_blank(), axis.ticks.y=element_blank(),axis.line=element_blank(), axis.text.x=element_text(angle=90, size=5, vjust=-0.05))+guides(color=F)+geom_point(data=plt[complete.cases(plt),], shape=1, size=4, color='black')

ggsave('~/Documents/Regression/GWAS_p_vs_GLM_p.jpg',height=10,width=4)