##plot SFS sig Han
library(data.table)
library(ggplot2)
library(reshape2)

Qual<-fread('~/Dropbox/LukeTemp/SubMeta.txt')
Qual$V1<-NULL

Han<- fread('~/genomes/genomes/Han90/83_Han_Sig.frq')

KGP<- fread('~/Documents/QualityPaper/SigVCF/1kGP_90Han_shared.frq')

plt<-unique(merge(Han, KGP, by=c('V1','V2'), all=T))

ggplot(plt, aes(x=V6.x, y=V6.y))+geom_bin2d(binwidth=c(1/83,1/83))+scale_fill_distiller(palette = "Spectral")+labs(x='high depth',y='1kGP')+theme_classic()

gt1<-fread('/Users/luke/Documents/QualityPaper/SigVCF/1kGP_all.genotypes.txt')
gt1.m<-melt(gt1, id.vars = c('#CHROM','POS'))
colnames(gt1.m)=c('CHROM','POS','Sample','kGT')


gtHan<-fread('~/genomes/genomes/Han90/83_Han_all.genotypes.txt')
gtHan$'#CHROM' <- as.integer(as.character(gsub('chr','',gtHan$'#CHROM')))
gtHan.m<-melt(gtHan, id.vars = c('#CHROM','POS'))
colnames(gtHan.m)=c('CHROM','POS','Sample','HanGT')

gt<-merge(gtHan.m, gt1.m, by=c('CHROM','POS','Sample'))
gt <-merge(gt, Qual, by.x='Sample',by.y='Name')
gt$diff <- (1-as.numeric(gt$HanGT == gt$kGT))

agg.gt<-aggregate(. ~ Sample, data = gt, sum)
agg.gt$diff <- agg.gt$diff/149*100

agg.pos<-aggregate(gt[,c('POS','average_quality_of_mapped_bases')], by=list(gt$HanGT,gt$POS),FUN= mean)

ggplot(agg.gt, aes(x=HanGT, y= kGT, color= average_quality_of_mapped_bases))+geom_point()+ scale_color_continuous(name = "Quality")+theme_classic()+labs(x='# of Mutations in Resequence',y='# of Mutations in 1kGP')

ggplot(agg.gt, aes(y=diff, x= average_quality_of_mapped_bases, color=Pop))+geom_point()+labs(y='% mismatches')+theme_classic()

ggplot(agg.gt, aes(y=HanGT, x= average_quality_of_mapped_bases, color=Pop))+geom_point()+theme_classic()

ggplot(agg.gt, aes(y=kGT, x= average_quality_of_mapped_bases, color=Pop))+geom_point()+theme_classic()

ggplot(agg.pos, aes(x=as.character(Group.1),y= average_quality_of_mapped_bases))+geom_boxp()