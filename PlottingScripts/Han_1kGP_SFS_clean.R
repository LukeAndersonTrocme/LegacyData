##plot SFS sig Han
library(data.table)
library(ggplot2)
library(reshape2)

#read Genotypes
kgp1='~/Documents/QualityPaper/SigVCF/1kGP_83Han_genomeWide.genotypes.txt'
kgp1.n=fread('~/Documents/QualityPaper/SigVCF/1kGP_83Han_genomeWide.header')

han1='~/genomes/genomes/Han90/90_Han_Chinese_SigPos.genotypes.txt'
han1.n=fread('~/genomes/genomes/Han90/90_Han_Chinese_SigPos.header')

gt1<-unique(fread(kgp1,col.names=colnames(kgp1.n)))
gt1.m<-melt(gt1, id.vars = c('#CHROM','POS','ID'))
colnames(gt1.m)=c('CHROM','POS','ID','Sample','kGT')

gtHan <-unique(fread(han1,col.names=colnames(han1.n)))
gtHan$'#CHROM' <- as.integer(as.character(gsub('chr','',gtHan$'#CHROM')))
gtHan.m<-melt(gtHan, id.vars = c('#CHROM','POS','ID'))
colnames(gtHan.m)=c('CHROM','POS','ID','Sample','HanGT')

#merge
gt<-merge(gtHan.m, gt1.m, by=c('CHROM','POS','ID','Sample'), all=T)
gt <-merge(gt, Qual, by.x='Sample',by.y='Name')
gt[is.na(gt)] <- 0

Concordance<-as.data.frame.matrix(table(gt$kGT, gt$HanGT))

AF<-aggregate(. ~ POS+ID, data = gt[,c('POS','ID','HanGT','kGT')], sum)
AF<-AF[-which((AF$HanGT==0)&(AF$kGT==0)),]
ggplot(AF, aes(x=HanGT, y=kGT))+geom_bin2d(bins=83)+scale_fill_distiller(palette = "Spectral")+labs(x='high depth',y='1kGP')+theme_classic()

ggsave('~/Documents/QualityPaper/Figures/Han83.jpg',height=5,width=6)

kgp2='~/Documents/QualityPaper/SigVCF/1kGP_genomeWide_SINGLEPOP_83.genotypes.txt'
kgp2.n=fread('~/Documents/QualityPaper/SigVCF/1kGP_genomeWide_SINGLEPOP_83.header')

han2='~/genomes/genomes/Han90/90_Han_Chinese_SinglePop.genotypes.txt'
han2.n=fread('~/genomes/genomes/Han90/90_Han_Chinese_SinglePop.header')

gt1<-unique(fread(kgp2,col.names=colnames(kgp2.n)))
gt1.m<-melt(gt1, id.vars = c('#CHROM','POS','ID'))
colnames(gt1.m)=c('CHROM','POS','ID','Sample','kGT')


gtHan <-unique(fread(han2,col.names=colnames(han2.n)))
gtHan$'#CHROM' <- as.integer(as.character(gsub('chr','',gtHan$'#CHROM')))
gtHan.m<-melt(gtHan, id.vars = c('#CHROM','POS','ID'))
colnames(gtHan.m)=c('CHROM','POS','ID','Sample','HanGT')

#merge
gt<-merge(gtHan.m, gt1.m, by=c('CHROM','POS','ID','Sample'), all=T)
gt <-merge(gt, Qual, by.x='Sample',by.y='Name')
gt[is.na(gt)] <- 0

Concordance<-as.data.frame.matrix(table(gt$kGT, gt$HanGT))

AF<-aggregate(. ~ POS+ID, data = gt[,c('POS','ID','HanGT','kGT')], sum)
AF<-AF[-which((AF$HanGT==0)&(AF$kGT==0)),]
ggplot(AF, aes(x=HanGT, y=kGT))+geom_bin2d(bins=83)+scale_fill_distiller(palette = "Spectral")+labs(x='high depth',y='1kGP')+theme_classic()
ggsave('~/Documents/QualityPaper/Figures/Han83_singlePop.jpg',height=5,width=6)