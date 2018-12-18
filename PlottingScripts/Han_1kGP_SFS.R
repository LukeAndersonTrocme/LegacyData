##plot SFS sig Han
library(data.table)
library(ggplot2)
library(reshape2)

awkPos<-fread('/Users/luke/Documents/QualityPaper/sig/Total_Sig_awk_POS.txt')

Qual<-fread('~/Dropbox/LukeTemp/SubMeta.txt')
Qual$V1<-NULL

Han<- fread('~/genomes/genomes/Han90/83_Han_Sig.frq')
Han$V1 <- as.numeric(Han$V1)
Han$V2 <- as.numeric(Han$V2)

					
KGP <- read.table('~/Documents/QualityPaper/SigVCF/1kGP_90Han_shared.frq',
					sep="\t",fill=TRUE,col.names=paste("V",1:8,sep=''), skip=1)
					
plt<-unique(merge(Han, KGP, by=c('V1','V2'), all=T))

plt[which(is.na(plt$V6.x)),]$V6.x <- 0
plt[which(is.na(plt$V6.y)),]$V6.y <- 0

ggplot(plt, aes(x=V6.x, y=V6.y))+geom_bin2d(binwidth=c(1/83,1/83))+scale_fill_distiller(palette = "Spectral")+labs(x='high depth',y='1kGP')+theme_classic()

ggsave('~/Documents/QualityPaper/Figures/Han83.jpg',height=5,width=5)

#gt1<-fread('/Users/luke/Documents/QualityPaper/SigVCF/1kGP_all.genotypes.txt')
gt1<-fread('~/Documents/QualityPaper/SigVCF/1kGP_83Han_genomeWide.genotypes.txt')
gt1.m<-melt(gt1, id.vars = c('#CHROM','POS'))
colnames(gt1.m)=c('CHROM','POS','Sample','kGT')


#gtHan<-fread('~/genomes/genomes/Han90/83_Han_all.genotypes.txt')
gtHan <-fread('~/genomes/genomes/Han90/90_Han_Chinese_SigPos.genotypes.txt')
gtHan$'#CHROM' <- as.integer(as.character(gsub('chr','',gtHan$'#CHROM')))
gtHan.m<-melt(gtHan, id.vars = c('#CHROM','POS'))
colnames(gtHan.m)=c('CHROM','POS','Sample','HanGT')

gt<-merge(gtHan.m, gt1.m, by=c('CHROM','POS','Sample'))
gt <-merge(gt, Qual, by.x='Sample',by.y='Name')
gt$diff <- (1-as.numeric(gt$HanGT == gt$kGT))
gt$Pop<-NULL

Concordance<-as.data.frame.matrix(table(gt$kGT, gt$HanGT))

AF<-aggregate(. ~ CHROM+POS, data = gt[,c('CHROM','POS','HanGT','kGT')], sum)
ggplot(AF, aes(x=HanGT/2, y=kGT/2))+geom_bin2d()+scale_fill_distiller(palette = "Spectral")+labs(x='high depth',y='1kGP')+theme_classic()


tableGT<-as.data.frame(table(gt$Sample,gt$diff))

tableGT  <- merge(tableGT, Qual, by.x='Var1',by.y='Name')
tableGT$prop <- tableGT$Freq/162*100
ggplot(tableGT[which(tableGT$Var2==0),], aes(x= average_quality_of_mapped_bases, y=prop))+geom_point()

agg.gt<-aggregate(. ~ Sample, data = gt[,c('Sample','diff')], sum)
agg.gt$diff <- agg.gt$diff/162*100

agg.pos<-aggregate(gt[,c('POS','average_quality_of_mapped_bases')], by=list(gt$HanGT,gt$POS),FUN= mean)

ggplot(agg.gt, aes(x=HanGT, y= kGT, color= average_quality_of_mapped_bases))+geom_point()+ scale_color_continuous(name = "Quality")+theme_classic()+labs(x='# of Mutations in Resequence',y='# of Mutations in 1kGP')

ggplot(agg.gt, aes(y=diff, x= average_quality_of_mapped_bases, color=Pop))+geom_point()+labs(y='% mismatches')+theme_classic()

ggplot(agg.gt, aes(y=HanGT, x= average_quality_of_mapped_bases, color=Pop))+geom_point()+theme_classic()

ggplot(agg.gt, aes(y=kGT, x= average_quality_of_mapped_bases, color=Pop))+geom_point()+theme_classic()

ggplot(agg.pos, aes(x=as.character(Group.1),y= average_quality_of_mapped_bases))+geom_boxplot()