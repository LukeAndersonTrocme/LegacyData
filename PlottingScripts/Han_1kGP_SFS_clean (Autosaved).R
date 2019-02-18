##plot SFS sig Han
library(data.table)
library(ggplot2)
library(reshape2)

#read Genotypes
kgp1 = fread('~/Documents/QualityPaper/sig/SigInHan_GenomeWide_gt.txt')
ColNames = read.table('/Users/luke/genomes/genomes/hg19/Genotypes/ColNames.txt', header=F)
names(kgp1) = as.character(unlist(ColNames))
#melt to long format
kgp = melt(kgp1, id=c('CHROM','POS','ID'))

#read Genotypes
han1=fread('~/Documents/QualityPaper/sig/Han90_SigInHan_GenomeWide_gt.txt')
HanNames=fread('~/genomes/genomes/Han90/SignificantSNPs_Han90_COLNAMES.txt',header=F)
names(han1) = as.character(unlist(HanNames))
#melt to long format
han = melt(han1, id=c('#CHROM','POS','ID'))
#fix chromosome names
han$CHROM = as.numeric(as.character(substring(han$'#CHROM',4)))
han$'#CHROM' <- NULL

#merge 1kGP and han resequenced data.
combo <- merge(han, kgp, by=c('CHROM','POS','variable'), all=F)
#to deal with INDELs we only count one instance of duplicate positions (i.e. pick max GT)
combo <- combo %>% group_by(CHROM,POS,variable) %>% summarize(Han.gt = max(value.x), Kgp.gt = max(value.y))
#to get allele frequencies we group by POS and summarize GT
agg <- combo %>% group_by(CHROM,POS) %>% summarize(hanAF = sum(Han.gt), kgpAF = sum(Kgp.gt)) 


ggplot(agg, aes(x=hanAF, y=kgpAF))+geom_bin2d()+scale_fill_distiller(palette = "Spectral")+labs(x='Resequenced Han', y='1000 Genomes Project')+ggtitle('Joint Frequency Plot of Suspicious variants')
ggsave('~/Documents/QualityPaper/Figures/Han_1kGP_SFS.jpg',height=10, width=10)

#concordance per sample
samples = fread('/Users/luke/Documents/PCAperPop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt')
#get quality per sample
combo = merge(combo, samples[,c('Name','average_quality_of_mapped_bases')], by.x='variable',by.y='Name')
#get difference between Han and kGP for each SNP
combo$concordance = as.factor(combo$Kgp.gt - combo$Han.gt)
levels(combo$concordance) <- c('Missing in kGP','Match','Missing in Han')

con <- combo %>% group_by(variable, average_quality_of_mapped_bases, concordance) %>% summarize(Count = n())
ggplot(con, aes(x= average_quality_of_mapped_bases, y = Count))+geom_point()+facet_wrap(.~concordance, scales='free')
combo$concor = combo$Han.gt - combo$Kgp.gt
con2 <- combo %>% group_by(variable, average_quality_of_mapped_bases) %>% summarize(concor = sum(concor))
ggplot(con2, aes(x= average_quality_of_mapped_bases, y = concor))+geom_point()

ggsave('~/Documents/QualityPaper/Figures/Concordance.jpg',height=4, width=12)










