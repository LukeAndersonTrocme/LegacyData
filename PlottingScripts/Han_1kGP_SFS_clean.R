##plot SFS sig Han
library(data.table)
library(ggplot2)
library(reshape2)

#read Genotypes
han1=fread('~/Documents/QualityPaper/sig/Han90_SigInHan_GenomeWide_gt.txt')
HanNames=fread('~/genomes/genomes/Han90/SignificantSNPs_Han90_COLNAMES.txt',header=F)
names(han1) = as.character(unlist(HanNames))
#melt to long format
han = melt(han1, id=c('#CHROM','POS','ID'))
#fix chromosome names
han$CHROM = as.numeric(as.character(substring(han$'#CHROM',4)))
han$'#CHROM' <- NULL

#read Genotypes
kgp1 = fread('/Users/luke/Documents/QualityPaper/sig/SigVar_1kGP_clean_GT.txt')
ColNames = read.table('/Users/luke/genomes/genomes/hg19/Genotypes/ColNames.txt', header=F)
names(kgp1) = as.character(unlist(ColNames))
#melt to long format
kgp = melt(kgp1, id=c('CHROM','POS','ID'))
#keep only the 83 resequenced individuals
kgp = kgp[which(kgp$variable %in% han$variable),]

#merge 1kGP and han resequenced data.
combo <- merge(han, kgp, by=c('CHROM','POS','variable'), all=T)
#missing sites replaced with 0
combo[is.na(combo)] <- 0 
#to deal with INDELs we only count one instance of duplicate positions (i.e. pick max GT)
combo <- combo %>% group_by(CHROM,POS,variable) %>% summarize(Han.gt = max(value.x), Kgp.gt = max(value.y))

#to get allele frequencies we group by POS and summarize GT
agg <- combo %>% group_by(CHROM,POS) %>% summarize(hanAF = sum(Han.gt), kgpAF = sum(Kgp.gt)) 
#remove sites missing in both groups
agg <- agg[which(agg$hanAF + agg$kgpAF > 1),]

INDEL = merge(Indels, agg, by=c('CHROM','POS'))
SNP = agg[which((agg$POS %in% INDEL$POS)==F),]

p1 =ggplot(SNP, aes(x=hanAF, y=kgpAF))+geom_bin2d(bins = 83)+scale_fill_distiller(palette = "Spectral")+labs(x='Resequenced Han', y='1000 Genomes Project')+ggtitle('Joint Frequency Plot of Suspicious SNPs')

p2 =ggplot(INDEL, aes(x=hanAF, y=kgpAF))+geom_bin2d(bins = 83)+scale_fill_distiller(palette = "Spectral")+labs(x='Resequenced Han', y='1000 Genomes Project')+ggtitle('Joint Frequency Plot of Suspicious INDELS')

plot_grid(p1,p2)
ggsave('~/Documents/QualityPaper/Figures/Han_1kGP_SFS_noSignles.jpg',height=10, width=20)


#concordance per sample
samples = fread('/Users/luke/Documents/PCAperPop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt')
#get quality per sample
combo = merge(combo, samples[,c('Name','average_quality_of_mapped_bases')], by.x='variable',by.y='Name')
#get difference between Han and kGP for each SNP
combo$concordance = as.factor(combo$Kgp.gt - combo$Han.gt)
levels(combo$concordance) <- c('Missing in kGP','Match','Missing in Han')

con <- combo %>% group_by(variable, average_quality_of_mapped_bases, concordance) %>% summarize(Count = n())
ggplot(con, aes(x= average_quality_of_mapped_bases, y = Count))+geom_point()+facet_wrap(.~concordance, scales='free')
combo$concor = combo$Kgp.gt - combo$Han.gt
con2 <- combo %>% group_by(variable, average_quality_of_mapped_bases) %>% summarize(concor = sum(concor))
ggplot(con2, aes(x= average_quality_of_mapped_bases, y = concor))+geom_point()

ggsave('~/Documents/QualityPaper/Figures/Concordance.jpg',height=4, width=12)










