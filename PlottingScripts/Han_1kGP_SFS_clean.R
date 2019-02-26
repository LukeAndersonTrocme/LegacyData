##plot SFS sig Han
library(data.table)
library(ggplot2)
library(reshape2)

#bcftools view -R /Users/luke/Documents/QualityPaper/sig/SignificantVariants_chr_CHROM_POS.txt /Users/luke/genomes/genomes/Han90/90_Han_Chinese_GenomeWide_sorted.vcf.gz  -Oz -o /Users/luke/genomes/genomes/Han90/SignificantSNPs_Han90_0.001.vcf.gz

#java -jar ~/bin/snpEff_latest_core/snpEff/SnpSift.jar extractFields /Users/luke/genomes/genomes/Han90/SignificantSNPs_Han90_0.001.vcf.gz CHROM POS AF "GEN[*].GT" | sed "s/0\/0/0/g ; s/0\/[1-9]/1/g ; s/[1-9]\/0/1/g ; s/[1-9]\/[1-9]/1/g ; /\//d" | tail -n +2 > /Users/luke/genomes/genomes/Han90/SignificantSNPs_Han90_0.001_genotype.txt


#read Genotypes
han1=fread('/Users/luke/genomes/genomes/Han90/SignificantSNPs_Han90_0.001_genotype.txt')
HanNames=fread('~/genomes/genomes/Han90/SignificantSNPs_Han90_COLNAMES.txt',header=F)
names(han1) = as.character(unlist(HanNames))
#melt to long format
han = melt(han1, id=c('#CHROM','POS','ID'))
#fix chromosome names
han$CHROM = as.numeric(as.character(substring(han$'#CHROM',4)))
han$'#CHROM' <- NULL

#to deal with INDELs we only count one instance of duplicate positions (i.e. pick max GT)
han <- han %>% group_by(CHROM,POS,variable) %>% summarize(han.gt = max(value))
#to get allele frequencies we group by POS and summarize GT
hanAF <- han %>% group_by(CHROM,POS) %>% summarize(hanAF = sum(han.gt))


#read Genotypes
kgp1 = fread('~/Documents/QualityPaper/sig/SigVar_0.001.genotypes.txt')
ColNames = read.table('/Users/luke/genomes/genomes/hg19/Genotypes/ColNames.txt', header=F)
names(kgp1) = as.character(unlist(ColNames))
#melt to long format
kgp = melt(kgp1, id=c('CHROM','POS','ID'))
#keep only the 83 resequenced individuals
kgp = kgp[which(kgp$variable %in% han$variable),]

#to deal with INDELs we only count one instance of duplicate positions (i.e. pick max GT)
kgp <- kgp %>% group_by(CHROM,POS,variable) %>% summarize(Kgp.gt = max(value))
#to get allele frequencies we group by POS and summarize GT
kgpAF <- kgp %>% group_by(CHROM,POS) %>% summarize(kgpAF = sum(Kgp.gt))
#keep only positions that are present in the 83 individuals
kgpAF <- kgpAF[which(kgpAF$kgpAF > 0),]


#merge 1kGP and han resequenced data.
combo <- merge(hanAF, kgpAF, by=c('CHROM','POS'), all.y=T)
#missing sites replaced with 0
combo[is.na(combo)] <- 0 

##read Indels
Indels <- fread('/Users/luke/Documents/QualityPaper/Misc/Reg_indels.POS', skip=1 ,col.names=c('CHROM','POS'))
Indels$Indel = 'yes'

INDEL = merge(Indels, combo, by=c('CHROM','POS'))
SNP = combo[which((combo$POS %in% INDEL$POS)==F),]

p1 =ggplot(SNP, aes(x=hanAF, y=kgpAF))+geom_bin2d(bins = 83)+scale_fill_distiller(palette = "Spectral")+labs(x='Resequenced Han', y='1000 Genomes Project')+ggtitle('Joint Frequency Plot of Suspicious SNPs')

p2 =ggplot(INDEL, aes(x=hanAF, y=kgpAF))+geom_bin2d(bins = 83)+scale_fill_distiller(palette = "Spectral")+labs(x='Resequenced Han', y='1000 Genomes Project')+ggtitle('Joint Frequency Plot of Suspicious INDELS')

plot_grid(p1,p2)
ggsave('~/Documents/QualityPaper/Figures/Han_1kGP_SFS.jpg',height=10, width=20)


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


####SINGLE POPULATION TEST

kgp1='~/Documents/QualityPaper/SigVCF/1kGP_genomeWide_SINGLEPOP_83.genotypes.txt'
kgp1.n=fread('~/Documents/QualityPaper/SigVCF/1kGP_genomeWide_SINGLEPOP_83.header')

han1='~/genomes/genomes/Han90/90_Han_Chinese_SinglePop.genotypes.txt'
han1.n=fread('~/genomes/genomes/Han90/90_Han_Chinese_SINGLEPOP.header')

gt1<-unique(fread(kgp1,col.names=colnames(kgp1.n)))
gt1.m<-melt(gt1, id.vars = c('#CHROM','POS','ID'))
colnames(gt1.m)=c('CHROM','POS','ID','Sample','kGT')

gtHan <-unique(fread(han1,col.names=colnames(han1.n)))
gtHan$'#CHROM' <- as.integer(as.character(gsub('chr','',gtHan$'#CHROM')))
gtHan.m<-melt(gtHan, id.vars = c('#CHROM','POS','ID'))
colnames(gtHan.m)=c('CHROM','POS','ID','Sample','HanGT')

#merge
gt<-merge(gtHan.m, gt1.m, by=c('CHROM','POS','ID','Sample'), all.y=T)
gt <-merge(gt, Qual, by.x='Sample',by.y='Name')
gt[is.na(gt)] <- 0
head(gt)

#to get allele frequencies we group by POS and summarize GT
aggAF <- gt %>% group_by(CHROM,POS) %>% summarize(hanAF = sum(HanGT), kgpAF = sum(kGT)) 
#remove sites missing in this subset of individuals
aggAF <- aggAF[which(aggAF$kgpAF > 0),]

ggplot(aggAF, aes(x=hanAF/2, y=kgpAF/2))+geom_bin2d(bins = 83)+scale_fill_distiller(palette = "Spectral")+labs(x='Resequenced Han', y='1000 Genomes Project')+ggtitle('Joint frequency plot of 296 suspicious SNPs\nfrom single population test')

ggsave('~/Documents/QualityPaper/Figures/Han_1kGP_SFS_singlePop.jpg',height=10, width=10)