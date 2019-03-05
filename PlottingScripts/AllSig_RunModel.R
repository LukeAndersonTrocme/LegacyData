library(dplyr)
library(data.table)
library(ggplot2)
library(glm2)
library(cowplot)

#read colors per pop (used in 1kG Nature)
colors<-read.table("~/Dropbox/LukeTemp/NamePopBigPop_color.html.txt", comment.char='>',col.names=c('Pop','hex','BigPop'))
MyColour <- as.character(colors $hex)
names(MyColour) <- colors $Pop

path = '/Users/luke/genomes/genomes/hg19/Genotypes/'

#read Genotypes
GT = fread('~/Documents/QualityPaper/sig/SigVar_Pop_0.01_genotype.txt')
ColNames = read.table('/Users/luke/genomes/genomes/hg19/Genotypes/ColNames.txt', header=F)
names(GT) = as.character(unlist(ColNames))

samples = fread('/Users/luke/Documents/PCAperPop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt')
samples$Pop = factor(samples$Pop, levels=c('FIN','GBR','CEU','IBS','TSI','CHS','CDX','CHB','JPT','KHV','GIH','STU','PJL','ITU','BEB','PEL','MXL','CLM','PUR','ASW','ACB','GWD','YRI','LWK','ESN','MSL'))

model = apply(GT[,-(1:3)], 1, function(x)
						glm2(x ~
							samples$Pop +
							samples$average_quality_of_mapped_bases))
#get coef of Qual per site
coef = lapply(model, function(x) coef(x)[[27]])				

dev = lapply(model, function(x) anova(x, test='Chi')[[3,2]])

plt = cbind(unlist(coef), GT[,c(1:3)])
plt = cbind(plt, unlist(dev))
plt$ID = as.numeric(as.character(plt$ID))
plt$pval = -log10(pchisq(plt$V2, 1, lower.tail=F))

#to deal with duplicate entries, merge based on Chr Pos and AF
SigPos = fread('/Users/luke/Documents/QualityPaper/sig/SignificantVariants.txt')
plt = merge(plt, SigPos, by = c('CHROM','POS','ID'))

names(plt) = c('Beta','CHROM','POS','AF','deviance','-log10(p)')

write.table(plt,'/Users/luke/Documents/QualityPaper/sig/SigVar_1kGP_Coef.txt', quote=F, row.names=F)

ggplot(plt, aes(x=ID, y=V1, color=log10P_0.01))+geom_point()+theme_bw()+labs(x='Allele Frequency',y='Beta')

##Compare over and under 30

mGT = melt(unique(GT), id = c('CHROM','POS','ID'))
mGT = merge(samples[,c('Name','Pop','average_quality_of_mapped_bases')], mGT, by.x='Name',by.y='variable')
mGT = merge(mGT, SigPos, by=c('CHROM','POS','ID'))

over30 = mGT[which(mGT$average_quality_of_mapped_bases > 30),]

under30 = mGT[which(mGT$average_quality_of_mapped_bases < 30),]

o30AF <- over30 %>% group_by(CHROM,POS,ID) %>% summarize(oAF = sum(value))
u30AF <- under30 %>% group_by(CHROM,POS,ID) %>% summarize(uAF = sum(value))

AF30 <- merge(u30AF, o30AF, by=c('CHROM','POS','ID'))

af1 = ggplot(AF30, aes(x=uAF, y=oAF))+geom_abline(intercept=0, slope = 2279/225, linetype=2, color='grey60')+geom_bin2d(bins=100)+scale_fill_distiller(palette = "Spectral")+labs(x='Samples with Q < 30',y='Samples with Q > 30')

af2 = ggplot(AF30[which(AF30$oAF <50),], aes(x=uAF, y=oAF))+geom_abline(intercept=0, slope = 2279/225, linetype=2, color='grey60')+geom_bin2d(bins=50)+scale_fill_distiller(palette = "Spectral")+labs(x='Samples with Q < 30',y='Samples with Q > 30')

plot_grid(af1, af2)

all = merge(all, AF30, by=c('CHROM','POS','rsID'))

af3=ggplot(all, aes(x=uAF, y=oAF, color=PC.y))+geom_abline(intercept=0, slope = 2279/225, linetype=2, color='grey60')+geom_point()+scale_color_gradient(low='green',high='blue')+labs(x='Samples with Q < 30',y='Samples with Q > 30',color='Beta')

af4=ggplot(all, aes(x=uAF, y=oAF, color=PC.x))+geom_abline(intercept=0, slope = 2279/225, linetype=2, color='grey60')+geom_point()+scale_color_gradient(low='green',high='blue')+labs(x='Samples with Q < 30',y='Samples with Q > 30',color='Deviance')

plot_grid(af3, af4)

