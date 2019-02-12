library(dplyr)
library(data.table)
library(ggplot2)
library(glm2)

path = '/Users/luke/genomes/genomes/hg19/Genotypes/'

GT = fread('~/Desktop/PublishedHits.txt')
ColNames = read.table(paste(path,"ColNames.txt", sep=''), header=F)
names(GT)= as.character(unlist(ColNames))
samples = fread('/Users/luke/Documents/PCAperPop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt')

logP = apply(GT[,-(1:3)], 1, function(x) 
				anova(
					glm2(x ~
							samples$PC1 + 
							samples$PC2 +
							samples$PC3 +
							samples$PC4 +
							samples$PC5 + 
							samples$average_quality_of_mapped_bases,
						family=binomial),
					test="Chi")[7,2])

model = apply(GT[,-(1:3)], 1, function(x)
					glm2(x ~
							samples$PC1 + 
							samples$PC2 +
							samples$PC3 +
							samples$PC4 +
							samples$PC5 + 
							samples$average_quality_of_mapped_bases,
						family=binomial))
						
exp.coef = apply(GT[,-(1:3)], 1, function(x)
					exp(coefficients(glm2(x ~
							samples$PC1 + 
							samples$PC2 +
							samples$PC3 +
							samples$PC4 +
							samples$PC5 + 
							samples$average_quality_of_mapped_bases,
						family=binomial)))[[7]])						

dat = cbind(GT[,c(1,2)], exp.coef)
rsid = fread('~/Desktop/SigPos.txt')
rsID <- merge(dat, rsid, by.x=c("CHROM",'POS'), by.y=c('rsID','CHROM'))

mGT = melt(GT,id=c('CHROM','POS','ID'))
mGT <- merge(samples, mGT, by.y='variable',by.x='Name')
cGT <- mGT %>% group_by(Pop) %>% summarize(var = sum(value))
keep <- cGT[which(cGT$var > 0),]$Pop