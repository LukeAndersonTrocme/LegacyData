library(dplyr)
library(data.table)
library(ggplot2)
library(glm2)
library(cowplot)

#read Genotypes
GT = fread('gzcat ~/Desktop/CHR22.Genotypes.txt.gz',nrows=10000)
GT = fread('~/Desktop/SigVar_1kGP_clean_GT.txt')
SigPos = fread('~/Desktop/Total_Sig_POSfreq.txt', col.names=c('CHROM','POS','rsID','ID'))

ColNames = read.table('~/Desktop/ColNames.txt', header=F)
names(GT) = as.character(unlist(ColNames))
samples = fread('~/Desktop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt') 
#/Users/luke/Documents/PCAperPop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt

Pop = apply(GT[,-(1:3)], 1, function(x)
						glm2(x ~
							samples$Pop +
							samples$average_quality_of_mapped_bases))

PC = apply(GT[,-(1:3)], 1, function(x)
						glm2(x ~
							samples$PC1 + 
							samples$PC2 +
							samples$PC3 +
							samples$PC4 +
							samples$PC5 + 
							samples$average_quality_of_mapped_bases))
							
PopPC = apply(GT[,-(1:3)], 1, function(x)
						glm2(x ~
							samples$Pop +
							samples$PC1 + 
							samples$PC2 +
							samples$PC3 +
							samples$PC4 +
							samples$PC5 + 
							samples$average_quality_of_mapped_bases))														
#get coef of Qual per site
coefPop = lapply(Pop, function(x) coef(x)[[27]])
coefPC = lapply(PC, function(x) coef(x)[[7]])
coefPopPC = lapply(PopPC, function(x) coef(x)[[32]])

coef = data.frame('Pop'=unlist(coefPop),
				  'PC' = unlist(coefPC),
				  'PopPC' = unlist(coefPopPC))

coef = cbind(coef, GT[,c(1:3)])
coef1 = merge(coef, SigPos, by=c('CHROM','POS','ID'))
pa=ggplot(coef1, aes(x=as.numeric(ID), y=Pop))+geom_point(shape=1)+theme_bw()+labs(x='Allele Frequency',y='Beta')
pb=ggplot(coef1, aes(x=as.numeric(ID), y=PC))+geom_point(shape=1)+theme_bw()+labs(x='Allele Frequency',y='Beta')
pc=ggplot(coef1, aes(x=as.numeric(ID), y=PopPC))+geom_point(shape=1)+theme_bw()+labs(x='Allele Frequency',y='Beta')
			  
p1=ggplot(coef1, aes(x=Pop, y=PC))+geom_point(shape=1, alpha=0.3)+theme_classic()+ggtitle('Coef of Q')+theme(plot.title = element_text(hjust = 0.5))
p2=ggplot(coef1, aes(x=Pop, y=PopPC))+geom_point(shape=1, alpha=0.3)+theme_classic()+ggtitle('Coef of Q')+theme(plot.title = element_text(hjust = 0.5))		  
p3=ggplot(coef1, aes(x=PC, y=PopPC))+geom_point(shape=1, alpha=0.3)+theme_classic()+ggtitle('Coef of Q')+theme(plot.title = element_text(hjust = 0.5))
plot_grid(p2,p3,p1, pa,pb,pc, nrow=2)		

devPop = lapply(Pop, function(x) anova(x, test='Chi')[[3,2]])
devPC = lapply(PC, function(x) anova(x, test='Chi')[[7,2]])
devPopPC = lapply(PopPC, function(x) anova(x, test='Chi')[[8,2]])

dev = data.frame('Pop'=unlist(devPop),
				  'PC' = unlist(devPC),
				  'PopPC' = unlist(devPopPC))
				  
dev1 = cbind(dev, GT[,c(1:3)])
dev1 = merge(dev1, SigPos, by=c('CHROM','POS','ID'))				  
d1=ggplot(dev1, aes(x=Pop, y=PC))+geom_abline(slope=1,intercept=0, color='blue')+geom_point(shape=1, alpha=0.3)+theme_classic()+ggtitle('Deviance of Q')+theme(plot.title = element_text(hjust = 0.5))
d2=ggplot(dev1, aes(x=Pop, y=PopPC))+geom_abline(slope=1,intercept=0, color='blue')+geom_point(shape=1, alpha=0.3)+theme_classic()+ggtitle('Deviance of Q')+theme(plot.title = element_text(hjust = 0.5))		  
d3=ggplot(dev1, aes(x=PC, y=PopPC))+geom_abline(slope=1,intercept=0, color='blue')+geom_point(shape=1, alpha=0.3)+theme_classic()+ggtitle('Deviance of Q')+theme(plot.title = element_text(hjust = 0.5))
plot_grid(p2,p3,p1, d2,d3,d1, pa,pb,pc, nrow=3)

ggplot(dev1, aes(x=Pop, y=PC))+geom_abline(slope=1,intercept=0, color='blue')+geom_point(shape=1, alpha=0.3)+theme_classic()+ggtitle('Deviance of Q')+theme(plot.title = element_text(hjust = 0.5))+xlim(c(0,10))+ylim(c(0,10))

q1=ggplot(dev1, aes(x=as.numeric(ID), y=PC))+geom_point(shape=1, alpha=0.3)+theme_classic()+ggtitle('Deviance of Q')+theme(plot.title = element_text(hjust = 0.5))

q2=ggplot(dev1, aes(x=as.numeric(ID), y=Pop))+geom_point(shape=1, alpha=0.3)+theme_classic()+ggtitle('Deviance of Q')+theme(plot.title = element_text(hjust = 0.5))

plot_grid(q1,q2)



coef = cbind(unlist(coef), GT[,c(1:3)])
coef = cbind(coef, unlist(dev))
coef$ID = as.numeric(as.character(coef$ID))
coef$pval = -log10(pchisq(coef$V2, 1, lower.tail=F))
names(coef) = c('Beta','CHROM','POS','AF','deviance','-log10(p)')
write.table(coef,'/Users/luke/Documents/QualityPaper/sig/SigVar_1kGP_Coef.txt', quote=F, row.names=F)

