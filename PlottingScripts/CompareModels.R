library(dplyr)
library(data.table)
library(ggplot2)
library(glm2)
library(cowplot)

#read Genotypes
GT = fread('gzcat ~/genomes/genomes/hg19/Genotypes/CHR22.Genotypes.txt.gz',nrows=50000)

ColNames = read.table('/Users/luke/genomes/genomes/hg19/Genotypes/ColNames.txt', header=F)
names(GT) = as.character(unlist(ColNames))

samples = fread('/Users/luke/Documents/PCAperPop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt') 

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

#get the deviance of Qual
devPop = lapply(Pop, function(x) anova(x, test='Chi')[[3,2]])
devPC = lapply(PC, function(x) anova(x, test='Chi')[[7,2]])
devPopPC = lapply(PopPC, function(x) anova(x, test='Chi')[[2,2]])

plt = data.frame('coefPopPC' = unlist(coefPopPC),
				  'coefPop'=unlist(coefPop),
				  'coefPC' = unlist(coefPC),
				  'devLF' = unlist(devLF),
				  'devPop'=unlist(devPop),
				  'devPC' = unlist(devPC))
				  

plt = cbind(plt, GT[,c(1:3)])
				  
n1<-ggplot(dev, aes(x=Pop, y=PC))+geom_point()			  
n2<-ggplot(dev, aes(x=LF, y=PC))+geom_point()
n3<-ggplot(dev, aes(x=LF, y=Pop))+geom_point()
plot_grid(n1,n2,n3, nrow=1)					  

dev = cbind(coef, dev)


all = merge(dev1, coef1, by=c('CHROM','POS','ID', 'rsID'))

pLF=ggplot(dev, aes(x=as.numeric(ID), y=LFd))+geom_point(shape=1)+theme_bw()+labs(x='Allele Frequency',y='Deviance')
pPop=ggplot(dev, aes(x=as.numeric(ID), y=Popd))+geom_point(shape=1)+theme_bw()+labs(x='Allele Frequency',y='Deviance')
pb=ggplot(coef, aes(x=as.numeric(ID), y=PC))+geom_point(shape=1)+theme_bw()+labs(x='Allele Frequency',y='Beta')
pp=ggplot(coef, aes(x=as.numeric(ID), y=LF))+geom_point(shape=1)+theme_bw()+labs(x='Allele Frequency',y='Beta')
plot_grid(pb,pp,pPop,pLF)
			  
p1=ggplot(coef1, aes(x=Pop, y=PC))+geom_point(shape=1, alpha=0.3)+theme_classic()+ggtitle('Coef of Q')+theme(plot.title = element_text(hjust = 0.5))
			  
d1=ggplot(dev1, aes(x=Pop, y=PC))+geom_abline(slope=1,intercept=0, color='blue')+geom_point(shape=1, alpha=0.3)+theme_classic()+ggtitle('Deviance of Q')+theme(plot.title = element_text(hjust = 0.5))

plot_grid(p1,d1, pa,pb, nrow=2)

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


#Freq Test

mGT = melt(unique(GT), id = c('CHROM','POS','ID'))
mGT = merge(samples[,c('Name','Pop','average_quality_of_mapped_bases')], mGT, by.x='Name',by.y='variable')
mGT = merge(mGT, SigPos, by=c('CHROM','POS','ID'))

over30 = mGT[which(mGT$average_quality_of_mapped_bases > 30),]

under30 = mGT[which(mGT$average_quality_of_mapped_bases < 30),]

o30AF <- over30 %>% group_by(CHROM,POS,rsID) %>% summarize(oAF = sum(value))
u30AF <- under30 %>% group_by(CHROM,POS,rsID) %>% summarize(uAF = sum(value))

AF30 <- merge(u30AF, o30AF, by=c('CHROM','POS','rsID'))

af1 = ggplot(AF30, aes(x=uAF, y=oAF))+geom_abline(intercept=0, slope = 2279/225, linetype=2, color='grey60')+geom_bin2d(bins=100)+scale_fill_distiller(palette = "Spectral")+labs(x='Samples with Q < 30',y='Samples with Q > 30')

af2 = ggplot(AF30[which(AF30$oAF <50),], aes(x=uAF, y=oAF))+geom_abline(intercept=0, slope = 2279/225, linetype=2, color='grey60')+geom_bin2d(bins=50)+scale_fill_distiller(palette = "Spectral")+labs(x='Samples with Q < 30',y='Samples with Q > 30')

plot_grid(af1, af2)

all = merge(all, AF30, by=c('CHROM','POS','rsID'))

af3=ggplot(all, aes(x=uAF, y=oAF, color=PC.y))+geom_abline(intercept=0, slope = 2279/225, linetype=2, color='grey60')+geom_point()+scale_color_gradient(low='green',high='blue')+labs(x='Samples with Q < 30',y='Samples with Q > 30',color='Beta')

af4=ggplot(all, aes(x=uAF, y=oAF, color=PC.x))+geom_abline(intercept=0, slope = 2279/225, linetype=2, color='grey60')+geom_point()+scale_color_gradient(low='green',high='blue')+labs(x='Samples with Q < 30',y='Samples with Q > 30',color='Deviance')

plot_grid(af3, af4)