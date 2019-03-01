library(dplyr)
library(data.table)
library(ggplot2)
library(glm2)
library(cowplot)
library(lfa)
library(tidyverse)

#read Genotypes
#GT = fread('gzcat ~/genomes/genomes/hg19/Genotypes/CHR22.Genotypes.txt.gz',nrows=50000)
GT = fread('/Users/luke/Documents/QualityPaper/sig/SigVar_0.001.genotypes.txt')


SigPos = fread('~/Documents/QualityPaper/sig/Total_Sig_POSfreq.txt', col.names=c('CHROM','POS','rsID','ID'))

ColNames = read.table('/Users/luke/genomes/genomes/hg19/Genotypes/ColNames.txt', header=F)
names(GT) = as.character(unlist(ColNames))
GT = merge(GT, SigPos, by=c('CHROM','POS','ID'))
GT$rsID = NULL


samples = fread('/Users/luke/Documents/PCAperPop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt') 

# Logistic factor analysis
GT_mtx <- as.matrix(GT[,-(1:3)])
#compute the first three factors including the intercept
LF <- lfa(GT_mtx, 4)
subset <- af(GT_mtx, LF)
subset <- cbind(GT[,c(1:3)], subset)

to_datalist <- function(i, data_genotype, data_iaf){
  
    data <- tibble(genotype=data_genotype[i,],
               iaf= subset[i,])
               
    return(data)       
}

loci <- seq(1,nrow(GT_mtx))
names(loci) <- loci

data_list <- map(loci, 
                 to_datalist,
                 data_genotype = GT_mtx,
                 data_iaf = subset)

fit_model_1 <- function(data){
	
	fit <- glm2(genotype ~ samples$average_quality_of_mapped_bases + offset(iaf), 
				family ="binomial", data=data)
	
	return(fit)
}

fits <- lapply(head(data_list,10000), fit_model_1)

Pop = apply(GT[c(1:10000),-(1:3)], 1, function(x)
						glm2(x ~
							samples$Pop +
							samples$average_quality_of_mapped_bases))

PC = apply(GT[c(1:10000),-(1:3)], 1, function(x)
						glm2(x ~
							samples$PC1 + 
							samples$PC2 +
							samples$PC3 +
							samples$PC4 +
							samples$PC5 + 
							samples$average_quality_of_mapped_bases))
																					
#get coef of Qual per site
coefLF = lapply(fits, function(x) coef(x)[[2]])
coefPop = lapply(Pop, function(x) coef(x)[[27]])
coefPC = lapply(PC, function(x) coef(x)[[7]])

coef = data.frame('LF' = unlist(coefLF),
				  'Pop'=unlist(coefPop),
				  'PC' = unlist(coefPC))
				  
q0<-ggplot(coef, aes(x=Pop, y=PC))+geom_point()
q1<-ggplot(coef, aes(x=LF, y=Pop))+geom_point()+xlim(c(-0.5,0.5))
q2<-ggplot(coef, aes(x=LF, y=PC))+geom_point()+xlim(c(-0.5,0.5))
plot_grid(q0,q1,q2, nrow=1)


makePlot<-function(pos){
	testing = data.frame(
			GT = unlist(GT[which(GT$POS == pos),-c(1:3)]),
			IAF = unlist(subset[which(subset$POS == pos),-c(1:3)]), 
			Qual = samples$average_quality_of_mapped_bases,
			Pop = samples$Pop)
			
	ggplot(testing, aes(x=Qual, y=GT, color=IAF))+
				geom_point()+facet_wrap(.~Pop, ncol=3)+
				guides(color=F)+geom_smooth(method = "glm", 
     			method.args = list(family = "binomial"), 
     			se = FALSE)
}

coef = cbind(coef, GT[c(1:10000),c(1:3)])
coef1 = merge(coef, SigPos, by=c('CHROM','POS','ID'))

devLF = lapply(fits, function(x) anova(x, test='Chi')[[2,2]])
devPop = lapply(Pop, function(x) anova(x, test='Chi')[[3,2]])
devPC = lapply(PC, function(x) anova(x, test='Chi')[[7,2]])

dev = data.frame('LFd' = unlist(devLF),
				  'Popd'=unlist(devPop),
				  'PCd' = unlist(devPC))
				  
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