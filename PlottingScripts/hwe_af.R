library(ggplot2)
library(colorspace)
library(cowplot)
library(data.table)

f1kGP='/Users/luke/genomes/genomes/1kGP_4bedFiltered/1kGP_GenomeWide_JPT_1.frq'
fnag='~/genomes/genomes/NAG_4bedFiltered/NAG_GenomeWide.4bed_filtered.freq.frq'
hweFile='/Users/luke/Documents/GWAS_Qual/Final_GWAS/1kGP_GenomeWide_JPT_sig6.hwe'

JPT_hwe <- fread('/Users/luke/genomes/genomes/1kGP_4bedFiltered/1kGP_GenomeWide.4bed_filtered_1.freq.hwe')
colnames(JPT_hwe)[c(1,2,8)]<-c('CHROM','POS','HWE_p')
#LogisticRegression
JPT_reg <- fread('~/Documents/QualityPaper/Misc/GenomeWide.Regression_JPT.csv')
colnames(JPT_reg)[c(1,2)]<-c('CHROM','POS')
JPT_reg$CHROM<-as.numeric(as.character(JPT_reg$CHROM))
JPT_reg$POS<-as.numeric(as.character(JPT_reg$POS))
JPT_reg$dev<-as.numeric(as.character(JPT_reg$dev))
JPT_reg<-JPT_reg[complete.cases(JPT_reg)]

#Allele Frequency files
JPT<-fread(f1kGP, fill=T, col.names=c('CHROM','POS','N_ALLELES','N_CHR','JPT_AF','JPT_MAF'), colClasses=c('numeric'))
NAG<-fread(fnag, fill=T,col.names=c('CHROM','POS','N_ALLELES','N_CHR','NAG_AF','NAG_MAF'), colClasses=c('numeric'))

#merge them together to get joint frequency spectrum
join<-merge(JPT[which(JPT$JPT_AF <= 1),],
 			NAG[which(NAG$NAG_AF <= 1),], 
 			by=c('CHROM','POS'), all=T)
join[is.na(join)]<-0

join=join[,c('CHROM','POS','JPT_MAF','JPT_AF','NAG_MAF','NAG_AF')]

join$CHROM<-as.numeric(as.character(join$CHROM))
join$POS<-as.numeric(as.character(join$POS))
join$JPT_AF<-as.numeric(as.character(join$JPT_AF))
join$JPT_MAF<-as.numeric(as.character(join$JPT_MAF))
join$NAG_AF<-as.numeric(as.character(join$NAG_AF))
join$NAG_MAF<-as.numeric(as.character(join$NAG_MAF))
join<-join[complete.cases(join),]

join <- merge(join, JPT_hwe, by=c('CHROM','POS'))
join <- merge(join, JPT_reg, by=c('CHROM','POS'))
join$devP<-pchisq(join$dev, 1, lower.tail=T)
join$log10_devP<- - pchisq(join$dev, 1, lower.tail=F, log.p=T)/log(10)

miss<-join[which(join$JPT_MAF > 0.01 & join$NAG_MAF == 0),]
nomiss<-join[ ! which(join$JPT_MAF > 0.01 & join$NAG_MAF == 0),]

dev1<-as.data.frame(quantile(
				miss$log10_devP,
				probs = seq(0, 1, 0.0001)))
				
dev2<-as.data.frame(quantile(
				nomiss$log10_devP,
				probs = seq(0, 1, 0.0001)))
				
dev3<-cbind(dev1,dev2)

dev=ggplot(dev3, 
		aes(y=dev3[[1]], 
		x=dev3[[2]]))+	
		geom_point()+
		geom_abline(slope=1,intercept=0,color='blue')+
		labs(x='not missing from Nagahama',
		y='missing from Nagahama')+
		ggtitle('Logistic Regression GWAS')

hwe1<-as.data.frame(quantile(miss$HWE_p,probs = seq(0, 1, 0.0001)))
hwe2<-as.data.frame(quantile(nomiss$HWE_p,probs = seq(0, 1, 0.0001)))
hwe3<-cbind(hwe1,hwe2)

hwe=ggplot(hwe3, aes(y=hwe3[[1]], x=hwe3[[2]]))+geom_point()+geom_abline(slope=1,intercept=0,color='blue')+labs(x='not missing from Nagahama',y='missing from Nagahama')+ggtitle('Hardy-Weinberg')

plot_grid(dev,hwe)

hist<-ggplot(miss, 
		aes(x=-log10(HWE_p)))+
		geom_histogram(binwidth=0.01)+
		xlim(c(0,10))+
		labs(x='Hardy Weinberg for SNPs\nmissing from Nagahama data')
	
hist1<-ggplot(nomiss, 
		aes(x=-log10(HWE_p)))+
		geom_histogram(binwidth=0.01)+
		xlim(c(0,10))+
		labs(x='Hardy Weinberg for SNPs\npresent in Nagahama data')
	
maf<-ggplot(miss, 
		aes(x=JPT_MAF,
			y=-log10(HWE_p)), 
			alpha=0.1)+
			geom_point()
			
maf1<-ggplot(nomiss[which(nomiss$CHROM>18),], 
			aes(x=JPT_MAF,
			y=-log10(HWE_p)), 
			alpha=0.1)+
			geom_point()

histr<-ggplot(miss, 
			aes(x=devP))+
			geom_histogram(binwidth=0.01)+
			xlim(c(0,10))+
			labs(x='GWAS for SNPs\nmissing from Nagahama data')
			
histr1<-ggplot(nomiss, 
			aes(x= devP))+
			geom_histogram(binwidth=0.01)+
			xlim(c(0,10))+
			labs(x='GWAS for SNPs\npresent in Nagahama data')
			
mafr<-ggplot(miss, 
			aes(x=JPT_MAF,
			y= devP), 
			alpha=0.1)+
			geom_point()+
			ylim(c(0,20))
			
mafr1<-ggplot(nomiss[which(nomiss$CHROM>18),], 
			aes(x=JPT_MAF,
			y= devP), 
			alpha=0.1)+
			geom_point()+
			ylim(c(0,20))

plot_grid(hist, hist1, histr, histr1,maf, maf1, mafr, mafr1, ncol=2)