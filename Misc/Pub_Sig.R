library(dplyr)
library(data.table)
library(ggplot2)
library(glm2)

#read colors per pop (used in 1kG Nature)
colors<-read.table("~/Dropbox/LukeTemp/NamePopBigPop_color.html.txt", comment.char='>',col.names=c('Pop','hex','BigPop'))
MyColour <- as.character(colors $hex)
names(MyColour) <- colors $Pop

path = '/Users/luke/genomes/genomes/hg19/Genotypes/'

GT = fread('~/Desktop/PublishedHits.txt')
ColNames = read.table(paste(path,"ColNames.txt", sep=''), header=F)
names(GT)= as.character(unlist(ColNames))
samples = fread('/Users/luke/Documents/PCAperPop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt')
samples$Pop = factor(samples$Pop, levels=c('FIN','GBR','CEU','IBS','TSI','CHS','CDX','CHB','JPT','KHV','GIH','STU','PJL','ITU','BEB','PEL','MXL','CLM','PUR','ASW','ACB','GWD','YRI','LWK','ESN','MSL'))

SigPub = fread('~/Documents/QualityPaper/sig/SignificantPublications.txt')


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
					predict(glm2(x ~
							samples$PC1 + 
							samples$PC2 +
							samples$PC3 +
							samples$PC4 +
							samples$PC5 + 
							samples$average_quality_of_mapped_bases,
						family=binomial)))

exp.coef = apply(GT[,-(1:3)], 1, function(x)
					exp(coefficients(glm2(x ~
							samples$PC1 + 
							samples$PC2 +
							samples$PC3 +
							samples$PC4 +
							samples$PC5 + 
							samples$average_quality_of_mapped_bases,
						family=binomial)))[[7]])

plt = cbind(samples, model)
plt = melt(plt, id=c('Name','Pop','average_quality_of_mapped_bases','PC1','PC2','PC3','PC4','PC5'))
plt = plt[order(plt$Pop),]
q = list()
for(p in seq(1,14)){

	sub = plt[which(plt$variable == paste('V',p,sep='')),]
	snp = GT[[p,2]]
	journal = SigPub[which(SigPub$POS == snp),]$JOURNAL
	pval = round(SigPub[which(SigPub$POS == snp),]$log10P_0.01,2)
	rsID = SigPub[which(SigPub$POS == snp),]$SNPS
	logOdds = round(exp.coef[[p]],2)
	
	title = paste(journal, rsID, '\n -log10(p) : ', pval, ' odds ratio : ', logOdds)
	
	q[[p]] = ggplot(sub,aes(x=value, y= average_quality_of_mapped_bases, color=Pop))+scale_color_manual(breaks= plt$Pop,values = MyColour)+guides(color=F)+geom_point()+ggtitle(title)+labs(x='Predicted',y='Quality')
}
l = ggplot(sub,aes(x=value, y= average_quality_of_mapped_bases, color=Pop))+scale_color_manual(breaks= plt$Pop,values = MyColour)+geom_point()+ggtitle(title)+labs(x='Predicted',y='Quality')+guides(colour = guide_legend(override.aes = list(shape = 15, size=5), ncol=6, title='Population'))

legend <- get_legend(l)

plot_grid(q[[5]],q[[7]],q[[9]],q[[2]],q[[13]],q[[6]],q[[10]],q[[8]],q[[12]],q[[11]],q[[3]],q[[14]],q[[4]],q[[1]],legend,ncol=3)

ggsave('~/Documents/QualityPaper/Figures/PublishedSNPs.jpg',height=20, width=18)					
		