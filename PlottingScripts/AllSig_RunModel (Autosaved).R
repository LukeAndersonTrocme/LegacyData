library(dplyr)
library(data.table)
library(ggplot2)
library(glm2)

#read colors per pop (used in 1kG Nature)
colors<-read.table("~/Dropbox/LukeTemp/NamePopBigPop_color.html.txt", comment.char='>',col.names=c('Pop','hex','BigPop'))
MyColour <- as.character(colors $hex)
names(MyColour) <- colors $Pop

Sig = fread('/Users/luke/Documents/QualityPaper/Misc/SignificantVariants.csv',col.names=c("CHROM", "POS", "AF", "Rawlog10P", "Repeat", "Indel", "rsID", "log10P_0.01", "P_0.01", "p"))

#read Genotypes
GT = fread('/Users/luke/Documents/QualityPaper/sig/SigVar_1kGP_clean_GT.txt')
ColNames = read.table('/Users/luke/genomes/genomes/hg19/Genotypes/ColNames.txt', header=F)
names(GT) = as.character(unlist(ColNames))
samples = fread('/Users/luke/Documents/PCAperPop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt')
samples$Pop = factor(samples$Pop, levels=c('FIN','GBR','CEU','IBS','TSI','CHS','CDX','CHB','JPT','KHV','GIH','STU','PJL','ITU','BEB','PEL','MXL','CLM','PUR','ASW','ACB','GWD','YRI','LWK','ESN','MSL'))

model = apply(GT[,-(1:3)], 1, function(x)
						glm2(x ~
							samples$PC1 + 
							samples$PC2 +
							samples$PC3 +
							samples$PC4 +
							samples$PC5 + 
							samples$average_quality_of_mapped_bases,))
#get coef of Qual per site
coef = lapply(model, function(x) coef(x)[[7]])				

dev = lapply(model, function(x) anova(x, test='Chi')[[7,2]])

coef = cbind(unlist(coef), GT[,c(1:3)])
coef = cbind(coef, unlist(dev))
coef$ID = as.numeric(as.character(coef$ID))
coef$pval = -log10(pchisq(coef$V2, 1, lower.tail=F))
names(coef) = c('Beta','CHROM','POS','AF','deviance','-log10(p)')
write.table(coef,'/Users/luke/Documents/QualityPaper/sig/SigVar_1kGP_Coef.txt', quote=F, row.names=F)

ggplot(coef, aes(x=ID, y=V1, color=pval))+geom_point()+theme_bw()+labs(x='Allele Frequency',y='Beta')

all = merge(coef, Sig, by=c('CHROM','POS','AF'))

ggplot(all, aes(x=AF, y= Beta, color= log10P_0.01))+geom_point()+theme_bw()+labs(x='Allele Frequency',y='Beta')
