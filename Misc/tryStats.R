library(data.table)
args = commandArgs(trailingOnly=T)
samples = fread('/Users/luke/Documents/PCAperPop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt')

Input <- fread(paste("gunzip -c ",args[1]), header=T)
#Input <- fread('gzcat ~/genomes/genomes/hg19/Genotypes/Chr22_GT.csv.gz', header=F, nrows = 10000)

chunks = split(Input, cumsum((1:nrow(Input)-1)%%1000==0))
i=1
Pop.Qual = list()
PC.Qual = list()
PopPC.Qual = list()
for(chunk in chunks){
	Pop.Qual[[i]]= apply(chunk, 1, function(x) 
					anova(
						glm(x ~ 
							samples$Pop +
							samples$average_quality_of_mapped_bases,
						family=binomial),
					test="Chi")[3,2])

	PC.Qual[[i]] = apply(chunk, 1, function(x) 
					anova(
						glm(x ~  
							samples$PC1 + 
							samples$PC2 + 
							samples$average_quality_of_mapped_bases,
						family=binomial),
					test="Chi")[4,2])					
	
	PopPC.Qual[[i]] = apply(chunk, 1, function(x) 
					anova(
						glm(x ~ 
							samples$Pop + 
							samples$PC1 + 
							samples$PC2 + 
							samples$average_quality_of_mapped_bases,
						family=binomial),
					test="Chi")[5,2])
			
	i=i+1
	}
	
Pop.Qual1 = unlist(Pop.Qual)
PC.Qual1 = unlist(PC.Qual)
PopPC.Qual1 = unlist(PopPC.Qual)

plt <-data.frame('popQual' = -log10(pchisq(Pop.Qual1, 1, lower.tail=F)),
				'PCQual' = -log10(pchisq(PC.Qual1, 1, lower.tail=F)), 
				'PopPCQual' = -log10(pchisq(PopPC.Qual1, 1, lower.tail=F)))
library(ggplot2)
library(cowplot)				
p1=ggplot(plt, aes(x=plt$popQual, y=plt$PopPCQual))+geom_point()
p2=ggplot(plt, aes(x=plt$PCQual, y=plt$PopPCQual))+geom_point()	

plot_grid(p1,p2)			