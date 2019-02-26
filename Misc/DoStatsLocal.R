library(data.table)
library(dplyr)
library(glm2)

args = commandArgs(trailingOnly=T)
#file with sample info
samples = fread('/Users/luke/Documents/PCAperPop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt')
#column names for genotype file
ColNames = read.table('/Users/luke/genomes/genomes/hg19/Genotypes/ColNames.txt', header=F)

fileName = paste("gzcat ",args[1], sep='')
print(fileName)
genotypes <- fread(fileName, sep='\t', header=F)
names(genotypes)= as.character(unlist(ColNames))

print('Read file')
chunks = split(genotypes, cumsum((1:nrow(genotypes)-1)%%10000==0))
i=1
for(chunk in chunks){
	print(i)
	#split data
	positions = chunk[,c(1:3)]
	gt = chunk[,-(1:3)]
	out = apply(gt, 1, function(x) 
				anova(
					glm2(x ~
							samples$Pop +
							samples$average_quality_of_mapped_bases,
						family=binomial),
					test="Chi")[3,2])
					
	positions$dev.p = -log10(pchisq(out, 1, lower.tail=F))

	#write results to file
	if(i == 1){
		write.table(positions,file=args[2], 
		quote=F, row.names=F, col.names=T, append=FALSE)}
	if(i > 1){
    	write.table(positions,file=args[2], 
    	quote=F, row.names=F, col.names=F, append=TRUE)}
    
    i=i+1
    	}
    	