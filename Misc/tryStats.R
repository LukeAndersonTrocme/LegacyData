chnk <- fread('~/Desktop/chr22.csv', header=T)
chnk$V1=NULL
samples = fread('/Users/luke/Documents/PCAperPop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt')
ptm <- proc.time()
outt <- apply(chunk, 1, function(x) anova(glm(x ~ samples$Pop + samples$PC1 + samples$PC2 + samples$average_quality_of_mapped_bases, family=binomial),test="Chi")[5,2])
print(proc.time() - ptm)