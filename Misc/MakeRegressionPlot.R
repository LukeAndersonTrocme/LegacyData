library(dplyr)
library(data.table)
library(ggplot2)

args = commandArgs(trailingOnly=T)
Chr = args[1]
Pos = args[2]

path = '/Users/luke/genomes/genomes/hg19/Genotypes/'
fName = paste("/Users/luke/Desktop/R",Pos,".temp.txt", sep='')
command = paste("gzcat ",path,"CHR",Chr,".Genotypes.txt.gz | awk '$2 == ",Pos,"' > ",fName,sep='')
print(command)
system(command)
GT <- fread(fName)
file.remove(fName)

ColNames = read.table(paste(path,"ColNames.txt", sep=''), header=F)
names(GT)= as.character(unlist(ColNames))
samples = fread('/Users/luke/Documents/PCAperPop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt')

mGT = melt(GT[,-(1:3)])
mGT <- merge(samples, mGT, by.y='variable',by.x='Name')
cGT <- mGT %>% group_by(Pop) %>% summarize(var = sum(value))
keep <- cGT[which(cGT$var > 0),]$Pop

ggplot(mGT[which(mGT$Pop %in% keep),], aes(x= average_quality_of_mapped_bases, y=value)) +
 	facet_wrap(.~Pop, ncol=1) + 
 	geom_smooth(method = "glm", 
     			method.args = list(family = "binomial"), 
     			se = FALSE) +
     geom_point() +
     theme_classic() + 
     labs(y='genotypes', title = paste(Chr," : ",Pos)) +
     theme(plot.title = element_text(hjust = 0.5))

ggsave(paste('~/Documents/QualityPaper/Misc/Rplot_LogisticReg/',Chr,'_',Pos,'.jpg',sep=''), height = 8, width=6)
