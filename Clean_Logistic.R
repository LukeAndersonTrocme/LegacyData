library(reshape2)
library(magrittr)
library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=T)

#Pop and Ind Names
NamePop<-read.table('~/genomes/genomes/PopNames/Name.BigPop.Pop.txt', col.names=c('Name', 'BigPop', 'Pop'))
NamePop<-NamePop[which(NamePop$Pop != 'NAG'),]

#Quality per individual
SubMeta <-read.table('~/Dropbox/LukeTemp/SubMeta.txt')
print('Welcome to the R script')
#Read Genotype Data
tab5rows <- read.table(args[1], header = F, nrows = 20)
classes <- sapply(tab5rows, class)
GT<-fread(args[1],header=F, colClasses= classes)
names(GT)=c('Chr','Pos', as.character(NamePop$Name))

print('Table is Loaded')
print(object.size(GT), units="Mb")

out<-data.frame()
for(pop in unique(SubMeta$Pop)){
print(pop)
#subset
SubGT<-GT[,c('Chr','Pos',as.character(SubMeta[which(SubMeta$Pop==pop),]$Name)), with=F]
print(object.size(SubGT), units="Mb")
#remove pop from GT table to save memory
#GT<-GT[,-which(names(GT) %in% SubMeta[which(SubMeta$Pop==pop),]$Name)]
print('subset')
SubMelt<-melt(SubGT,id=c('Chr','Pos'))
print('melt')
rm(SubGT)
GT.Qual<-droplevels(merge(SubMeta, SubMelt, by.y='variable',by.x='Name'))
print(object.size(GT.Qual), units="Mb")
print('merge')
rm(SubMelt)
GT.Qual$value<-as.numeric(as.character(gsub(2,1, GT.Qual$value))) #collapse GT
dfSums<-GT.Qual %>% group_by(Pos) %>% summarize(SumOver=length(table(value))) %>% as.data.frame
GT.Qual<-merge(GT.Qual,dfSums,by=c('Pos'))
GT.Qual<-GT.Qual[which(GT.Qual$SumOver==2),]
#for each position, run glm
print('stats')
stat<-with(GT.Qual, 
		by(GT.Qual, GT.Qual $Pos, 
			function(x) glm(value ~ average_quality_of_mapped_bases, data=x, family=binomial)))
#extract P value
GT.Qual$p=''					
GT.Qual$p<--log10(as.numeric(sapply(stat,function(x) summary(x)$coef[2,4])))
print('writing')
write.table(GT.Qual,file=paste(args[2],pop,'.txt',sep=''), quote=F, row.names=F)
rm(GT.Qual)				
}