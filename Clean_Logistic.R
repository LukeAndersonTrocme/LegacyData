library(reshape2)
library(magrittr)
args = commandArgs(trailingOnly=T)

#Pop and Ind Names #~/genomes/genomes/PopNames/Name.BigPop.Pop.txt
NamePop<-read.table('~/Documents/PhD/data/Name.BigPop.Pop.txt', col.names=c('Name', 'BigPop', 'Pop'))
NamePop<-NamePop[which(NamePop$Pop != 'NAG'),]

#Quality per individual
SubMeta <-read.table('~/Dropbox/LukeTemp/SubMeta.txt')
print('Welcome to the R script')
#Read Genotype Data
tab5rows <- read.table(args[1], header = F, nrows = 20)
classes <- sapply(tab5rows, class)
GT<-read.table(args[1],header=F, colClasses= classes, nrow=100)
names(GT)=c('Chr','Pos', as.character(NamePop$Name))
print('Table is Loaded')
#Melt, reshape
melt.GT<-melt(GT,id=c('Chr','Pos'))
rm(GT)
print('Table is Melted')
#Merge with Quality
GT.Qual<-merge(SubMeta,melt.GT, by.y='variable',by.x='Name')
#rm(melt.GT)
#make it a factor
GT.Qual$Pop<-as.factor(GT.Qual$Pop)
GT.Qual$Pos<-as.factor(GT.Qual$Pos)
####
GT.Qual$value<-gsub(2,1, GT.Qual$value) #collapse GT
####
GT.Qual $value<-as.numeric(as.character(GT.Qual $value))
#Get number of alleles per position per population
dfSums <-as.data.frame(table(GT.Qual$value,GT.Qual$Pop,GT.Qual$Pos))
names(dfSums)<-c('value','Pop','Pos','Freq')
#dfSums<-GT.Qual %>% group_by(Pos,Pop) %>% summarize(SumOver=length(table(value,Pos,Pop))) %>% as.data.frame
GT.Qual<-merge(GT.Qual,dfSums,by=c('Pos','Pop'))
GT.Qual<-GT.Qual[which(GT.Qual$Freq>0),]
GT.Qual<-GT.Qual[,c('Pos','Pop','Name','average_quality_of_mapped_bases','Chr','value.x')]
test<-droplevels(GT.Qual[which((GT.Qual$Pop=='JPT')&(GT.Qual$Pos=='16554752')),])
GT.Qual$p=''
#For each pop
for(pop in unique(GT.Qual$Pop)){
	print(pop)
	#subset
gtQual<-droplevels(GT.Qual[which(GT.Qual$Pop==pop),])
gtQual$value<-as.numeric(as.character(gtQual$value.x))
#for each position, run glm
stat<-with(gtQual, 
		by(gtQual, gtQual$Pos, 
			function(x) glm(value.x ~ average_quality_of_mapped_bases, data=x, family=binomial)))
#extract P value					
GT.Qual[which(GT.Qual$Pop==pop),]$p<-sapply(stat,function(x) summary(x)$coef[2,4])				
}
GT.Qual$Plog10<--log10(as.numeric(GT.Qual$p))

write.table(GT.Qual,file=args[2], quote=F, row.names=F)