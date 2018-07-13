library(reshape2)
args = commandArgs(trailingOnly=T)

#Pop and Ind Names
NamePop<-read.table('~/genomes/genomes/PopNames/Name.BigPop.Pop.txt', col.names=c('Name', 'BigPop', 'Pop'))
NamePop<-NamePop[which(NamePop$Pop != 'NAG'),]

#Quality per individual
Meta<-read.table("~/Dropbox/LukeTemp/1000GenomesMetaData.txt")
Meta$SUBMISSION.DATE<-as.Date(Meta$SUBMISSION.DATE, "%Y-%m-%d")  
SubMeta<-unique(Meta[c('Name','Pop', 'average_quality_of_mapped_bases')])
print('Welcome to the R script')
#Read Genotype Data
tab5rows <- read.table(args[1], header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
GT<-read.table(args[1],header=F, colClasses= classes, nrow=100000)
names(GT)=c('Chr','Pos', as.character(NamePop$Name))
print('Table is Loaded')
#Melt, reshape
melt.GT<-melt(GT,id=c('Chr','Pos'))
print('Table is Melted')
#Merge with Quality
GT.Qual<-merge(SubMeta,melt.GT, by.y='variable',by.x='Name')
#make it a factor
GT.Qual$Pop<-as.factor(GT.Qual$Pop)
GT.Qual$Pos<-as.factor(GT.Qual$Pos)
####
GT.Qual$value<-gsub(2,1, GT.Qual$value) #collapse GT
####
GT.Qual $value<-as.numeric(as.character(GT.Qual $value))
#Get number of alleles per position per population
dfSums<-GT.Qual %>% group_by(Pos,Pop) %>% summarize(SumOver=length(table(value))) %>% as.data.frame
GT.Qual<-merge(GT.Qual,dfSums,by=c('Pos','Pop'))
GT.Qual<-GT.Qual[which(GT.Qual$SumOver==2),]

GT.Qual$p=''
#For each pop
for(pop in unique(GT.Qual$Pop)){
	print(pop)
	#subset
gtQual<-droplevels(GT.Qual[which(GT.Qual$Pop==pop),])
gtQual$value<-as.numeric(as.character(gtQual$value))
#for each position, run glm
stat<-with(gtQual, 
		by(gtQual, gtQual$Pos, 
			function(x) glm(value ~ Quality, data=x, family=binomial)))
#extract P value					
GT.Qual[which(GT.Qual$Pop==pop),]$p<-sapply(stat,function(x) summary(x)$coef[2,4])				
}
GT.Qual$Plog10<--log10(as.numeric(GT.Qual$p))
write.table(GT.Qual,file=args[2], quote=F, row.names=F)