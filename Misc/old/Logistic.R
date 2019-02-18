library(reshape2)
library(ggplot2)
library(cowplot)
args = commandArgs(trailingOnly=T)

#Pop and Ind Names
NamePop<-read.table('~/genomes/genomes/PopNames/Name.BigPop.Pop.txt', col.names=c('Name', 'BigPop', 'Pop'))
NamePop<-NamePop[which(NamePop$Pop != 'NAG'),]

#Quality per individual
Meta<-read.table("~/Dropbox/LukeTemp/1000GenomesMetaData.txt")
Meta$SUBMISSION.DATE<-as.Date(Meta$SUBMISSION.DATE, "%Y-%m-%d")  
SubMeta<-unique(Meta[c('Name','Pop', 'average_quality_of_mapped_bases')])
SubMeta$BinQual<-0
SubMeta[which(SubMeta$average_quality_of_mapped_bases<30),]$BinQual<-1
##Transformed Quality Score
NorMeta<-read.table('~/genomes/genomes/1000GenomesQual_INT.pheno',col.names=c('Name','ID','Quality'))
SubMeta<-merge(SubMeta,NorMeta[,c('Name','Quality')], by='Name')

#Read Genotype Data
GT<-read.table(args[1],header=F)
names(GT)=c('Chr','Pos', as.character(NamePop$Name))

#Melt, reshape
melt.GT<-melt(GT,id=c('Chr','Pos'))
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
GT.Qual$intercept=''
GT.Qual$beta=''
GT.Qual$p=''


for(pop in unique(GT.Qual$Pop)){
gtQual<-droplevels(GT.Qual[which(GT.Qual$Pop==pop),])
gtQual$value<-as.numeric(as.character(gtQual$value))
print(pop)
stat<-with(gtQual, 
		by(gtQual, gtQual$Pos, 
			function(x) glm(value ~ Quality, data=x, family=binomial)))
			
GT.Qual[which(GT.Qual$Pop==pop),]$intercept<-sapply(stat,function(x) coef(x)[1])
GT.Qual[which(GT.Qual$Pop==pop),]$beta<-sapply(stat,function(x) coef(x)[2])			
GT.Qual[which(GT.Qual$Pop==pop),]$p<-sapply(stat,function(x) summary(x)$coef[2,4])				
}
GT.Qual$Pos<-as.numeric(as.character(GT.Qual$Pos))
GT.Qual$beta<-as.numeric(as.character(GT.Qual$beta))
GT.Qual$intercept<-as.numeric(as.character(GT.Qual$intercept))
GT.Qual$Plog10<--log10(as.numeric(GT.Qual$p))

ggplot(GT.Qual, aes(x=intercept, fill=Pop))+geom_density()+facet_wrap(.~Pop, ncol=2, scales='free')+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank())+guides(fill=FALSE)+geom_vline(xintercept=0,color='grey70')+ggtitle('Density of Intercepts for \nLogistic Regression of Genotype by Quality')
ggsave('~/Documents/Regression/JPT_sig6_LogisticPerPop_interceptDensity.jpg')

ggplot(GT.Qual, aes(x=intercept,y=beta, color=Pop))+geom_point()+facet_wrap(.~Pop, ncol=2, scales='free')+guides(color=FALSE)+geom_vline(xintercept=0,color='grey70')+ggtitle('Intercept by Beta')
ggsave('~/Documents/Regression/JPT_sig6_LogisticPerPop_interceptBeta.jpg')

ggplot(GT.Qual, aes(x=Plog10, fill=Pop))+geom_density()+facet_wrap(.~Pop, ncol=2)+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank())+guides(fill=FALSE)+geom_vline(xintercept=0,color='grey70')+ggtitle('Density of Intercepts for \nLogistic Regression of Genotype by Quality')
ggsave('~/Documents/Regression/JPT_sig6_Plog10PerPop_interceptDensity.jpg')

ggplot(GT.Qual, aes(x=Pop,y=Plog10, fill=Pop))+geom_violin()+guides(fill=FALSE)
ggsave('~/Documents/Regression/JPT_sig6_Plog10PerPop.jpg',height=8,width=11)


ggplot(GT.Qual[which((GT.Qual$Pop=='JPT')&(GT.Qual$Pos== 100033562)),], aes(x= average_quality_of_mapped_bases, y=value))+geom_point()+scale_y_continuous(breaks=c(0,1))+labs(y='Genotype',x='Quality')+geom_smooth(method = "glm",     method.args = list(family = "binomial"), se = FALSE,size=0.5) 




##Trash

if(FALSE){
	#Statistical Test
statTest <- function(arg1){
	sub.glm=glm(formula=value ~ average_quality_of_mapped_bases, data=arg1, family=binomial)
	return(sub.glm)
}

stat<-with(GT.Qual, 
		by(GT.Qual, GT.Qual$Pos, 
			function(x) glm(value ~ average_quality_of_mapped_bases + Pop, data=x, family=binomial)))
			
GT.Qual$p<-sapply(stat,function(x) coef(x)[2])

#### HAMMING DISTANCE
Hamming<-GT.Qual %>% group_by(Pos,Pop) %>% summarize(Distance=sum(value != BinQual)) %>% as.data.frame
Hamming<-merge(Hamming, unique(GT.Qual[,c('Chr','Pos','Pop')]), by=c('Pos','Pop'))
ggplot(Hamming,aes(x=reorder(Pos,Chr),y=Distance,color=Pop))+geom_point()+facet_wrap(.~Pop)+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

ggplot(Hamming,aes(x=Distance,fill=Pop))+geom_density()+facet_wrap(.~Pop,ncol=2)


library(GGally)
Hsub<-droplevels(Hamming[which(Hamming$Pop %in% c('JPT','YRI','PUR','GBR','CEU')),])
ggpairs(Hsub[,c('Distance','Pop')], aes(colour = Pop, alpha = 0.4))

ggplot(Hamming,aes(x=Distance, fill=Pop))+geom_density()+facet_wrap(.~Pop,ncol=2)

}