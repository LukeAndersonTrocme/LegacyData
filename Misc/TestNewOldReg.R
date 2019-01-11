library(data.table)

f1<-'~/Documents/GenomeWideRegression/GenomeWide.MeanDev_ALLPOPS.csv'
f2<-'~/Documents/QualityPaper/Misc/Reg_ALL_SITES_5col.txt'

new<-fread(f1)
new $Chr<-as.numeric(as.character(new $Chr))
new $Pos<-as.numeric(as.character(new $Pos))
new $SumDev <-as.numeric(as.character(new $SumDev))
new $Count<-as.numeric(as.character(new $Count))
new$MeanDev<-NULL

old<-fread(f2)
old $Chr<-as.numeric(as.character(old $Chr))
old $Pos<-as.numeric(as.character(old $Pos))
old $SumDev<-as.numeric(as.character(old $SumDev))
old $Count<-as.numeric(as.character(old $Count))
old$MeanDev<-NULL

test<-unique(merge(old[which(old$Chr==1),], new[which(new$Chr==1),], by=c('Chr','Pos')))

test$p.new<-pchisq(test$SumDev.y, test$Count.y, lower.tail=F)

test$p.old<-pchisq(test$SumDev.x, test$Count.x, lower.tail=F)

test$diff<-test$Count.y-test$Count.x

library(ggplot2)

ggplot(test, aes(x= -log10(p.new), y= -log10(p.old)))+geom_point()

chr22<-list.files(path="/Volumes/gravel/luke_projects/1000Genomes/FinalRegression/",pattern="CHR22.Regression",full.names=T)


weird<-data.frame()
for(f in chr22){
	print(f)
newReg<-fread(f)
newReg $Chr<-as.numeric(as.character(newReg $Chr))
newReg $Pos<-as.numeric(as.character(newReg $Pos))
newReg $dev<-as.numeric(as.character(newReg $dev))
name<-newReg$Pop[[1]]

oldReg<-fread(paste('/Volumes/gravel/luke_projects/1000Genomes/Regression/CHR22.Regression_',name,'.csv', sep=''))
oldReg $Chr<-as.numeric(as.character(oldReg $Chr))
oldReg $Pos<-as.numeric(as.character(oldReg $Pos))
oldReg $dev<-as.numeric(as.character(oldReg $dev))

test2<-unique(merge(oldReg, newReg, by=c('Chr','Pos')))
test2$p.new<-pchisq(test2$dev.y, 1, lower.tail=F)
test2$p.old<-pchisq(test2$dev.x, 1, lower.tail=F)

test2$diff <- -log10(test2$p.new) - -log10(test2$p.old)

Diff <- test2[which(test2$diff > 5),]
if(nrow(Diff) > 0){
format <- fread(paste('gzcat ~/Documents/Regression/CHR22.Format_',name,'.csv.gz',sep=''))

format <- format[which(format$Pos %in% Diff$Pos), c('Chr','Pos','value','average_quality_of_mapped_bases')]

format <- aggregate(. ~ Pos + value, data=format, mean)

Diff <- merge(Diff, format, by='Pos')

weird <- rbind(weird, Diff)
}
}
ggplot(test2, 
	aes(x= -log10(p.new), 
		y= -log10(p.old)))+
	geom_abline(slope=1, intercept=0, color='blue')+
	geom_point(shape=1)+
	ggtitle(oldReg$Pop[[1]])
	
ggsave(paste('~/Documents/QualityPaper/Misc/CheckReg/Chr22',name,'.jpg',sep=''))
}

format<-