library(data.table)
library(ggplot2)
f1<-'/Users/luke/Documents/PythonRegression/Regression_Chr_19.csv'
f2<-'/Users/luke/Documents/QualityPaper/Misc/Reg_chr19.txt'

new<-fread(f1)
new $CHROM <-as.numeric(as.character(new $CHROM))
new $POS <-as.numeric(as.character(new $POS))
new $AF<-as.numeric(as.character(new $AF))
new $pop_PCs <-as.numeric(as.character(new $pop_PCs))
new $pop_PCs_Q <-as.numeric(as.character(new $pop_PCs_Q))
new$dev <- new$pop_PCs - new $pop_PCs_Q


old<-fread(f2,col.names=c('CHROM', 'POS', 'SumDev', 'MeanDev', 'Count'))
old $CHROM <-as.numeric(as.character(old $CHROM))
old $POS<-as.numeric(as.character(old $POS))
old $SumDev<-as.numeric(as.character(old $SumDev))
old $Count<-as.numeric(as.character(old $Count))
old$MeanDev<-NULL

newR <- fread('/Volumes/gravel/luke_projects/1000Genomes/NewR_Regression/CHR19.Regression.csv', col.names='NEWdev')
newPos <- fread('~/genomes/genomes/hg19/Genotypes/Chr19_Pos.csv')
newR<- cbind(newR, head(newPos, nrow(newR)))

newR $CHROM <-as.numeric(as.character(newR $CHROM))
newR $POS <-as.numeric(as.character(newR $POS))
newR $AF <- as.numeric(as.character(newR $AF))
newR$NEWdev <- as.numeric(as.character(newR $NEWdev))

test<-unique(merge(old[which(old$CHROM==19),], new[which(new$CHROM==19),], by=c('CHROM','POS')))

test <-unique(merge(test, newR, by=c('CHROM','POS')))

test$pythonRegression<--log10(pchisq(test$dev, 1, lower.tail=F))
test$R.Regression<--log10(pchisq(test$NEWdev, 1, lower.tail=F))
test$OldMethod<--log10(pchisq(test$SumDev, test$Count, lower.tail=F))


plt1=ggplot(test[which(test$AF.x < 0.02)], aes(x= pythonRegression, y= R.Regression, color=AF.x))+geom_point()

plt2=ggplot(test[which(test$AF.x < 0.02)], aes(x= pythonRegression, y= OldMethod, color=AF.x))+geom_point()
plt3=ggplot(test[which(test$AF.x < 0.02)], aes(x= R.Regression, y= OldMethod, color=AF.x))+geom_point()
plot_grid(plt1,plt2,plt3, nrow=1)

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