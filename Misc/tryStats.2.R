library(data.table)
library(dplyr)

args = commandArgs(trailingOnly=T)
samples = fread('/Users/luke/Documents/PCAperPop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt')
samples$Q = samples$average_quality_of_mapped_bases > 30
#get positions of individuals from each population
pops <- split(rownames(samples), samples$Pop)

Input <- fread(paste("gunzip -c ",args[1]), header=T)
#Input <- fread('gzcat /Users/luke/genomes/genomes/hg19/Genotypes/Chr22_GT_Pos.csv.gz ', header=F, nrows = 10000)
newInput = fread('gzcat ~/genomes/genomes/hg19/Genotypes/CHR22.Genotypes.txt.gz',sep='\t',nrows = 1000)
ColNames = read.table('/Users/luke/genomes/genomes/hg19/Genotypes/ColNames.txt', header=F)
names(newInput)= as.character(unlist(ColNames))

#Get chromsomes and positions
ChrPos = newInput[,c(1,2,3)]
names(ChrPos) = c('CHROM','POS','AF')
#newInput[,c(1:3)]=NULL

chunks = split(newInput, cumsum((1:nrow(newInput)-1)%%1000==0))

i=1
meta = list()
multi = list()
for(chunk in chunks){
	print(i)

	#Do the new test, takes pops as covariate
	PopPC.Qual = apply(chunk[,-(1:3)], 1, function(x) 
					anova(
						glm(x ~ 
							samples$Pop + 
							samples$PC1 + 
							samples$PC2 + 
							samples$average_quality_of_mapped_bases,
						family=binomial),
					test="Chi")[5,2])
					
	PopPC5.Qual = apply(chunk[,-(1:3)], 1, function(x) 
					anova(
						glm(x ~ 
							samples$Pop + 
							samples$PC1 + 
							samples$PC2 +
							samples$PC3 + 
							samples$PC4 +
							samples$PC5 +
							samples$average_quality_of_mapped_bases,
						family=binomial),
					test="Chi")[8,2])					
	PC5.Qual = apply(chunk[,-(1:3)], 1, function(x) 
					anova(
						glm(x ~ 
							samples$PC1 + 
							samples$PC2 +
							samples$PC3 + 
							samples$PC4 +
							samples$PC5 +
							samples$average_quality_of_mapped_bases,
						family=binomial),
					test="Chi")[7,2])


	PC.Qual = apply(chunk[,-(1:3)], 1, function(x) 
					anova(
						glm(x ~  
							samples$PC1 + 
							samples$PC2 + 
							samples$average_quality_of_mapped_bases,
						family=binomial),
					test="Chi")[4,2])
	
	logF = apply(chunk[,-(1:3)], 1, function(x) 
					logistf(x ~ 
							samples$PC1 + 
					   		samples$PC2 + 
					   		samples$average_quality_of_mapped_bases)$prob[[4]])				   				
	
	combined <-data.frame('CHROM' = chunk[,1], 'POS' = chunk[,2],
				'PopPC5.Qual' = -log10(pchisq(PopPC5.Qual, 1, lower.tail=F)),
				'PCQual' = -log10(pchisq(PC.Qual, 1, lower.tail=F)), 
				'PopPCQual' = -log10(pchisq(PopPC.Qual, 1, lower.tail=F)),
				'logF'=logF)
	
	multi[[i]] = combined
	
	i=i+1				
}

meta1 = bind_rows(meta)
multi1 = bind_rows(multi)	

plt = merge(combined, metaOut, by=c('CHROM','POS'))
	
p1 = ggplot(plt, aes(x=plt$PopPC5.Qual, y=plt$PopPCQual))+geom_point()+geom_abline(slope=1, intercept=0, color='blue')
p2 = ggplot(plt, aes(x=plt$PCQual, y=plt$PopPCQual))+geom_point()+geom_abline(slope=1, intercept=0, color='blue')	
p3 = ggplot(plt, aes(x=plt$logF.p, y=plt$PopPCQual))+geom_point()+geom_abline(slope=1, intercept=0, color='blue')
p4 = ggplot(plt, aes(y=plt$metaP, x=plt$PopPCQual))+geom_point()+geom_abline(slope=1, intercept=0, color='blue')
p5=ggplot(plt, aes(y=plt$metaP, x=plt$PCQual))+geom_point()+geom_abline(slope=1, intercept=0, color='blue')
p6=ggplot(plt, aes(y=plt$PC5.Qual, x=plt$PopPCQual))+geom_point()+geom_abline(slope=1, intercept=0, color='blue')
p7=ggplot(plt, aes(y=plt$PC5.Qual, x=plt$PCQual))+geom_point()+geom_abline(slope=1, intercept=0, color='blue')
plot_grid(p1,p4,p2,p3, p5)	
	
plot_grid(p5, p4, p6,p7)

	
	
##concatenate the data
m1=as.data.frame(meta[[1]],col.names=unique(samples$Pop))
m2=as.data.frame(meta[[2]],col.names=unique(samples$Pop))
m3=as.data.frame(meta[[3]],col.names=unique(samples$Pop))
m4=as.data.frame(meta[[4]],col.names=unique(samples$Pop))
m<-rbind(m1,m2)
m<-rbind(m, m3)
m<-rbind(m, m4)
m.Pos = bind_rows(positions)
#add chrom and position for merging
m=cbind(m.Pos,m)
m$NewTest = -log10(pchisq(unlist(PopPC.Qual), 1, lower.tail=F))
m=melt(m, id=c('CHROM','POS','ID','NewTest'))
m$popAF <- unlist(counts)
names(m)[5:6] = c('Pop','deviance')

chr22 = fread('/Users/luke/Documents/QualityPaper/Misc/Chr22_MetaReg_10000.csv')
chr22 $Chr<-as.numeric(as.character(chr22 $Chr))
chr22 $Pos<-as.numeric(as.character(chr22 $Pos))
chr22 $dev<-as.numeric(as.character(chr22 $dev))

preCombo = merge(chr22, m, by.x=c('Chr','Pos', 'Pop'),by.y=c('CHROM','POS','Pop'))
preComboAgg <- preCombo %>% group_by(Chr,Pos) %>% summarize(SumDevOLD = sum(dev),SumDevNew = sum(deviance), Count = n())
preComboAgg $newP = -log10(pchisq(preComboAgg $SumDevNew, preComboAgg$Count, lower.tail=F))
preComboAgg $oldP = -log10(pchisq(preComboAgg $SumDevOLD, preComboAgg$Count, lower.tail=F))
preComboAgg$Diff = abs(preComboAgg $newP - preComboAgg $oldP)
plot2 = ggplot(preComboAgg, aes(oldP, newP))+geom_point()+geom_abline(slope=1, intercept=0, color='blue')


m <- m[which(m $popAF > 3),]

mSub <- m %>% group_by(CHROM,POS) %>% summarize(SumDevNew = sum(deviance), Count = n())
newVSold = merge(mSub, m[,c('CHROM','POS','NewTest')], by=c('CHROM','POS'))
newVSold$OldTest = -log10(pchisq(newVSold $SumDevNew, newVSold $Count, lower.tail=F))
plot3 = ggplot(newVSold, aes(y= OldTest, x=NewTest))+geom_point()+geom_abline(slope=1, intercept=0, color='blue')+xlim(c(0,3.6))

chr22agg <- chr22 %>% group_by(Chr,Pos) %>% summarize(SumDevOLD = sum(dev), Count = n())
Mcombo = merge(chr22agg, mSub, by.x=c('Chr','Pos'), by.y=c('CHROM','POS'))
Mcombo $newP = -log10(pchisq(Mcombo $SumDevNew, 1, lower.tail=F))
Mcombo $oldP = -log10(pchisq(Mcombo $SumDevOLD, 1, lower.tail=F))
Mcombo$diff = abs(Mcombo$Count.x - Mcombo$Count.y)
plot1 = ggplot(Mcombo, aes(oldP, newP, color=diff))+geom_point()


newVSold = merge(m, preComboAgg, by.x=c('CHROM','POS'),by.y=c('Chr','Pos'))

plot4 = ggplot(newVSold, aes(y=oldP, x= NewTest))+geom_point()+geom_abline(slope=1, intercept=0, color='blue')+xlim(c(0,3.6))

plot_grid(plot2, plot3, ncol=1)

metaCount = apply(m, 1, function(c)sum(c>0.0001))
metaSum = rowSums(m)
metaDF = data.frame('metaSum'=metaSum, 'metaCount'=metaCount)
pos <- fread('~/genomes/genomes/hg19/Genotypes/Chr22_Pos.csv')
metaDF <-cbind(metaDF, head(pos,1000))

#metaDev = -log10(pchisq(metaSum, metaCount, lower.tail=F))
	
Pop.Qual1 = unlist(Pop.Qual)
PC.Qual1 = unlist(PC.Qual)
PopPC.Qual1 = unlist(PopPC.Qual)
PopPC.Q1 = unlist(PopPC.Q)

plt <-data.frame('metaDev' = -log10(pchisq(metaSum, metaCount, lower.tail=F)),
				'PopPCQ' = -log10(pchisq(PopPC.Q1, 1, lower.tail=F)),
				'popQual' = -log10(pchisq(Pop.Qual1, 1, lower.tail=F)),
				'PCQual' = -log10(pchisq(PC.Qual1, 1, lower.tail=F)), 
				'PopPCQual' = -log10(pchisq(PopPC.Qual1, 1, lower.tail=F)))
				
plt <- cbind(plt, m.Pos)
				
library(ggplot2)
library(cowplot)				
p1 = ggplot(plt, aes(x=plt$popQual, y=plt$PopPCQual))+geom_point()+xlim(c(0,20))+ylim(c(0,20))+geom_abline(slope=1, intercept=0, color='blue')
p2 = ggplot(plt, aes(x=plt$PCQual, y=plt$PopPCQual))+geom_point()+xlim(c(0,20))+ylim(c(0,20))+geom_abline(slope=1, intercept=0, color='blue')	
p3 = ggplot(plt, aes(x=plt$PopPCQ, y=plt$PopPCQual))+geom_point()+xlim(c(0,20))+ylim(c(0,20))+geom_abline(slope=1, intercept=0, color='blue')
p4 = ggplot(plt, aes(x=plt$SumDevNew.p, y=plt$PopPCQual))+geom_point()+xlim(c(0,20))+ylim(c(0,20))+geom_abline(slope=1, intercept=0, color='blue')

plot_grid(p1,p4,p2,p3)

########################
m2<-as.data.frame(meta2)
names(m2)=unique(samples$Pop)	
metaCount2 = apply(m2, 1, function(c)sum(c>0.0001))
metaSum2 = rowSums(m2)
metaDF2 = data.frame('metaSum'=metaSum2, 'metaCount'=metaCount2)
pos <- fread('~/genomes/genomes/hg19/Genotypes/Chr22_Pos.csv')
metaDF2 <-cbind(metaDF2, head(pos,1000))	



#####
#TEST META ANALYSIS

	singlePop = list() #list of deviances for each pop
	count = list() #allele frequency per pop
	j=1
	for (p in unique(samples$Pop)){
		print(p)
		PopNames = samples[which(samples$Pop == p),]$Name
		PopGeno <- chunk %>% select(PopNames)
		count[[j]] = rowSums(PopGeno) #get AF
		samp = samples[which(samples$Pop == p),]
		singlePop[[j]] = apply(PopGeno,1, function(x) 
					anova(
						glm(x ~
							samp$average_quality_of_mapped_bases,
						family=binomial),
					test="Chi")[2,2])	
		j=j+1
	}
	singlePop = as.data.frame(singlePop, col.names=paste(unique(samples$Pop),'.dev',sep=''))
	singlePop = cbind(chunk[,c(1:3)], singlePop)
	sMelt = melt(singlePop, id=c('CHROM','POS','ID'))
	count = as.data.frame(count, col.names=paste(unique(samples$Pop),'.freq',sep=''))
	count = cbind(chunk[,c(1:3)], count)
	cMelt = melt(count, id=c('CHROM','POS','ID'))
	sMelt$count = cMelt$value
	sMelt <- sMelt[which(sMelt $count > 3),]
	metaOut <- sMelt %>% group_by(CHROM,POS) %>% summarize(SumDev = sum(value), Count = n())
	metaOut$metaP = -log10(pchisq(metaOut$SumDev, metaOut$Count, lower.tail=F))
	
	meta[[i]]=metaOut	
					