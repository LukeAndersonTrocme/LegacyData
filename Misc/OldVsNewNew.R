library(data.table)
library(ggplot2)
library(cowplot)
###other


chr22= fread('/Users/luke/Desktop/Chr22_5PC_Regression.txt')
chr22$ID<-as.numeric(as.character(chr22 $ID))
chr10=fread('~/Desktop/CHR10.Regression.aa.txt')
chr22.glm2 = fread('/Users/luke/Desktop/Chr22_5PC_Regression.GLM2.txt')
gl = merge(chr22, chr22.glm2, by=c('CHROM','POS','ID'))
ggplot(gl, aes(dev.p.x, dev.p.y))+geom_point()+labs(x='GLM - 5PCs', y='GLM2 - 5PCs')


All.old= fread('~/Documents/QualityPaper/Misc/Reg_ALL_SITES_5col.txt')
All.old.10 <- All.old[which(All.old$Chr=="10"),]
All.old <- All.old[which(All.old$Chr=="22"),]

All.old.10 $Chr<-as.numeric(as.character(All.old.10 $Chr))
All.old.10 $Pos<-as.numeric(as.character(All.old.10 $Pos))
All.old.10 $SumDev <-as.numeric(as.character(All.old.10 $SumDev))
All.old.10 $Count<-as.numeric(as.character(All.old.10 $Count))
All.old.10 $old.P = -log10(pchisq(All.old.10 $SumDev, All.old.10 $Count,lower.tail=F))

All.old $Chr<-as.numeric(as.character(All.old $Chr))
All.old $Pos<-as.numeric(as.character(All.old $Pos))
All.old $SumDev <-as.numeric(as.character(All.old $SumDev))
All.old $Count<-as.numeric(as.character(All.old $Count))
All.old $old.P = -log10(pchisq(All.old $SumDev, All.old $Count,lower.tail=F))

plt = unique(merge(chr10, All.old.10, by.x=c('CHROM','POS'), by.y=c('Chr','Pos')))
plt $ID<-as.numeric(as.character(plt $ID))
plt=plt[which(complete.cases(plt)==T),]

p1=
ggplot(plt,aes(x=dev.p, y=old.P))+geom_point(shape=1)+geom_abline(slope=1,intercept=0, color='blue')+ggtitle('200,000 snps from Chr 10')+guides(color=F)+labs(x='5PC glm2 -log10(p)', y='meta analysis -log10(p)')

p2=ggplot(plt[which(plt$ID > 0.0026),],aes(x=dev.p, y=old.P, color=ID))+geom_point()+geom_abline(slope=1,intercept=0, color='blue')+ggtitle('Remove sites below 0.26%')+guides(color=F)

p3=ggplot(plt[which(plt$ID < 0.0026),],aes(x=dev.p, y=old.P, color = ID))+geom_point()+geom_abline(slope=1,intercept=0, color='blue')+ggtitle('The removed sites')+xlim(c(0,170))+ylim(c(0,130))+guides(color=F)

plot_grid(p1,p2,p3,nrow=1)

f1 = ggplot(plt, aes(x=ID, y=dev.p))+geom_point(shape=1)+geom_vline(xintercept=0.0026)+ggtitle('Model with 5PCs')+labs(x='Allele Frequency',y='-log10(p)')
f2 = ggplot(plt, aes(x=ID, y=old.P))+geom_point(shape=1)+geom_vline(xintercept=0.0026)+ggtitle('Meta Analysis')+labs(x='Allele Frequency',y='-log10(p)')
plot_grid(f1,f2)

theoretical <- 1:nrow(plt[which(plt$ID > 0.0026),])/(nrow(plt[which(plt$ID > 0.0026),])+1)

ggplot(plt[which(plt$ID > 0.0026),]) +geom_point(aes(x = -log10(theoretical), y = sort(plt[which(plt$ID > 0.0026),]$old.P, decreasing=T)), colour = "grey80", shape=1)+ geom_point(aes(x = -log10(theoretical), y = sort(plt[which(plt$ID > 0.0026),]$dev.p, decreasing=T)), colour = "blue", shape=1) +geom_abline(linetype=3)+theme_bw()+ylab('-log10(p-values)')
  
  
#####
GT = fread('~/Desktop/R18010772.temp.txt')
names(GT) = as.character(unlist(ColNames))
mGT = melt(GT, id=c('CHROM','POS','ID'))
mGT = merge(mGT, samples, by.x='variable', by.y='Name')

current = glm(mGT$value ~ mGT $Pop + mGT $PC1 + mGT $PC2 + mGT $average_quality_of_mapped_bases, family=binomial)
pred.current = predict(current)

full = glm(mGT$value ~ mGT $PC1 + mGT $PC2 + mGT $PC3 + mGT $PC4 + mGT$PC5 + mGT $average_quality_of_mapped_bases, family=binomial)
pred.full = predict(full)

full.c = glm(mGT$value ~ mGT $PC1 + mGT $PC2 + mGT $PC3 + mGT $PC4 + mGT$PC5 + mGT $average_quality_of_mapped_bases, family=binomial,control = list(maxit = 10000))
pred.c=predict(full.c)

PCflip = glm(mGT$value ~ mGT$PC5 + mGT $PC4 + mGT $PC3 + mGT $PC2 + mGT $PC1 + mGT $average_quality_of_mapped_bases , family=binomial)


logF = logistf(mGT$value ~ mGT $PC1 + mGT $PC2 + mGT $PC3 + mGT $PC4 + mGT$PC5 + mGT $average_quality_of_mapped_bases)