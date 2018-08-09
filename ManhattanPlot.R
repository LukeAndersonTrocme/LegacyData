##Manhattan
library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(multtest)

dir='/Volumes/gravel/luke_projects/1000Genomes/Regression/'
out='/Volumes/gravel/luke_projects/1000Genomes/MeanDev/'
if(F){
for(chrom in seq(1,22)){
	print(chrom)
	fileNames = list.files(path=dir,pattern=paste("CHR",chrom,".Regression",sep=''), full.names = T)
	if(length(fileNames) != 0){
	Reg = do.call(rbind, lapply(fileNames, function(x) read.table(x, header=T, sep=' ')))
	
	print(head(Reg))
	
	Reg<-Reg[which(Reg$Pos != 'Pos'),]
	Reg$Chr<-as.numeric(as.character(Reg$Chr))
	Reg$Pos<-as.numeric(as.character(Reg$Pos))
	Reg$dev<-as.numeric(as.character(Reg$dev))
	
	Reg <- Reg %>% group_by(Chr,Pos) %>% summarize(SumDev = sum(dev), MeanDev = mean(dev), Count = n())
	
	write.table(Reg,paste(out,chrom,'_MeanDev.csv',sep=''))
	}
}
}
fileNames = list.files(path= out,pattern='*_MeanDev.csv', full.names = T)
	
Reg = do.call(rbind, lapply(fileNames, function(x) read.table(x, header=T)))
Reg$plog10 <- - pchisq(Reg$SumDev, Reg$Count, lower.tail=F, log.p=T)/log(10)
Reg$p <- pchisq(Reg$SumDev, Reg$Count, lower.tail=F)

sig.6<-Reg[which((Reg$plog10 >= 6)&(Reg$plog10 < 20)&(Reg$Count > 1)),]
sig.20<-Reg[which((Reg$plog10 >= 20)&(Reg$Count > 1)),]
sig.20$plog10=20
NotSig<-Reg[which((Reg$plog10 < 6)&(Reg$Count > 1)),]


ggplot(NotSig, aes(x=Pos, y = plog10, color=as.factor(Chr)))+
facet_grid(~Chr, scales='free_x', space='free_x', switch='x')+
scale_fill_manual (values=getPalette(colourCount))+
scale_y_continuous(expand=c(0,0))+
geom_point(alpha=0.3, size=1)+theme_classic()+
geom_point(data=sig.6, aes(x=Pos, y = plog10,),color='black',shape=3)+
geom_point(data=sig.20, aes(x=Pos, y = plog10,),color='black',shape=1)+
labs(y='-log10(p)', x='Chromosome')+
theme(plot.title = element_text(hjust = 0.5), 
plot.subtitle = element_text(hjust = 0.5), 
axis.text.x=element_blank(), 
axis.ticks.x=element_blank(),
axis.line.x=element_blank(), 
strip.background = element_blank(),
strip.text.x = element_text(size = 6))+
guides(color=F)+expand_limits(y=c(0,21))
ggsave('~/Documents/QualityPaper/ManhattanPlot.jpg', height=5, width=10)

##RICK Adjustment
procedures <- c( "TSBH")
adjusted <- mt.rawp2adjp(Reg$p, procedures)
adj<-as.data.frame(adjusted$adj[order(adjusted$index), ])
Reg$adjusted<--log10(adj$TSBH_0.05)
Reg$adjustedP<-adj$TSBH_0.05
adj.6<-Reg[which((Reg$adjusted >= -log10(0.05))&(Reg$adjusted <20)&(Reg$Count > 1)),]
adj.20<-Reg[which((Reg$adjusted >= 20)&(Reg$Count > 1)),]
adj.20$adjusted=20
adj.ns<-Reg[which((Reg$adjusted < -log10(0.05))&(Reg$Count > 1)),]


ggplot(adj.ns, aes(x=Pos, y = adjusted, color=as.factor(Chr)))+
facet_grid(~Chr, scales='free_x', space='free_x', switch='x')+
scale_fill_manual (values=getPalette(colourCount))+
scale_y_continuous(expand=c(0,0))+
geom_point(alpha=0.3, size=1)+theme_classic()+
geom_point(data=adj.6, aes(x=Pos, y = adjusted,),color='black',shape=3)+
geom_point(data=adj.20, aes(x=Pos, y = adjusted,),color='black',shape=1)+
labs(y='adjusted -log10(p)', x='Chromosome')+
theme(plot.title = element_text(hjust = 0.5), 
plot.subtitle = element_text(hjust = 0.5), 
axis.text.x=element_blank(), 
axis.ticks.x=element_blank(),
axis.line.x=element_blank(), 
strip.background = element_blank(),
strip.text.x = element_text(size = 6))+
guides(color=F)+expand_limits(y=c(0,21))
ggsave('~/Documents/QualityPaper/ManhattanPlot_adjusted.jpg', height=5, width=10)


gwas<-fread('~/Documents/GWAS_Qual/Meta/INT_.meta')
names(gwas)<-c('Chr','Pos','SNP','A1','A2','N','P','P(R)','BETA','BETA(R)','Q','I')
gwas$log10P<--log10(gwas$P)
gwas$log10PR<--log10(gwas$'P(R)')
gwas<-gwas[,c('Chr','Pos','log10P','log10PR')]


combo<-merge(Reg[,c('Chr','Pos','plog10')], gwas, by=c('Chr','Pos'))

ggplot(combo, aes(x=plog10, y=log10P))+geom_point(shape=1)+labs(x='Logistic Regression', y='GWAS Meta Analysis')+theme_classic()+geom_vline(xintercept=8, color='blue')+geom_hline(yintercept=8, color='blue')
ggsave(file='/Users/luke/Documents/Regression/Logistic_GWAS_pval.jpg',height=6,width=8)


write.table(sig.6[,c('Chr','Pos')],file='/Users/luke/Documents/Regression/Sig.6.Pos', quote=F, col.names=F, row.names=F, sep='\t')


t<-Reg[which((Reg$Count>3)&(Reg$Chr==22)&(Reg$p<0.1)),]
t<-t[order(t$p),]
t$index <- seq(length=nrow(t))
ggplot(t, aes(x=index, y=p))+geom_point(shape=1)+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+geom_line(aes(x=index,y=0.1/nrow(t)), color='blue')+xlim(c(0,10000))


ggplot(Reg[which(Reg$p<0.01),], aes(x=p))+geom_histogram(bins=500)+theme_classic()
ggplot(Reg[which(Reg$p>0.05),], aes(x=p))+geom_histogram(bins=1000)+theme_classic()

ggplot(Reg, aes(x=SumDev))+geom_histogram()+facet_wrap(~Count,ncol=1)

ggplot(Reg[which((Reg$Count==26)),], aes(x=SumDev))+geom_histogram(bins=100)+xlim(c(0,50))

gwasCat<-fread('/Users/luke/Documents/GWAS_Qual/gwas_catalog_v1.0-associations_e91_r2018-03-13_cols.txt')
gwasCat$CHR_ID<-as.numeric(as.character(gwasCat$CHR_ID))
gwasCat$CHR_POS<-as.numeric(as.character(gwasCat$CHR_POS))
Pub<-merge(gwasCat, sig.6, by.x=c('CHR_ID','CHR_POS'), by.y=c('Chr','Pos'))


adjusted<-