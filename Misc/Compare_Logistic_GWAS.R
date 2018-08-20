##Compare GWAS to Logistic Regression
library(ggplot2)


GWAS<-read.table('~/Documents/Regression/Chr22_GBR_GWAS.assoc.linear')
names(GWAS)<-c('Chr','Snp','Pos','GT','V1','V2','V3','V4','V5','V6','V7','Pval')
GWAS$Pval<--log10(GWAS$Pval)
GWAS<-GWAS[,c('Chr','Pos','Pval')]
Reg<-read.table('~/Documents/Regression/CHR22.Regression_GBR.csv', header=T)
Reg$p<-as.numeric(as.character(Reg$p))
both<-merge(GWAS, Reg, by=c('Chr','Pos'))

ggplot(both, aes(x=p,y=Pval))+geom_point()+theme_classic()+labs(x='Logistic Regression',y='GWAS',main='GBR Chromosome 22')

ggsave(file='~/Desktop/Chr22_JPT_GWAS_vs_Regression.jpeg',height=8,width=8)