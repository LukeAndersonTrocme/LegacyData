##Compare GWAS to Logistic Regression
library(ggplot2)


GWAS<-read.table('~/Desktop/CHR22_GWAS_JPT.assoc.linear')
names(GWAS)<-c('Chr','Snp','Pos','GT','A','B','C','D','E','F','G','Pval')
GWAS$Pval<--log10(GWAS$Pval)
GWAS<-GWAS[,c('Chr','Pos','Pval')]
Reg<-read.table('~/Desktop/CHR22.Regression_JPT.csv', header=T)
Reg$p<-as.numeric(as.character(Reg$p))
both<-merge(GWAS, Reg, by=c('Chr','Pos'))

ggplot(both, aes(x=p,y=Pval))+geom_point()+theme_classic()+labs(x='Logistic Regression',y='GWAS',main='JPT Chromosome 22')

ggsave(file='~/Desktop/Chr22_JPT_GWAS_vs_Regression.jpeg',height=8,width=8)