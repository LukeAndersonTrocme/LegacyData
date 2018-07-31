##All Pops Logistic
library(ggplot2)
library(cowplot)

all<-read.table('~/Documents/Regression/CHR22.Regression_ALL.csv', header=T)
all$Chr<-as.numeric(as.character(all$Chr))
all$Pos<-as.numeric(as.character(all$Pos))
all$p<-as.numeric(as.character(all$p))

all<-all[complete.cases(all),]

ggplot(all, aes(x=Pos, y=p, color=Pop))+geom_point()+facet_wrap(~Pop, ncol=2)

JPT<-all[which(all$Pop=='JPT'),]
plots<-list()
for(pop in unique(all$Pop)){
sub<-all[which(all$Pop==pop),]
m<-merge(sub,JPT, by='Pos')
plots[[pop]]<-ggplot(m, aes(x=p.x, y=p.y))+geom_point()+labs(x=pop,y='JPT')
}


plot_grid(plots[['ACB']], plots[['ASW']], plots[['BEB']], plots[['CDX']], plots[['CEU']], plots[['CHB']], plots[['CHS']], plots[['CLM']], plots[['ESN']], plots[['FIN']], plots[['GBR']], plots[['GIH']], plots[['GWD']], plots[['IBS']], plots[['ITU']], plots[['KHV']], plots[['LWK']], plots[['MSL']], plots[['MXL']], plots[['PEL']], plots[['PJL']], plots[['PUR']], plots[['STU']], plots[['TSI']], plots[['YRI']])

ggsave('~/Desktop/PvalComparison.jpg',height=10,width=10)

agg<-aggregate(.~Pos, data=all, prod)