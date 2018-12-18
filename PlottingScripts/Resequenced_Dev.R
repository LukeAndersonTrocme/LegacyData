library(data.table)
library(ggplot2)
library(cowplot)

CHB.kgp<-fread('/Users/luke/Documents/Regression/GenomeWideRegression_CHB.csv')
CHS.kgp<-fread('/Users/luke/Documents/Regression/GenomeWideRegression_CHS.csv')

CHB.kgp$Pos <- as.numeric(as.character(CHB.kgp$Pos))
CHS.kgp$Pos <- as.numeric(as.character(CHS.kgp $Pos))
CHB.kgp$Chr <- as.numeric(as.character(CHB.kgp$Chr))
CHS.kgp$Chr <- as.numeric(as.character(CHS.kgp $Chr))
CHB.kgp$dev <- as.numeric(as.character(CHB.kgp$dev))
CHS.kgp$dev <- as.numeric(as.character(CHS.kgp $dev))

CHB.kgp<-CHB.kgp[complete.cases(CHB.kgp),]
CHS.kgp<-CHS.kgp[complete.cases(CHS.kgp),]


CHB.han<-fread('~/genomes/genomes/Han90/90_Han_Chinese_GenomeWide_CHB_Regression.txt')
CHS.han<-fread('~/genomes/genomes/Han90/90_Han_Chinese_GenomeWide_CHS_Regression.txt')

CHB.han$Chr <- as.numeric(as.character(gsub('chr','',CHB.han$Chr)))
CHS.han$Chr <- as.numeric(as.character(gsub('chr','',CHS.han$Chr)))
CHB.han$Pos <- as.numeric(as.character(CHB.han$Pos))
CHS.han$Pos <- as.numeric(as.character(CHS.han$Pos))
CHB.han$dev <- as.numeric(as.character(CHB.han$dev))
CHS.han$dev <- as.numeric(as.character(CHS.han$dev))

CHB.han<-CHB.han[complete.cases(CHB.han),]
CHS.han<-CHS.han[complete.cases(CHS.han),]

CHB<-merge(CHB.han, CHB.kgp, by=c('Chr', 'Pos', 'Pop'), all =T)
CHS<-merge(CHS.han, CHS.kgp, by=c('Chr', 'Pos', 'Pop'), all= T)
CHB[is.na(CHB)] <- 0
CHS[is.na(CHS)] <- 0

CHS$p.han <- -log10(pchisq(CHS$dev.x, 1, lower.tail=F))
CHS$p.kgp <- -log10(pchisq(CHS$dev.y, 1, lower.tail=F))

CHB$p.han <- -log10(pchisq(CHB$dev.x, 1, lower.tail=F))
CHB$p.kgp <- -log10(pchisq(CHB$dev.y, 1, lower.tail=F))

p1=ggplot(CHB, aes(x=p.han, y=p.kgp))+geom_point(shape=1)+theme_classic()+labs(x='Resequence',y='1kGP', main='CHB')+theme(plot.title = element_text(hjust = 0.5))
p2=ggplot(CHS, aes(x=p.han, y=p.kgp))+geom_point(shape=1)+theme_classic()+labs(x='Resequence',y='1kGP', main='CHB')+theme(plot.title = element_text(hjust = 0.5))

tot=plot_grid(p1,p2, ncol=1)

ggsave('~/Documents/QualityPaper/Figures/ResequencePvals.jpg', height=10, width=5)

