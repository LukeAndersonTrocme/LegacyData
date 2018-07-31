##Deviations
library(dplyr)
library(data.table)
library(ggplot2)

Reg<-fread('/Users/luke/Documents/Regression/CHR19-22.Regression.csv')
Reg<-Reg[which(Reg$Pos != 'Pos'),]
Reg<-Reg[which(Reg$Chr != '14'),]
Reg$Chr<-as.numeric(as.character(Reg$Chr))
Reg$Pos<-as.numeric(as.character(Reg$Pos))
Reg$dev<-as.numeric(as.character(Reg$p))
Reg$Pop=NULL
Reg$p=NULL

Meta <- Reg %>% group_by(Chr,Pos) %>% summarize(SumDev = sum(dev), MeanDev = mean(dev), Count = n())

Meta$plog10 <- - pchisq(Meta$SumDev, Meta$Count, lower.tail=F, log.p=T)/log(10)

ggplot(Meta, aes(x=Pos, y= plog10, color=Chr))+geom_point(shape=1)+facet_grid(~Chr, scales='free_x', space= 'free_x')+theme_classic()+guides(color=F)+geom_hline(yintercept=8, color='blue')
ggsave('/Users/luke/Documents/Regression/Manhattan_Chr19-22.jpg', height=8, width=11)

ggplot(Meta[which(Meta$Count > 1),], aes(x=Pos, y= MeanDev, color=Chr))+geom_point(shape=1)+facet_grid(Count~Chr, scales='free_x', space= 'free_x')+theme_classic()+guides(color=F)
ggplot(Meta[which((Meta$Count > 1)&(Meta$MeanDev > 10)),], aes(x=Pos, y= MeanDev, color=Chr))+geom_point(shape=1)+facet_grid(Count~Chr, scales='free_x', space= 'free_x')+theme_classic()+guides(color=F)

Meta[which((Meta$Count > 2)&(Meta$plog10 > 20)),]




gwas<-fread('~/Documents/GWAS_Qual/Meta/INT_19.22.meta')
#gwas<-fread('/Users/luke/Documents/GWAS_Qual/Meta/NoNormal_19.22.meta')
names(gwas)<-c('Chr','Pos','SNP','A1','A2','N','P','P(R)','BETA','BETA(R)','Q','I')

gwas$log10P<--log10(gwas$P)
gwas$log10PR<--log10(gwas$'P(R)')
gwas<-gwas[,c('Chr','Pos','log10P','log10PR')]
head(combo)

combo<-merge(Meta[,c('Chr','Pos','plog10')], gwas, by=c('Chr','Pos'))
ggplot(combo, aes(x=plog10, y=log10P))+geom_point(shape=1)+labs(x='Logistic Regression', y='GWAS Meta Analysis')+theme_classic()+geom_vline(xintercept=8, color='blue')+geom_hline(yintercept=8, color='blue')
ggsave(file='/Users/luke/Documents/Regression/Logistic_GWAS_pval.jpg',height=6,width=8)

L<-combo[which((combo$plog10 > 8)&(combo$log10P < 8)),]

G<-combo[which((combo$plog10 < 8)&(combo$log10P > 8)),]

GL<-combo[which((combo$plog10 > 8)&(combo$log10P > 8)&(combo$plog10 <10)&(combo$log10P < 10)),]

topGL<-read.table('/Users/luke/Documents/Regression/CHR19_GL.txt',sep=',')
ggplot(topGL, aes(x=V6, y=V4))+geom_point()+facet_grid(V5~V2)+geom_smooth(method = "glm",     method.args = list(family = "binomial"), se = FALSE,size=0.5)+theme_classic()
ggsave('/Users/luke/Documents/Regression/CHR19_GL.jpg')

topL<-read.table('/Users/luke/Documents/Regression/CHR19_L.txt',sep=',')
ggplot(topL, aes(x=V6, y=V4))+geom_point()+facet_grid(V5~V2)+geom_smooth(method = "glm",     method.args = list(family = "binomial"), se = FALSE,size=0.5)+theme_classic()
ggsave('/Users/luke/Documents/Regression/CHR19_L.jpg')

topG<-read.table('/Users/luke/Documents/Regression/CHR19_G.txt',sep=',')
ggplot(topG, aes(x=V6, y=V4))+geom_point()+facet_grid(V5~V2)+geom_smooth(method = "glm",     method.args = list(family = "binomial"), se = FALSE,size=0.5)+theme_classic()
ggsave('/Users/luke/Documents/Regression/CHR19_G.jpg')

topL20<-read.table('/Users/luke/Documents/Regression/CHR20_L.txt',sep=',')
ggplot(topL20, aes(x=V6, y=V4))+geom_point()+facet_grid(V5~V2)+geom_smooth(method = "glm",     method.args = list(family = "binomial"), se = FALSE,size=0.5)+theme_classic()
ggsave('/Users/luke/Documents/Regression/CHR20_L.jpg')