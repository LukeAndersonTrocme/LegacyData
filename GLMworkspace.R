#########WorkSpace
library(lme4)
library(lattice)

pos5<-head(unique(GT.Qual$Pos),10)
jj<-GT.Qual[which((GT.Qual$Pos=='7262382')&(GT.Qual$Pop=='JPT')),]
jj<-GT.Qual[which(GT.Qual$Pos=='7262382'),]
j5<-GT.Qual[which(GT.Qual$Pos %in% pos5),]
j5$Pos<-as.factor(as.numeric(j5$Pos))
j5$value<-gsub(2,1,j5$value)


glm.fit<-glm(as.numeric(value)-1~ average_quality_of_mapped_bases+Pop, data=jj, family=binomial)

lmer.fit<-lme4::glmer(as.numeric(value)-1~ average_quality_of_mapped_bases+(average_quality_of_mapped_bases|Pop), data=jj, family=binomial)

lmer.fit1<-lme4::glmer(as.numeric(value)-1~ average_quality_of_mapped_bases+(1|Pop), data=jj, family=binomial)

glmer.all<-lme4::glmer(as.numeric(value) ~ average_quality_of_mapped_bases + (1|Pos) +(1|Pop), data=GT.Qual, family=binomial)
str(glmer1 <- ranef(glmer.all, condVar = TRUE))

byPop<-glmer1$Pop
byPop$Pop<-rownames(byPop)
names(byPop)<-c('Intercept','Pop')

byPos<-glmer1$Pos
byPos$Pos<-rownames(byPos)
names(byPos)<-c('Intercept','Pos')

p1=ggplot(byPos,aes(x=Intercept,y=reorder(Pos,-Intercept)))+geom_vline(xintercept=0, color='grey70')+geom_point(shape=1)+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+labs(y='SNPs')+ggtitle('Generalized linear mixed model (logit)\nRandom Effects')

p2=ggplot(byPop,aes(x=Intercept,y=reorder(Pop,-Intercept)))+geom_vline(xintercept=0, color='grey70')+geom_point(shape=1)+labs(y='Populations')+ggtitle('Generalized linear mixed model (logit)\nRandom Effects')
sub<-glmer.pos[which(glmer.pos$'(Intercept)'< -3),]
subSub<-merge(sub,GT.Qual, by='Pos')
subSub$Title<-paste(subSub$Pos,'\n', round(subSub$'(Intercept)', digits=3))

sub1<-ggplot(subSub[which(subSub$Pop=='JPT'),], aes(x= average_quality_of_mapped_bases, y=as.numeric(value)))+geom_point()+stat_smooth(method="glm", se=F, method.args = list(family="binomial"))+facet_wrap(~Title)+theme_classic()+scale_y_continuous(breaks=c(0,1))+labs(x='Quality',y='Genotype')

more<-glmer.pos[which(glmer.pos$'(Intercept)'>8),]
subMore<-merge(more,GT.Qual, by='Pos')
subMore $Title<-paste(subMore $Pos,'\n', round(subMore $'(Intercept)', digits=3))

sub2<-ggplot(subMore[which(subMore $Pop=='JPT'),], aes(x= average_quality_of_mapped_bases, y=as.numeric(value)))+geom_point()+stat_smooth(method="glm", se=F, method.args = list(family="binomial"))+facet_wrap(~Title)+theme_classic()+scale_y_continuous(breaks=c(0,1))+labs(x='Quality',y='Genotype')

pf<-plot_grid(p1,p2,sub1,sub2,nrow=2)
ggsave('~/Documents/Regression/JPT_sig6_GLMM.jpg')


glmer.pos<-ranef(glmer.all)$Pos
glmer.pos$Pos<-row.names(glmer.pos)


lmer<-ranef(lmer.fit)$Pop
lmer$Pop<-rownames(lmer)
names(lmer)=c('Int','Qual')
ggplot(lmer,aes(x=Int,y=Qual))+geom_point()
ggplot(lmer,aes(y=Qual,x=Pop))+geom_point()

str(rr1 <- ranef(lmer.fit1, condVar = TRUE))
dotplot(rr1)

str(rr2 <- ranef(lmer.fit1, condVar = TRUE))
dotplot(rr2)



ggsave('~/Documents/QualityPaper/GLM_binomial_perPop.jpg',heigh=10,width=5)


ggplot(jj, aes(x= average_quality_of_mapped_bases, y=reorder(variable, average_quality_of_mapped_bases)))+geom_point()+facet_wrap(~Pop, ncol=2)+theme_classic()+labs(x='Quality',y='Genotype')+theme(axis.text=element_blank())
#########WorkSpace


########HERE

##Postion is Covariate
jpt.yul.glm<-glm(value ~ average_quality_of_mapped_bases + PopPos, data= JPT.YUL, family=binomial)

coef.jptyul<-tail(as.data.frame(coef(summary(jpt.yul.glm))[1,4]),nrow(as.data.frame(coef(jpt.yul.glm)))-2)
coef.jptyul$Pos<-gsub('PopPos','',rownames(coef.jptyul))
names(coef.jptyul)<-c('P','PopPos')
coef.jptyul$Plog10<--log10(coef.jptyul$P)
jptYul<-merge(JPT.YUL[,c(Pop)], coef.jptyul, by='PopPos')
jptYul<-melt(jptYul,)
ggplot(jptYul, aes(x=Coef, y=reorder(Pos,as.numeric(Coef))))+geom_point()+theme_classic()+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())



top5<-JPT.GT[which((GT.Qual$variable %in% head(unique(GT.Qual$variable),2))&(GT.Qual$Pos %in% head(unique(GT.Qual$Pos),5))),]
top5<-droplevels(top5)
rownames(top5)<-paste(top5$variable,top5$Pos, sep='-')

model.matrix(~-1+average_quality_of_mapped_bases + Pos,top5)

##OFFSET
#jpt.glm<-glm(value ~ average_quality_of_mapped_bases + offset(AF), data= JPT.GT, family=binomial)

##Postion is Covariate
jpt.glm<-glm(value ~ average_quality_of_mapped_bases + PopPos, data= JPT.GT, family=binomial)

##Position Covariate + Populations


coef.jpt<-tail(as.data.frame(coef(jpt.glm)),nrow(as.data.frame(coef(jpt.glm)))-2)
coef.jpt$Pos<-rownames(coef.jpt)
names(coef.jpt)<-c('Coef','Pos')
ggplot(coef.jpt, aes(x=Coef, y=reorder(Pos,as.numeric(Coef))))+geom_point()+theme_classic()+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())


JPT.GT<-head(GT.Qual[which(GT.Qual$Pop == 'JPT'),],200)
JPT.YUL<-droplevels(GT.Qual[which((GT.Qual$Pos %in% head(unique(GT.Qual$Pos),200)) & ((GT.Qual$Pop == 'JPT')|(GT.Qual$Pop == 'YRI'))),])
JPT.YUL$PopPos<-paste(JPT.YUL$Pos,JPT.YUL$Pop,sep=':')