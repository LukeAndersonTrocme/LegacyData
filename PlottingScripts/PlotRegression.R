library(ggplot2)

#pop names, used to read/write files
pops<-read.table('/Users/luke/Documents/MutSpect/WrapUpPlots/NamePop.txt')
#read colors per pop (used in 1kG Nature)
colors<-read.table("~/Dropbox/LukeTemp/NamePopBigPop_color.html.txt", comment.char='>',col.names=c('Pop','hex','BigPop'))
MyColour <- as.character(colors $hex)
names(MyColour) <- colors $Pop

snp<-read.table('~/Documents/QualityPaper/Misc/31384880.txt',sep=',')
dev<-read.table('~/Documents/QualityPaper/Misc/31384880.dev.txt')
dev$p<-- pchisq(dev$V4, 1, lower.tail=F, log.p=T)/log(10)
sig<-merge(snp, dev, by.x=c('V1','V2','V5'), by.y=c('V1','V2','V3'))
sig$Name<-paste(sig$V5, '\n -log10(p) :',round(sig$p, digits=2),'\n deviance :',round(sig$V4.y, digits=2))
sig<-sig[order(sig$p),]
sig$Name<-factor(sig$Name, levels=unique(sig$Name))

ggplot(sig, aes(x=V6, y=V4.x, color=V5))+geom_point(shape=1)+facet_wrap(~Name,ncol=1,strip.position='right')+geom_smooth(method = "glm",     method.args = list(family = "binomial"), se = FALSE,size=0.5)+theme_classic()+labs(x='Average quality of mapped bases', y='Presence of Alternate Allele')+scale_y_continuous(breaks=c(0,1))+scale_color_manual(breaks= sig$V5,values = MyColour)+guides(color=F)+theme(strip.text.y = element_text(angle = 0))

ggsave('~/Documents/QualityPaper/Figures/RegressionPlot.jpg',height=5,width=5)
ggsave('~/Documents/QualityPaper/Figures/RegressionPlot.tiff',height=5,width=5)


sig6<-read.table('~/Documents/QualityPaper/Significant6SNPs.txt', header=T)

snp<-read.table('~/Documents/QualityPaper/Misc/33264017.txt',sep=',', fill=T)
dev<-read.table('~/Documents/QualityPaper/Misc/33264017.dev.txt')
dev$p<-- pchisq(dev$V4, 1, lower.tail=F, log.p=T)/log(10)
dev<-dev[which(dev$p>6),]
sig<-merge(snp, dev, by.x=c('V2','V5'), by.y=c('V2','V3'))
sig$Name<-paste(sig$V5, '\n -log10(p) :',round(sig$p, digits=2),'\n deviance :',round(sig$V4.y, digits=2))
sig<-sig[order(sig$p),]
sig$Name<-factor(sig$Name, levels=unique(sig$Name))

ggplot(sig, aes(x=V6, y=V4.x, color=V5))+geom_point(shape=1)+facet_wrap(~Name,ncol=1, strip.position='right')+geom_smooth(method = "glm",     method.args = list(family = "binomial"), se = FALSE,size=0.5)+theme_classic()+labs(x='Average quality of mapped bases', y='Presence of Alternate Allele')+scale_y_continuous(breaks=c(0,1))+theme(strip.text.y = element_text(angle = 0))+geom_vline(xintercept=30, color='grey60',linetype=3)+scale_color_manual(breaks= sig$V5,values = MyColour)+guides(color=F)

ggsave('~/Documents/QualityPaper/Figures/RegressionPlot_mostSig2.jpg',height=11,width=8)
ggsave('~/Documents/QualityPaper/Figures/RegressionPlot_mostSig2.tiff',height=11,width=8)

weird<-fread('~/Documents/QualityPaper/Misc/LowAFsig_Format.txt')
ggplot(weird, aes(x=V6, y=V4))+facet_grid(.~V2)+geom_point()+coord_flip()+geom_smooth(method = "glm",     method.args = list(family = "binomial"), se = FALSE,size=0.5)+theme(axis.text.x=element_blank())+labs(x='Genotypes', y='Quality')


snp<-read.table('/Users/luke/Documents/QualityPaper/sig/57694260.format.txt',sep=',')
dev<-read.table('/Users/luke/Documents/QualityPaper/sig/57694260.deviance.txt')
dev$p<-- pchisq(dev$V4, 1, lower.tail=F, log.p=T)/log(10)
dev<-dev[which(dev$p>6),]
sig<-merge(snp, dev, by.x=c('V1','V2','V5'), by.y=c('V1','V2','V3'))
sig$Name<-paste(sig$V5, '\n -log10(p) :',round(sig$V4.y, digits=2))

ggplot(snp, aes(x=V6, y=V4, color=V5))+geom_point(shape=1)+facet_wrap(~V5,nrow=3)+geom_smooth(method = "glm",     method.args = list(family = "binomial"), se = FALSE,size=0.5)+theme_classic()+labs(x='Quality', y='Genotype')+scale_y_continuous(breaks=c(0,1))+geom_vline(xintercept=30, color='grey60',linetype=3)+scale_color_manual(breaks= snp$V5,values = MyColour)+guides(color=F)

ggsave('~/Documents/QualityPaper/RegressionPlot.jpg',height=5,width=5)