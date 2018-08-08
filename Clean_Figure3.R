#load libraries
library(ggplot2)
library(cowplot)
library(data.table)
#pop names, used to read/write files
pops<-read.table('/Users/luke/Documents/MutSpect/WrapUpPlots/NamePop.txt')
#read colors per pop (used in 1kG Nature)
colors<-read.table("~/Dropbox/LukeTemp/NamePopBigPop_color.html.txt", comment.char='>',col.names=c('Pop','hex','BigPop'))
MyColour <- as.character(colors $hex)
names(MyColour) <- colors $Pop

#read the list of sig snps
#bad=read.table('~/Documents/1FinalGWAS/1kGP_GenomeWide.4bed_filtered_sig6.POS.txt',col.names=c('Chr', 'Pos'))
dir='/Volumes/gravel-1/luke_projects/1000Genomes/Regression/'
SFS=data.frame()
for(f in seq(1,nrow(pops))){
pop=pops[f,]
print(as.character(pop))
fileNames = list.files(path=dir, pattern=paste(pop), full.names = T)
Reg = do.call(rbind, lapply(fileNames, function(x) fread(x)))
print('Read')
Reg<-Reg[which(Reg$Pos != 'Pos'),]
Reg$Chr<-as.numeric(as.character(Reg$Chr))
Reg$Pos<-as.numeric(as.character(Reg$Pos))
Reg$dev<-as.numeric(as.character(Reg$dev))
Reg$p=NULL
print('Clean')
Reg$plog10 <- - pchisq(Reg$dev, 1, lower.tail=F, log.p=T)/log(10)
sig.6<-Reg[which(Reg$plog10 > 6),]
rm(Reg)
write.table(sig.6[,c('Chr','Pos','plog10', 'Pop')], paste('~/Documents/Regression/',pop,'_sig6.csv', sep=','))
#SFS=rbind(SFS,sig.6[,c('Chr','Pos','plog10', 'Pop')])
}


SFS<-unique(SFS)
NumPop<-as.data.frame(table(SFS$Chr,SFS$Pos))
UnPop<-NumPop[which(NumPop$Freq==1),]
NumPop<-NumPop[which(NumPop$Freq>1),]
m=merge(SFS,NumPop, by.y=c('Var1','Var2'), by.x=c('Chr','Pos'))
n=as.data.frame(table(m$Pop))
plt=merge(m,n,by.x='Pop',by.y='Var1')
plt <- plt[with(plt, order(-Freq.x, as.numeric(Chr),as.numeric(Pos))),]
plt$Pos <- factor(plt$Pos, levels = unique(plt$Pos))
plt$nrow<-seq(length=nrow(plt))
plt$Pcat<-4
plt[which(plt$Plog10 <8),]$Pcat<-2
plt[which((plt$Plog10 >8) & (plt$Plog10 <10)),]$Pcat<-3

xaxis<-plt[!duplicated(plt[,c('Freq.x')]), ]
xaxis<-xaxis[which(xaxis$Freq.x<11),]
xaxis<-rbind(xaxis,tail(plt,1))

p1=ggplot(plt, aes(x=Pos,y= reorder(Pop,-Freq.y), color=Pop, size=Pcat))+geom_point(shape=3)+theme_classic()+scale_size(breaks=c(2,3,4), range=c(1,3),labels=c('6-8','8-10','10+'), guide = guide_legend(direction = "horizontal"), name='-log10(p)')+scale_color_manual(breaks= plt$Pop,values = MyColour)+ylab('Populations')+ggtitle('Overlap of SNPs found to be Significant in Quality-GWAS')+theme(legend.position = c(0.85,0.95),legend.box.background = element_rect(colour = "black"), plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.margin = unit(c(1,1,0,1), "cm"))+guides(color=F)

p2=ggplot(plt, aes(x=Pos,y= Freq.x, group=1))+geom_hline(yintercept=1, color='grey', linetype=2)+geom_hline(yintercept=5, color='grey', linetype=2)+geom_hline(yintercept=10, color='grey', linetype=2)+geom_line()+theme_classic()+theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0,1,1,1), "cm"))+labs(x='SNPs ranked by frequency of occurence and genomic position', y='Frequency\nof Occurence')+scale_y_continuous(limits=c(0,16),breaks=c(1,5,10,15))+scale_x_discrete(expand=c(0.01,0),breaks=xaxis$Pos, labels=xaxis$nrow)

p3=plot_grid(p1, p2, ncol=1, align="v", rel_heights = c(3, 1))
ggsave('~/Documents/QualityPaper/SNPOverlap6.jpg',p3, height=8, width=8)


#+geom_vline(xintercept = which(levels(plt$Var1) %in% '3358271'),color='grey', linetype=2)
#+geom_vline(xintercept = which(levels(plt$Var1) %in% '3358271'),color='grey', linetype=2)

UnPop<-NumPop[which(NumPop$Freq==1),]
m1=merge(SFS, UnPop, by.y=c('Var1','Var2'), by.x=c('Chr','Pos'))
n1=as.data.frame(table(m1$Pop))
plt1=merge(m1,n1,by.x='Pop',by.y='Var1')
plt1 <- plt1[with(plt1, order(as.numeric(Chr),as.numeric(Pos))),]
plt1$Pos <- factor(plt1$Pos, levels = unique(plt1$Pos))
xaxis<-plt1[!duplicated(plt1[,c('Chr')]), ]

ggplot(plt1, aes(x=Pos,y= reorder(Pop,-Freq.y), color=Pop, alpha=0.02*Plog10))+geom_point(shape=1)+theme_classic()+scale_color_manual(breaks= plt$Pop,values = MyColour)+labs(title='SNPs found to be Significant in Quality-GWAS',subtitle='40,384 SNPs included in this plot', y='Populations',x='Chromosomes')+theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),plot.margin = unit(c(1,1,1,1), "cm"))+guides(color=F)+scale_x_discrete(expand=c(0.01,0), breaks= xaxis$Pos, labels=xaxis$Chr)
ggsave('/Users/luke/Documents/GWAS_Qual/SNP_Uni_p6.tiff',p3, height=7, width=14)