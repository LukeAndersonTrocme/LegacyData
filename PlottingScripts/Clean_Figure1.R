library(ggplot2)
library(colorspace)
library(cowplot)
library(data.table)

if(F){
#f1kGP='~/genomes/genomes/1kGP_4bedFiltered/1kGP_GenomeWide.4bed_filtered.freq.frq'
f1kGP='/Users/luke/genomes/genomes/1kGP_4bedFiltered/1kGP_GenomeWide_JPT_1.frq'
fnag='~/genomes/genomes/NAG_4bedFiltered/NAG_GenomeWide.4bed_filtered.freq.frq'
hweFile='/Users/luke/Documents/GWAS_Qual/Final_GWAS/1kGP_GenomeWide_JPT_sig6.hwe'

#Allele Frequency files
JPT<-fread(f1kGP, fill=T, col.names=c('CHROM','POS','N_ALLELES','N_CHR','JPT_AF','JPT_MAF'), colClasses=c('numeric'))
NAG<-fread(fnag, fill=T,col.names=c('CHROM','POS','N_ALLELES','N_CHR','NAG_AF','NAG_MAF'), colClasses=c('numeric'))

#merge them together to get joint frequency spectrum
join<-merge(JPT[which(JPT$JPT_AF <= 1),],
 			NAG[which(NAG$NAG_AF <= 1),], 
 			by=c('CHROM','POS'), all=T)
join[is.na(join)]<-0
rm(NAG)
rm(JPT)
#REMOVE  HWE SITES
NAG_hwe<-fread('/Users/luke/genomes/genomes/NAG_4bedFiltered/NAG_GenomeWide.4bed_filtered_sig6.freq.hwe')
JPT_hwe <- fread('/Users/luke/genomes/genomes/1kGP_4bedFiltered/1kGP_GenomeWide.4bed_filtered_sig6.freq.hwe')
hwe=merge(NAG_hwe[,c('V1','V2')], JPT_hwe[,c('V1','V2')], by=c('V1','V2'), all=T)

join= join[! paste(join$CHROM, join$POS) %in% c(paste(hwe$V1,hwe$V2)), ]
join=join[,c('CHROM','POS','JPT_MAF','JPT_AF','NAG_MAF','NAG_AF')]

join$CHROM<-as.numeric(as.character(join$CHROM))
join$POS<-as.numeric(as.character(join$POS))
join$JPT_AF<-as.numeric(as.character(join$JPT_AF))
join$JPT_MAF<-as.numeric(as.character(join$JPT_MAF))
join$NAG_AF<-as.numeric(as.character(join$NAG_AF))
join$NAG_MAF<-as.numeric(as.character(join$NAG_MAF))
join<-join[complete.cases(join),]

write.table(join, '~/Documents/QualityPaper/join.txt')
} #make the joint frequency file
join=fread('~/Documents/QualityPaper/Misc/join.txt')

##JPT Manhattan
if(f){
dir='/Volumes/gravel/luke_projects/1000Genomes/Regression/'
out='/Volumes/gravel/luke_projects/1000Genomes/MeanDev/'

jptNames = list.files(path=dir,pattern='_JPT_dev10.csv', full.names = T)
jReg = do.call(rbind, lapply(jptNames, function(x) fread(x, header=T, sep=' ')))
jReg<-jReg[which(jReg$Pos != 'Pos'),]
jReg$Chr<-as.numeric(as.character(jReg$Chr))
jReg$Pos<-as.numeric(as.character(jReg$Pos))
jReg$dev<-as.numeric(as.character(jReg$dev))

jReg$plog10 <- - pchisq(jReg$dev, 1, lower.tail=F, log.p=T)/log(10)
jReg$p <- pchisq(jReg$dev, 1, lower.tail=F)

sig.6<-jReg[which(jReg$plog10 >= 6),]
sig.4<-jReg[which(jReg$plog10 >= 4),]
NotSig<-jReg[which(jReg$plog10 < 6),]

write.table(sig.6[,c('Chr','Pos')],'~/Documents/Regression/JPT_sig6.Pos', quote=F, col.names=F, row.names=F)
write.table(sig.4[,c('Chr','Pos')],'~/Documents/Regression/JPT_sig4.Pos', quote=F, col.names=F, row.names=F)
}
#sig.6<-fread('~/Documents/Regression/JPT_sig6.Pos', col.names=c('CHR','POS'))
sig.6$sig<-'Significant'

#CONTEXT Mutation Spectrum enrichment
#test SNPs
SNP<-read.table('~/Documents/Regression/JPT_sig6.Context', sep='\t', header=F, col.names=c('Chr','Pos','Ref','Alt','Flip','Context'))
#Random SNPs
rSNP<-read.table('~/Documents/GWAS_Qual/Final_GWAS/1kGP_GenomeWide_JPT_INT_randomNOTsig6.Context.txt', sep='\t', header=F, col.names=c('Chr','Pos','Ref','Alt','Flip','Context'))
#function to get reverse complement
ReverseComp <- function(x)
        chartr("ATGC","TACG",
        sapply(lapply(strsplit(x, NULL), rev), paste, collapse=""))
#get reverse complement of context, ref and alt
SNP$Rev<-ReverseComp(as.character(SNP$Context))
SNP$RevRef<-ReverseComp(as.character(SNP$Ref))
SNP$RevAlt<-ReverseComp(as.character(SNP$Alt))
#only get one half (to fold over)
AC<-SNP[which(SNP$Ref %in% c('A','C')),]
TG<-SNP[which(SNP$Ref %in% c('T','G')),]
#fold over
TG$Context<-TG$Rev
TG$Ref<-TG$RevRef
TG$Alt<-TG$RevAlt
#bind them together
SNP<-rbind(AC,TG)
#combine it all together to get a context
SNP$Mut<-paste(SNP$Ref,'->', substr(SNP$Context,1,1),SNP$Alt,substr(SNP$Context,3,3),sep='')

#do the same thing as above but for the random SNPs
rSNP$Rev<-ReverseComp(as.character(rSNP$Context))
rSNP$RevRef<-ReverseComp(as.character(rSNP$Ref))
rSNP$RevAlt<-ReverseComp(as.character(rSNP$Alt))

AC<-rSNP[which(rSNP$Ref %in% c('A','C')),]
TG<-rSNP[which(rSNP$Ref %in% c('T','G')),]
TG$Context<-TG$Rev
TG$Ref<-TG$RevRef
TG$Alt<-TG$RevAlt

rSNP<-rbind(AC,TG)
rSNP$Mut<-paste(rSNP$Ref,'->', substr(rSNP$Context,1,1),rSNP$Alt,substr(rSNP$Context,3,3),sep='')

#get the counts of each bin
plt<-as.data.frame(table(SNP$Mut))
plt$Start<-substr(plt$Var1,4,4)
plt$End<-substr(plt$Var1,6,6)
plt$Mut<-paste(substr(plt$Var1,1,3), substr(plt$Var1,5,5))

rplt<-as.data.frame(table(rSNP$Mut))
rplt$Start<-substr(rplt$Var1,4,4)
rplt$End<-substr(rplt$Var1,6,6)
rplt$Mut<-paste(substr(rplt$Var1,1,3), substr(rplt$Var1,5,5))

#SFS
join<-merge(join, sig.6, by.x=c('CHROM', 'POS'), by.y=c('CHR','POS'), all.x=T)
join<-join[which((join$JPT_MAF!=0)&(join$NAG_MAF!=0)),] #remove fixed
jj<-join[complete.cases(join),]

SFS=ggplot(data= join, aes(x=NAG_MAF, y=JPT_MAF))+
geom_bin2d(binwidth=c(1/884, 1/104))+geom_point(data=jj, size=0.8, alpha=0.7, shape=3)+
scale_fill_distiller(palette = "Spectral",trans = "log", breaks=c(1,10,100,1000,10000,100000,1000000,10000000),labels=c('1','10', '100', '1,000', '10,000', '100,000', '1,000,000','10,000,000'))+
scale_x_continuous(expand=c(0.01,0.01))+
scale_y_continuous(expand=c(0.01,0.01))+
theme_classic()+theme(axis.line=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
xlab(label="Nagahama")

HIST=ggplot(jj, aes(x= JPT_MAF))+geom_histogram(binwidth=0.01, fill='grey')+theme_classic()+scale_x_continuous(limits=c(0,1),expand=c(0.01,0.01))+scale_y_reverse(expand=c(0.01,0.01), breaks=c(10,50,100))+geom_vline(xintercept=0.05, linetype=3)+coord_flip()+xlab('')+ylab('')+theme(axis.line=element_blank(),axis.text.x=element_text(angle=90))+xlab(label="1000 Genomes Project")
B=plot_grid(HIST,SFS,nrow=1,rel_widths=c(1,5), align = 'h', labels=c('','B'))


sig=ggplot(plt, aes(x=End, y=Start, fill=Freq))+
facet_grid(Mut~., switch="y")+geom_tile()+theme_classic()+
theme(axis.title=element_blank(),
axis.line=element_blank(),
legend.position="bottom",
strip.placement = "outside",
strip.background=element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1))+
scale_fill_gradient2(low='red', mid='white', high='blue',
guide = guide_legend(title = 'Count',title.position = "top"),
limits=c(0,max(plt$Freq)))+labs(subtitle='Suspicious SNPs')

notsig=ggplot(rplt, aes(x=End, y=Start, fill=Freq))+
facet_grid(Mut~., switch="y")+geom_tile()+theme_classic()+
theme(axis.title=element_blank(),
axis.line=element_blank(),
legend.position="bottom",
strip.placement = "outside",
strip.background=element_blank(),strip.text=element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1))+
scale_fill_gradient2(low='red', mid='white', high='blue',
guide = guide_legend(title = 'Count',title.position = "top"),
limits=c(0,max(plt$Freq)))+labs(subtitle='Control SNPs')

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(sig)

            
A=plot_grid(plot_grid(
sig + theme(legend.position="none"),
notsig + theme(legend.position="none"),
rel_widths=c(1,0.77)),
mylegend,
nrow=2,rel_heights=c(1,0.1),labels=c('A'))

AB=plot_grid(A,B, rel_widths=c(1,2.7))

Reg<-fread('/Volumes/gravel/luke_projects/1000Genomes/Regression/GenomeWide.Regression_JPT.csv')

Reg<-Reg[which(Reg$Pos != 'Pos'),]
Reg$Chr<-as.numeric(as.character(Reg$Chr))
Reg$Pos<-as.numeric(as.character(Reg$Pos))
Reg$dev<-as.numeric(as.character(Reg$dev))
Reg$Pop=NULL
Reg$p=NULL

Reg$plog10 <- - pchisq(Reg$dev, 1, lower.tail=F, log.p=T)/log(10)
Reg$p<-pchisq(Reg$dev, 1, lower.tail=F)
###########
colnames(Reg)[c(1,2)]<-c('CHROM','POS')

RegJoin<-left_join(Reg, join, by=c('CHROM','POS'))
#RegJoin<-RegJoin[which(RegJoin$JPT_MAF > 1/208),]

ggplot(RegJoin, aes(sample=RegJoin$dev))+stat_qq(distribution=stats::qchisq, dparams=list(df=1))+geom_qq_line(distribution=stats::qchisq, dparams=list(df=1))

ggsave('~/Documents/QualityPaper/Figures/QQplot_JPT.jpg', height=5, width=7)

ggplot(Reg, aes(sample=Reg$dev))+stat_qq(distribution=stats::qchisq, dparams=list(df=1))+geom_qq_line(distribution=stats::qchisq, dparams=list(df=1))
ggsave('~/Documents/QualityPaper/Figures/QQplot_JPT_withSingles.jpg', height=5, width=7)

single<-RegJoin[which(RegJoin$JPT_MAF == 0.00480769),]
ggplot(single, aes(sample= single$dev))+stat_qq(distribution=stats::qchisq, dparams=list(df=1))+geom_qq_line(distribution=stats::qchisq, dparams=list(df=1))
ggsave('~/Documents/QualityPaper/Figures/QQplot_JPT_justSingles.jpg', height=5, width=7)

double<-RegJoin[which(RegJoin$JPT_MAF == 2*0.00480769),]
ggplot(double, aes(sample= double$dev))+stat_qq(distribution=stats::qchisq, dparams=list(df=1))+geom_qq_line(distribution=stats::qchisq, dparams=list(df=1))
ggsave('~/Documents/QualityPaper/Figures/QQplot_JPT_justDouble.jpg', height=5, width=7)

triple<-RegJoin[which(RegJoin$JPT_MAF == 3*0.00480769),]
ggplot(triple, aes(sample= triple$dev))+stat_qq(distribution=stats::qchisq, dparams=list(df=1))+geom_qq_line(distribution=stats::qchisq, dparams=list(df=1))
ggsave('~/Documents/QualityPaper/Figures/QQplot_JPT_justTriple.jpg', height=5, width=7)
##########
sig.6<-Reg[which(Reg$plog10 > 6),]
NotSig<-Reg[which(Reg$plog10 < 6),]

GWAS=ggplot(NotSig, aes(x=Pos, y = plog10, color=as.factor(Chr)))+
facet_grid(~Chr, scales='free_x', space='free_x', switch='x')+
scale_fill_manual (values=getPalette(colourCount))+
scale_y_continuous(expand=c(0,0))+
geom_point(alpha=0.3, size=1)+theme_classic()+
geom_point(data=sig.6, aes(x=Pos, y = plog10,),color='black',shape=3)+
labs(y='-log10(p)', x='Chromosome')+
theme(plot.title = element_text(hjust = 0.5), 
plot.subtitle = element_text(hjust = 0.5), 
axis.text.x=element_blank(), 
axis.ticks.x=element_blank(),
axis.line.x=element_blank(), 
strip.background = element_blank(),
strip.text.x = element_text(size = 6))+
guides(color=F)
#geom_hline(yintercept = 8, color='red')+geom_hline(yintercept = 6, color='blue')+

fig1=plot_grid(AB, GWAS, ncol=1, rel_heights=c(2,1), labels=c('','C'))
ggsave('~/Documents/QualityPaper/Figure1.jpg',fig1, height=11,width=11)
