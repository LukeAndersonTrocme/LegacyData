library(ggplot2)
library(colorspace)
library(cowplot)
f1kGP='~/Documents/GWAS_Qual/Final_GWAS/1kGP_GenomeWide_JPT.frq' #/Users/luke/genomes/genomes/1kGP_4bedFiltered/1kGP_GenomeWide.4bed_filtered.freq.frq
fnag='/Users/luke/genomes/genomes/NAG_4bedFiltered/NAG_GenomeWide.4bed_filtered.freq.frq'

tab5rows <- read.table(fnag, header = TRUE, nrows = 5, row.names=NULL)
classes <- sapply(tab5rows, class)
#Allele Frequency from JPT
JPT<-read.table(f1kGP, fill=T,skip=1,  skipNul = TRUE,colClasses = classes)
colnames(JPT)<-paste("JPT", colnames(JPT), sep = "_")
#Allele Frequency from NAG
NAG<-read.table(fnag, fill=T,skip=1,  skipNul = TRUE,colClasses = classes)
colnames(NAG)<-paste("NAG", colnames(NAG), sep = "_")

#merge them together
join<-merge(JPT[which(JPT$JPT_V6 <= 1),],
 			NAG[which(NAG$NAG_V6 <= 1),], 
 			by.x=c("JPT_V1","JPT_V2"), 
 			by.y=c("NAG_V1","NAG_V2"), all=T)
join[is.na(join)]<-0
#REMOVE  HWE SITES
join= join[! paste(join $JPT_V1, join $JPT_V2) %in% c(paste(hwe$V1,hwe$V2)), ]

join$sum<-join$JPT_V6 + join$NAG_V6
j<-join[which((join$sum > 0.01) & (join$sum < 2)),]

j<-j[,c('JPT_V1', 'JPT_V2','JPT_V6', 'NAG_V6', 'sum')]
write.table(j,'~/Documents/MutSpect/WrapUpPlots/SFS_NAG_JPT_4bed.plot.txt')
#j<-read.table('~/Documents/MutSpect/WrapUpPlots/SFS_NAG_JPT_4bed.plot.txt')

tab5rows <- read.table('/Users/luke/Documents/1FinalGWAS/1kGP_Chr2.4bed_filtered.assoc.linear', header = T, nrows = 5, row.names=NULL)
classes <- sapply(tab5rows, class)
GWAS_JPT<-read.table('~/Documents/GWAS_Qual/Final_GWAS/1kGP_GenomeWide_JPT_INT_sig6.assoc.linear', header=F,colClasses = classes, col.names=names(tab5rows))
GWAS_JPT$Plog10= -log10(GWAS_JPT$P)
sig6<-GWAS_JPT[which(GWAS_JPT$Plog10 > 6),]


pltspect<-merge(j, sig6, by.x=c('JPT_V1', 'JPT_V2'), by.y=c('CHR','BP'), all.x=T)

jj<-pltspect[complete.cases(pltspect),]

SNP<-read.table('/Users/luke/Documents/1FinalGWAS/1kGP_GenomeWide.4bed_filtered_sig6.Context.txt', sep='\t', header=F, col.names=c('Chr','Pos','Ref','Alt','Flip','Context'))

rSNP<-read.table('~/Documents/GWAS_Qual/Final_GWAS/1kGP_GenomeWide_JPT_INT_randomNOTsig6.Context.txt', sep='\t', header=F, col.names=c('Chr','Pos','Ref','Alt','Flip','Context'))

ReverseComp <- function(x)
        chartr("ATGC","TACG",
        sapply(lapply(strsplit(x, NULL), rev), paste, collapse=""))

SNP$Rev<-ReverseComp(as.character(SNP$Context))
SNP$RevRef<-ReverseComp(as.character(SNP$Ref))
SNP$RevAlt<-ReverseComp(as.character(SNP$Alt))

AC<-SNP[which(SNP$Ref %in% c('A','C')),]
TG<-SNP[which(SNP$Ref %in% c('T','G')),]
TG$Context<-TG$Rev
TG$Ref<-TG$RevRef
TG$Alt<-TG$RevAlt

SNP<-rbind(AC,TG)
SNP$Mut<-paste(SNP$Ref,'->', substr(SNP$Context,1,1),SNP$Alt,substr(SNP$Context,3,3),sep='')

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

plt<-as.data.frame(table(SNP$Mut))
plt$Start<-substr(plt$Var1,4,4)
plt$End<-substr(plt$Var1,6,6)
plt$Mut<-paste(substr(plt$Var1,1,3), substr(plt$Var1,5,5))

rplt<-as.data.frame(table(rSNP$Mut))
rplt$Start<-substr(rplt$Var1,4,4)
rplt$End<-substr(rplt$Var1,6,6)
rplt$Mut<-paste(substr(rplt$Var1,1,3), substr(rplt$Var1,5,5))

SFS=ggplot(data= pltspect, aes(x=NAG_V6, y=JPT_V6))+
geom_bin2d(binwidth=c(1/884, 1/104))+geom_point(data=jj, size=0.8, alpha=0.7, shape=3)+
scale_fill_distiller(palette = "Spectral",trans = "log", breaks=c(1,10,100,1000,10000,100000,1000000),labels=c('1','10', '100', '1,000', '10,000', '100,000', '1,000,000'))+
scale_x_continuous(expand=c(0.01,0.01))+
scale_y_continuous(expand=c(0.01,0.01))+
theme_classic()+theme(axis.line=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
xlab(label="NAG")

HIST=ggplot(jj, aes(x= JPT_V6))+geom_histogram(binwidth=0.01, fill='grey')+theme_classic()+scale_x_continuous(limits=c(0,1),expand=c(0.01,0.01))+scale_y_reverse(expand=c(0.01,0.01), breaks=c(0,50,100))+geom_vline(xintercept=min(jj$JPT_V6), linetype=3)+coord_flip()+xlab('')+ylab('')+theme(axis.line=element_blank(),axis.text.x=element_text(angle=90))+xlab(label="1000 Genomes Project")
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

AB=plot_grid(A,B, rel_widths=c(1,2.5))
#ggsave('~/Documents/QualityPaper/NAG_JPT_SFS_GenomeWide_frq.jpg',p3, height=5,width=7)

fname='~/Documents/GWAS_Qual/Final_GWAS/1kGP_GenomeWide_JPT_noNA_Dup.assoc.linear'
tab5rows <- read.table(fname, header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
assoc<-read.table(fname, header=TRUE, colClasses= classes)
assoc$Plog10=-log10(assoc$P)
sig.6<-assoc[which(assoc$Plog10 > 6),]
NotSig<-assoc[which(assoc$Plog10 < 6),]

GWAS=ggplot(NotSig, aes(x=BP, y = Plog10, color=as.factor(CHR)))+
facet_grid(~CHR, scales='free_x', space='free_x', switch='x')+
scale_fill_manual (values=getPalette(colourCount))+
scale_y_continuous(expand=c(0,0))+
geom_point(alpha=0.3, size=1)+theme_classic()+
geom_point(data=sig.6, aes(x=BP, y = Plog10,),color='black',shape=3)+
geom_hline(yintercept = 8, color='red')+
geom_hline(yintercept = 6, color='blue')+
labs(y='-log10(p)', x='Chromosome')+
theme(plot.title = element_text(hjust = 0.5), 
plot.subtitle = element_text(hjust = 0.5), 
axis.text.x=element_blank(), 
axis.ticks.x=element_blank(),
axis.line.x=element_blank(), 
strip.background = element_blank(),
strip.text.x = element_text(size = 6))+
guides(color=F)

fig1=plot_grid(AB, GWAS, ncol=1, rel_heights=c(2,1), labels=c('','C'))
ggsave('~/Documents/QualityPaper/Figure1.jpg',fig1, height=10,width=10)
