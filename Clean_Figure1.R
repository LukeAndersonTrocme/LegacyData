library(ggplot2)
library(colorspace)
library(cowplot)

join<-read.table('~/Documents/MutSpect/WrapUpPlots/SFS_NAG_JPT.plot.txt')

join$sum<-join$JPT_V6 + join$NAG_V6
j<-join[which((join$sum > 0.01) & (join$sum < 2)),]
j<-j[,c('JPT_V1', 'JPT_V2','JPT_V6', 'NAG_V6', 'sum')]

pltspect<-merge(j, sig6, by.x=c('JPT_V1', 'JPT_V2'), by.y=c('CHR','BP'), all.x=T)

jj<-pltspect[complete.cases(pltspect),]

SNP<-read.table('~/Documents/GWAS_Qual/GWAS_Data/CI/GenomeWide_JPT_Sig6.Context.txt', sep='\t', header=F, col.names=c('Chr','Pos','Ref','Alt','Flip','Context'))

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

plt<-as.data.frame(table(SNP$Mut))
plt$Start<-substr(plt$Var1,4,4)
plt$End<-substr(plt$Var1,6,6)
plt$Mut<-paste(substr(plt$Var1,1,3), substr(plt$Var1,5,5))

m<-nrow(SNP)/96
plt$Freq<-plt$Freq/m

SFS=ggplot(data= pltspect, aes(x=NAG_V6, y=JPT_V6))+
geom_bin2d(binwidth=c(1/100, 1/100))+geom_point(data=jj, size=0.8, alpha=0.7, shape=3)+
scale_fill_distiller(palette = "Spectral",trans = "log", breaks=c(1,10,100,1000,10000,100000,1000000),labels=c('1','10', '100', '1,000', '10,000', '100,000', '1,000,000'))+
scale_x_continuous(expand=c(0.01,0.01))+
scale_y_continuous(expand=c(0.01,0.01))+
theme_classic()+theme(axis.line=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
xlab(label="NAG")

HIST=ggplot(jj, aes(x= JPT_V6))+geom_histogram(binwidth=0.01, fill='grey')+theme_classic()+scale_x_continuous(limits=c(0,1),expand=c(0.01,0.01))+scale_y_reverse(expand=c(0.01,0.01), breaks=c(0,50,100))+geom_vline(xintercept=min(jj$JPT_V6), linetype=3)+coord_flip()+theme(axis.line=element_blank(), axis.title.x=element_blank(),axis.text.x=element_text(angle=90))+xlab(label="1000 Genomes Project")

SPECT=ggplot(plt, aes(x=End, y=Start, fill=Freq))+
facet_grid(Mut~., switch="y")+geom_tile()+theme_classic()+
theme(plot.title = element_blank(),
plot.subtitle = element_text(hjust = 0.5), 
axis.title=element_blank(),
legend.position="bottom",
strip.placement = "outside",
strip.background=element_blank())+
scale_fill_gradient(low='white', high='blue',
guide = guide_legend(title = 'Enrichment',title.position = "top"))

p3=plot_grid(HIST,SFS,SPECT,nrow=1,rel_widths=c(1,5,1.5), align = 'h', labels=c('A','','B'))
ggsave('~/Documents/QualityPaper/NAG_JPT_SFS_GenomeWide_frq.jpg',p3, height=5,width=7)

fname='/Users/luke/Documents/GWAS_Qual/GWAS_Data/CI/GenomeWide_JPT_GWAS_Qual_ci.assoc.linear'
tab5rows <- read.table(fname, header = TRUE, nrows = 5, row.names=NULL)
classes <- sapply(tab5rows, class)
assoc<-read.table(fname, header=T,colClasses= classes)
assoc$Plog10=-log10(assoc$P)


GWAS=ggplot(assoc, aes(x=BP, y = Plog10, color=as.factor(CHR)))+
facet_grid(~CHR, scales='free_x', space='free_x', switch='x')+
scale_fill_manual (values=getPalette(colourCount))+
scale_y_continuous(expand=c(0,0))+
geom_point(alpha=0.3, size=1)+theme_classic()+
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

fig1=plot_grid(p3, GWAS, ncol=1, rel_heights=c(3,2), labels=c('','C'))
ggsave('~/Documents/QualityPaper/Figure1.jpg',fig1, height=7,width=7)
