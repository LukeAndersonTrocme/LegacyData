library(data.table)
library(ggplot2)
library(stringr)
library(cowplot)
GP<-fread('/Users/luke/Documents/QualityPaper/sig/GenomeWide_sig0.01.frq', skip=1,fill=T,col.names=c('CHROM','POS','count','AC','AF', 'AN'))
GP$AC1<-round(GP$AC*(1-GP$AF), digits=0)

HM<-fread('/Users/luke/genomes/genomes/HapMap/hapmap_3.3.b37_Significant0.1SNPs_format.frq', sep='\t',fill=T, col.names=c('CHROM','POS','rsID','AC','AF', 'AN'))
HM<-HM[which(complete.cases(HM)==T),]


plot<-merge(GP, HM, by=c('CHROM','POS'))

A=ggplot(plot, aes(x=1-AF.x, y=AF.y))+geom_point(shape=1)+geom_point(data=plot[which(plot$rsID=='rs6057648'),], color='blue')+labs(x='Allele Frequency in HapMap', y='Allele Frequency in 1kGP', title='Significant SNPs found in HapMap and 1kGP')+theme_classic()

GN<-fread('/Users/luke/genomes/GnomAD/gnomad_AF.frq')
GN$AF_raw<-as.numeric(as.character(GN$AF_raw))
GN$AC_raw<-as.numeric(as.character(GN$AC_raw))
ggplot(GN, aes(x=reorder(POS,CHROM),y=AF_raw))+geom_point(shape=1)+theme_classic()+theme(axis.text.x=element_blank())

GNm<-melt(GN[,-c('AF_raw')], id=c('CHROM','POS','ID'))
GNm$value<-as.numeric(as.character(GNm$value))
GNm<-GNm[which(GNm$value!=0),]
ggplot(GNm, aes(x=reorder(POS,CHROM),y=value, color=variable))+geom_point(alpha=0.4)+theme_classic()+theme(axis.text.x=element_blank())+labs(y='Frequency',x='Position')

plt<-merge(GP, GN, by=c('CHROM','POS'))

B=ggplot(plt, aes(y=1-AF, x=AF_raw))+geom_point(shape=1)+labs(x='Allele Frequency in GnomAD', y='Allele Frequency in 1kGP', title='Significant SNPs found in GnomAD and 1kGP')+theme_classic()
C=plot_grid(A,B)

ggsave('/Users/luke/Documents/QualityPaper/Figures/Hap_GnomAD.jpg',height=7,width=11)