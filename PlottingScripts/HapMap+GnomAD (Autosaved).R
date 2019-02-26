library(data.table)
library(ggplot2)
library(stringr)
library(cowplot)

#GP<-fread('~/Documents/QualityPaper/sig/Total_Sig_POSfreq.txt',fill=TRUE,col.names=c('CHROM','POS','rsID', 'AF'))
GP<-fread('/Users/luke/Documents/QualityPaper/sig/SigVar_0.001.frq',fill=TRUE,col.names=c('CHROM','POS','rsID', 'AF'))
GP$CHROM<-as.numeric(as.character(GP$CHROM))
GP$POS<-as.numeric(as.character(GP$POS))



HM<-fread('/Users/luke/genomes/genomes/HapMap/hapmap_3.3.b37_sig_rsID.frq', sep='\t',fill=T, col.names=c('CHROM','POS','rsID','AF', 'AN'))
HM<-HM[which(complete.cases(HM)==T),]


plt1<-merge(GP, HM, by=c('CHROM','POS','rsID'))

ggplot(plt1, aes(x=AF.y, y=AF.x))+geom_point(shape=1)+geom_point(data=plt1[which(plt1$rsID=='rs6057648'),], color='blue')+geom_point(data=plt1[which(plt1$rsID=='rs301'),], color='green')+labs(x='Allele Frequency in HapMap', y='Allele Frequency in 1kGP', title='Significant SNPs found in HapMap and 1kGP')+theme_classic()+xlim(c(0,1))+ylim(c(0,1))

ggsave('/Users/luke/Documents/QualityPaper/Figures/Hap.jpg',height=7,width=7)


GN<-fread('/Users/luke/genomes/GnomAD/gnomad_AF.frq')
GN$AF_raw<-as.numeric(as.character(GN$AF_raw))
GN$AC_raw<-as.numeric(as.character(GN$AC_raw))
ggplot(GN, aes(x=reorder(POS,CHROM),y=AF_raw))+geom_point(shape=1)+theme_classic()+theme(axis.text.x=element_blank())

GNm<-melt(GN[,-c('AF_raw')], id=c('CHROM','POS','ID'))
GNm$value<-as.numeric(as.character(GNm$value))
GNm<-GNm[which(GNm$value!=0),]

plt<-merge(GP, GN, by=c('CHROM','POS'))

B=ggplot(plt, aes(y=AF, x=AF_raw))+geom_point(shape=1)+labs(x='Allele Frequency in GnomAD', y='Allele Frequency in 1kGP', title='Significant SNPs found in GnomAD and 1kGP')+theme_classic()
C=plot_grid(A,B)



nrow(merge(plt, HM, by=c('CHROM','POS','rsID')))