library(data.table)
library(ggplot2)
library(stringr)
library(cowplot)

rsID<-fread('~/Documents/QualityPaper/sig/Total_Sig_rsID.txt', header=F, col.names=c('rsID'))
pos<-fread('~/Documents/QualityPaper/sig/Total_Sig_POS.txt',col.names=c('CHROM','POS'))
tot<-cbind(rsID,pos)
tot $CHROM<-as.numeric(as.character(tot $CHROM))
tot $POS<-as.numeric(as.character(tot $POS))

GP<-fread('~/Documents/QualityPaper/sig/Total_Sig_FREQ.frq',fill=TRUE,col.names=c('CHROM','POS','count','AC','AF', 'AN'))
GP$AC1<-round(GP$AC*(1-GP$AF), digits=0)
GP$CHROM<-as.numeric(as.character(GP$CHROM))
GP$POS<-as.numeric(as.character(GP$POS))

GP<-merge(GP,tot, by=c('CHROM','POS'))


HM<-fread('/Users/luke/genomes/genomes/HapMap/hapmap_3.3.b37_sig_rsID.frq', sep='\t',fill=T, col.names=c('CHROM','POS','rsID','AF', 'AN'))
HM<-HM[which(complete.cases(HM)==T),]


plt1<-merge(GP, HM, by=c('rsID'))
plt1 $AF.x<-as.numeric(as.character(plt $AF.x))
plt1 $AF.y<-as.numeric(as.character(plt $AF.y))

A=ggplot(plt1, aes(x=1-AF.x, y=AF.y))+geom_point(shape=1)+geom_point(data=plt1[which(plt1$rsID=='rs6057648'),], color='blue')+geom_point(data=plt1[which(plt1$rsID=='rs301'),], color='green')+labs(x='Allele Frequency in HapMap', y='Allele Frequency in 1kGP', title='Significant SNPs found in HapMap and 1kGP')+theme_classic()+xlim(c(0,1))+ylim(c(0,1))

GN<-fread('/Users/luke/genomes/GnomAD/gnomad_AF.frq')
GN$AF_raw<-as.numeric(as.character(GN$AF_raw))
GN$AC_raw<-as.numeric(as.character(GN$AC_raw))
ggplot(GN, aes(x=reorder(POS,CHROM),y=AF_raw))+geom_point(shape=1)+theme_classic()+theme(axis.text.x=element_blank())

GNm<-melt(GN[,-c('AF_raw')], id=c('CHROM','POS','ID'))
GNm$value<-as.numeric(as.character(GNm$value))
GNm<-GNm[which(GNm$value!=0),]

plt<-merge(GP, GN, by=c('CHROM','POS'))

B=ggplot(plt, aes(y=1-AF, x=AF_raw))+geom_point(shape=1)+labs(x='Allele Frequency in GnomAD', y='Allele Frequency in 1kGP', title='Significant SNPs found in GnomAD and 1kGP')+theme_classic()
C=plot_grid(A,B)

ggsave('/Users/luke/Documents/QualityPaper/Figures/Hap_GnomAD.jpg',height=7,width=11)