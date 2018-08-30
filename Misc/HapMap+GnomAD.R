library(data.table)
library(ggplot2)
library(stringr)

GP<-fread('/Users/luke/Documents/QualityPaper/sig/GenomeWide_sig0.01_inHapMap.frq', fill=T)
HM<-fread('/Users/luke/genomes/genomes/HapMap/hapmap_3.3.b37_Significant0.1SNPs.frq', sep='\t')

HM$AF<-str_split_fixed(HM$V4, ";", 3)[,2]
HM$AF<-as.numeric(as.character(str_split_fixed(HM$AF, "=", 2)[,2]))

plot<-merge(GP, HM, by.x=c('CHROM','POS'), by.y=c('V1','V2'))

ggplot(plot, aes(x=V6, y=AF))+geom_point()+geom_point(data=plot[which(plot$V3=='rs6057648'),], color='blue')+labs(x='Allele Frequency in HapMap', y='Allele Frequency in 1kGP', title='Significant SNPs found in HapMap and 1kGP')+theme_classic()

GN<-fread('/Users/luke/genomes/GnomAD/gnomad_AF.frq', colClasses=c('numeric','numeric','character','numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric'))

ggplot(GN, aes(x=reorder(POS,CHROM),y=AF_raw))+geom_point(shape=1)+theme_classic()+theme(axis.text.x=element_blank())

GNm<-melt(GN[,-c('AF_raw')], id=c('CHROM','POS','ID'))
GNm$value<-as.numeric(as.character(GNm$value))
GNm<-GNm[which(GNm$value!=0),]
ggplot(GNm, aes(x=reorder(POS,CHROM),y=value, color=variable))+geom_point(alpha=0.4)+theme_classic()+theme(axis.text.x=element_blank())+labs(y='Frequency',x='Position')