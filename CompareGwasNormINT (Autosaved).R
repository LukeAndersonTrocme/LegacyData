##Compare GWAS p-values
library(ggplot2)
library(reshape2)

colors<-read.table("~/Dropbox/LukeTemp/NamePopBigPop_color.html.txt", comment.char='>',col.names=c('Pop','hex','BigPop'))
MyColour <- as.character(colors $hex)
names(MyColour) <- colors $Pop

pub<-read.table('~/Desktop/22SNPs_catalogue.txt', sep='\t')

pops<-list.files(path="/Users/luke/Documents/SigSnps",pattern=".frq$", full.names=T)

Freq<-data.frame()
for(f in seq(1,length(pops))){
pop<-unique(read.table(pops[f], row.names=NULL,header=T))
names(pop)=c('Chr', 'POS', 'N_ALLELES','N_CHR', 'Y.ALLELE.FREQ','X.ALLELE.FREQ')
pop$Pop=substr(pops[f],31,33)
Freq<-rbind(Freq, pop)
}

Freq$X.ALLELE.FREQ<-gsub('T:','',Freq$X.ALLELE.FREQ)
Freq$X.ALLELE.FREQ<-gsub('C:','',Freq$X.ALLELE.FREQ)
Freq$X.ALLELE.FREQ<-gsub('G:','',Freq$X.ALLELE.FREQ)
Freq$X.ALLELE.FREQ<-gsub('A:','',Freq$X.ALLELE.FREQ)

perPop<-read.table('/Users/luke/Documents/SigSnps/PerPop.txt', col.names=c('Chr','POS','Pop'))
perPop$Sig<-'Yes'


Norm<-read.table('~/Documents/GWAS_Qual/PublishedSNPs_GWAS.txt')

Int<-read.table('~/Documents/GWAS_Qual/PublishedSNPs_INT_GWAS.txt')

Norm$V1<-gsub('/Users/luke/Documents/GWAS_Qual/Final_GWAS/1kGP_GenomeWide_', '',Norm$V1)
Int$V1<-gsub('/Users/luke/Documents/GWAS_Qual/Final_GWAS/1kGP_GenomeWide_', '',Int$V1)

Norm$V1<-gsub('.assoc.linear:', '',Norm$V1)
Int$V1<-gsub('_INT.assoc.linear:', '',Int$V1)
Int<-Int[c('V1','V2', 'V3','V4','V10')]
Norm <-Norm[c('V1','V2', 'V3','V4','V10')]
Int$V10<--log10(Int$V10)
Norm$V10<--log10(Norm$V10)

both<-merge(Norm, Int, by=c('V1','V2', 'V3','V4'))
names(both)=c('Pop','Chr','SNP','POS','Pval_Norm','Pval_INT')
both<-merge(both, pub, by.x='SNP', by.y='V5')
both<-both[,c('Pop','Chr','SNP','POS','Pval_Norm','Pval_INT', 'V3')]

moth<-melt(both,id=c('Pop','Chr','SNP','POS', 'V3'))

ordr<-c('FIN', 'GBR', 'CEU', 'IBS', 'TSI', 'CHS', 'CDX', 'CHB', 'JPT', 'KHV', 'GIH', 'STU', 'PJL', 'ITU', 'BEB', 'PEL', 'MXL', 'CLM', 'PUR', 'ASW', 'ACB', 'GWD', 'YRI', 'LWK', 'ESN', 'MSL')

moth$ordered <- factor(moth$Pop, levels =ordr)

moth<-merge(moth,perPop, by=c('Chr', 'POS','Pop'),all=T)
moth<-merge(moth,Freq, by=c('Chr', 'POS','Pop'),all=T)



ggplot(moth, aes(x= ordered, y=value, color= ordered,group=1,shape=variable))+geom_hline(yintercept=0, color='grey70', linetype=1)+geom_hline(yintercept=6, color='grey70', linetype=3)+geom_point()+facet_wrap(.~V3,ncol=1)+theme_classic()+theme(axis.title=element_blank(), axis.ticks.y=element_blank(),axis.line=element_blank(),axis.text.y=element_text(size=5),axis.text.x=element_text(angle=90, size=5, vjust=-0.05))+scale_color_manual(values=MyColour)+scale_y_continuous(breaks=c(0,6))+scale_shape_manual(labels=c('Raw','INT'), values=c(0,3), name='pvalue')+guides(color=F)

ggsave('~/Documents/QualityPaper/CompareGWAS.jpg', height=10, width=4)

ggplot(moth, aes(x= ordered, y=as.numeric(X.ALLELE.FREQ), color= ordered,group=1))+geom_hline(yintercept=0, color='grey70', linetype=1)+geom_hline(yintercept=0.5, color='grey70', linetype=3)+geom_point(shape=2)+facet_wrap(.~V3,ncol=1)+theme_classic()+theme(axis.title=element_blank(), axis.ticks.y=element_blank(), axis.line=element_blank(), axis.text.y=element_text(size=5), axis.text.x=element_text(angle=90, size=5, vjust=-0.05))+scale_color_manual(values=MyColour)+scale_shape_manual(labels=c('Raw','INT'), values=c(0,3), name='pvalue')+guides(color=F)+geom_point(data=moth[complete.cases(moth),], shape=1, size=4, color='black')+labs(x='Population',y='Allele Frequency')+scale_y_continuous(breaks=c(0,0.5))

ggsave('~/Documents/QualityPaper/AFperPublication.jpg', height=10, width=4)
