##Hardy GWAS qq plot
library(ggplot2)
Pops<-read.table('~/genomes/genomes/PopNames/NamePop.txt')
for(f in seq(1,26,1)){
pop=Pops[f,]
print(pop)

fname=paste('~/Documents/GWAS_Qual/Final_GWAS/1kGP_GenomeWide_hwe_GWAS_',pop,'_ND.txt',sep='')
t5<-read.table(fname,header=T,nrows=5)
classes <- sapply(t5, class)
t <- read.table(fname, header = TRUE, colClasses = classes, col.names=c('CHR','BP','HWE','GWAS'))
t<-t[which((t$HWE != 1) & (t$GWAS != 1)),]

t$P.HWE<--log10(t$HWE)
t$P.GWAS<--log10(t$GWAS)

ggplot(t, aes(x=P.HWE, y=P.GWAS))+geom_bin2d(binwidth=0.15)+theme_classic()+labs(x='Hardy-Weinberg -log10(p)',y='GWAS -log10(p)',title=pop)+scale_fill_gradient(trans='log',breaks=c(1,10,100,1000,10000,100000,1000000,10000000))
ggsave(paste('~/Documents/GWAS_Qual/HWEvsGWAS_',pop,'.jpg',sep=''))
}
