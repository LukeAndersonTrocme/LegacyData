##rank-based inverse normal transformation

pops<-read.table('~/genomes/genomes/PopNames/Name.BigPop.Pop.txt')
pops<-pops[which(pops$V3 != 'NAG'),]
pheno<-read.table('~/genomes/genomes/1000GenomesQual.pheno')

final.pheno=data.frame()
for(p in unique(pops$V3)){
#subset table
print(p)
sub=pheno[which(pheno$V1 %in% pops[which(pops$V3 == p),]$V1),]
#INT transformation
#https://www.biostars.org/p/80597/
sub$INT=qnorm((rank(sub$V3,na.last="keep")-0.5)/sum(!is.na(sub$V3)))
final.pheno=rbind(final.pheno, sub)
}

write.table(final.pheno[,c('V1','V2','INT')], file='~/genomes/genomes/1000GenomesQual_INT.pheno',row.names=F, col.names=F, quote=F)