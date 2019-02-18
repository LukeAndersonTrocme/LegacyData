##Figure 1 Rejected

#f1kGP='~/genomes/genomes/1kGP_4bedFiltered/1kGP_GenomeWide_Frequency.txt'#'/Users/luke/genomes/genomes/1kGP_4bedFiltered/1kGP_GenomeWide_JPT.frq'
#'/Users/luke/genomes/genomes/NAG_4bedFiltered/NAG_GenomeWide.4bed_filtered.freq.frq'#~/genomes/genomes/NAG_4bedFiltered/NAG_GenomeWide.4bed_filtered.frq'

#col.names=c('rsID','CHROM', 'A1','A2','JPT_AF','NCHR','Chr','POS'))
#Allele Frequency from NAG
NAG<-fread(fnag, fill=T,col.names=c('CHROM','POS','N_Allele', 'N_Chr','NAG_AF','NAG_MAF'))
NAG $NAG_AF<-as.numeric(as.character(NAG $NAG_AF))
NAG $NAG_MAF<-as.numeric(as.character(NAG $NAG_MAF))
NAG $CHROM<-as.numeric(as.character(NAG $CHROM))
NAG $POS<-as.numeric(as.character(NAG $POS))
NAG<-NAG[complete.cases(NAG),]


if(F){
##PYTHON
JPT=pd.read_table('/Users/luke/genomes/genomes/1kGP_4bedFiltered/1kGP_GenomeWide_JPT_less.frq',
                  skiprows=1, names=['CHROM','POS','AF_JPT'])
JPT['AF_JPT']=JPT['AF_JPT'].str[2:]
NAG=pd.read_table('/Users/luke/genomes/genomes/NAG_4bedFiltered/NAG_GenomeWide.4bed_filtered_less.frq',
                  skiprows=1, names=['CHROM','POS','AF_NAG'])
NAG['AF_NAG']=NAG['AF_NAG'].str[2:]
join=pd.concat([JPT, NAG])
join.to_csv('~/Documents/QualityPaper/joinedAF.csv')
}


join<-fread('~/Documents/QualityPaper/joinedAF.csv', colClasses=c('numeric','numeric','numeric','numeric','numeric'))
join$V1<-NULL
join$AF_JPT<-as.numeric(as.character(join$AF_JPT))
join$AF_NAG<-as.numeric(as.character(join$AF_NAG))
join$CHROM<-as.numeric(as.character(join$CHROM))
join$POS<-as.numeric(as.character(join$POS))
join<-join[complete.cases(join),]



join$sum<-join$JPT_AF + join$NAG_AF
j<-join[which((join$sum > 0.01) & (join$sum < 2)),]

j<-j[,c('JPT_CHROM', 'JPT_POS','JPT_MAF', 'NAG_MAF', 'sum')]
write.table(j,'~/Documents/MutSpect/WrapUpPlots/SFS_NAG_JPT_4bed.plot.txt')
#j<-read.table('~/Documents/MutSpect/WrapUpPlots/SFS_NAG_JPT_4bed.plot.txt')



##Figure 3
SFS=data.frame()
for(f in seq(1,nrow(pops))){
pop=pops[f,]
print(as.character(pop))
fileNames = list.files(path=dir, pattern=paste(pop), full.names = T)
Reg = do.call(rbind, lapply(fileNames, function(x) fread(x)))
print('Read')
Reg<-Reg[which(Reg$Pos != 'Pos'),]
Reg$Chr<-as.numeric(as.character(Reg$Chr))
Reg$Pos<-as.numeric(as.character(Reg$Pos))
Reg$dev<-as.numeric(as.character(Reg$dev))
Reg$p=NULL
print('Clean')
Reg$plog10 <- - pchisq(Reg$dev, 1, lower.tail=F, log.p=T)/log(10)
sig.6<-Reg[which(Reg$plog10 > 6),]
rm(Reg)
write.table(sig.6[,c('Chr','Pos','plog10', 'Pop')], paste('~/Documents/Regression/',pop,'_sig6.csv', sep=','))
#SFS=rbind(SFS,sig.6[,c('Chr','Pos','plog10', 'Pop')])
}


