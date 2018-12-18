Reg<- fread('/Users/luke/Documents/QualityPaper/Misc/AllRegressionsOver10.txt')
Reg<-Reg[which(Reg$Pos != 'Pos'),]
Reg$Chr<-as.numeric(as.character(Reg$Chr))
Reg$Pos<-as.numeric(as.character(Reg$Pos))
Reg$dev<-as.numeric(as.character(Reg$dev))
Reg$plog10 <- - pchisq(Reg$dev, 1, lower.tail=F, log.p=T)/log(10)
Reg<-Reg[which(Reg$plog10 >= 6),]
RegCH<-Reg[which(Reg$Pop == 'CHS' | Reg$Pop == 'CHB'),]

write.table(RegCH[,c('Chr','Pos')],'/Users/luke/Documents/QualityPaper/sig/CHINESE_pos.txt', quote=F, col.names=F, row.names=F)

Overlap<- unique(RegCH[which(RegCH$Pos %in% awkPos$V2),c('Chr','Pos')])

KGP2 <- unique(merge(KGP, Overlap, by.x=c('V1','V2'), by.y=c('Chr','Pos')))
KGP2 <- unique(merge(KGP1, Reg, by.x=c('V1','V2'), by.y=c('Chr','Pos')))



KGP1 <- read.table('/Users/luke/Documents/QualityPaper/SigVCF/1kGP_GenomeWide_83Han.frq',
					sep="\t",fill=TRUE,col.names=paste("V",1:8,sep=''), skip=1)