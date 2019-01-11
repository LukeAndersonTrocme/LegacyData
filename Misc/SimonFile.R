#SIMON FILE
sb<-fread('~/Dropbox/LukeTemp/SubMeta.txt')
sb$V1<-NULL



PCs<-list.files(path="~/Documents/PCAperPop/",pattern=".eigenvec",full.names=T)
all_PCs<-data.frame()
for(p in PCs){
	pc<-fread(p)
	
	all_PCs<-rbind(all_PCs,pc[,c(1:4)])
}

s=merge(sb,all_PCs,by.x='Name',by.y='V1')
s$V2<-NULL


PC <-fread('/Users/luke/Documents/PCAperPop/ALL_Pops.eigenvec')
s=merge(sb,PC[,c('V1','V3','V4','V5','V6','V7')],by.x='Name',by.y='V1')
colnames(s)=c('Name','Pop','average_quality_of_mapped_bases','PC1','PC2','PC3','PC4','PC5')

write.table(s, file='~/Documents/PCAperPop/Name_Pop_Qual_PC1_PC2_PC3_PC4_PC5.txt', quote=F, row.names=F)