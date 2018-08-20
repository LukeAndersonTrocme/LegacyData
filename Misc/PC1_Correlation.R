library(ggplot2)
library(gridExtra)
#read colors per pop (used in 1kG Nature)
colors<-read.table("~/Dropbox/LukeTemp/NamePopBigPop_color.html.txt", comment.char='>',col.names=c('Pop','hex','BigPop'))
MyColour <- as.character(colors $hex)
names(MyColour) <- colors $Pop


#PC1 for each Pop by date
Meta<-read.table("~/genomes/genomes/1000GenomesMetaData.txt")
Meta$SUBMISSION.DATE<-as.Date(Meta$SUBMISSION.DATE, "%Y-%m-%d")
##Name, BigPop, Pop
NamePop<-read.table("~/bin/smaller_mut_spectrum_pipeline/1000genomes_phase3_sample_IDs_BigPop_NAG.txt")
NamePop<-NamePop[,c(1,2,3)]
names(NamePop)<-c('Name', 'BigPop', 'Pop')
##Name of OutLiers
CutOffIDs<-read.table("/Users/luke/Documents/MutSpect/OutputFiles/PCA_cutoff/CutOffID.txt")
CutOffIDs$Out<-'OutLiers'
NamePop <-merge(NamePop,CutOffIDs,by.x="Name", by.y="V1", all=T)
NamePop[is.na(NamePop)] <- 'Normal'

##MutationSpectrum Data
count<-read.table("~/Documents/MutSpect/FilteredData/2017-12-02/10/10-10-weekender/2017-12-02.chr10.MutSpect/files/_derived_each_lineage_chr10_nosingle.txt", stringsAsFactors=F, fill=T)

count[1,1] = c("AAA_C")
colnames(count)[1]<- "Mut_type"

#get the frequency of each mutation. 
tot.count <- sapply(seq(2,ncol(count),2), function(i) {
  rowMeans(count[,c(i, i+1)], na.rm=T)})
  
colnames(tot.count)<-NamePop$Pop #adds the POP name to the column

tot.count <- prop.table(as.matrix(tot.count),2)
tot.count <- tot.count[ , apply(tot.count, 2, var) != 0]

pop.count<-tot.count[, grep('JPT', colnames(tot.count))]
colnames(pop.count)<-NamePop[which(NamePop$Pop == 'JPT'),]$Name
pop.pca = prcomp(t(pop.count), scale=T)
pop.pca.scores = as.data.frame(pop.pca$x)[,1:5]
pop.pca.scores$Name<-rownames(pop.pca.scores)

plot<-merge(pop.pca.scores, Meta[,c('Name','Pop','BigPop',
									'SUBMISSION.DATE', 
									'average_quality_of_mapped_bases',
									'CENTER_NAME')], by='Name')
									
									
r1=lm(plot$average_quality_of_mapped_bases ~ plot$PC1, plot)
t1=paste("Adj R2 = ",signif(summary(r1)$adj.r.squared, 3),
                     " P =",signif(summary(r1)$coef[2,4], 3))

ggplot(plot, aes(y= average_quality_of_mapped_bases, x=plot$PC1))+
geom_point()+geom_smooth(method = "lm", se = FALSE,size=0.5)+
theme_classic()+labs(title=paste("PC1 vs Quality for JPT"), subtitle=t1, x='PC1',y='Average quality of mapped bases')+
theme(plot.title = element_text(hjust= 0.5),plot.subtitle = element_text(hjust= 0.5))

ggsave('~/Documents/QualityPaper/PC1_Correlation.jpg',height=4,width=8)