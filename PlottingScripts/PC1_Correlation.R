library(ggplot2)
library(gridExtra)
library(reshape2)
library(cowplot)
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
  
colnames(tot.count)<-NamePop$BigPop #adds the POP name to the column

tot.count <- prop.table(as.matrix(tot.count),2)
tot.count <- tot.count[ , apply(tot.count, 2, var) != 0]

pop.count<-tot.count[, grep('EAS', colnames(tot.count))]

EAS <- c('JPT','CHS','CHB','CDX','KHV')
colnames(pop.count)<-NamePop[which(NamePop$BigPop == 'EAS'),]$Name
pop.count<-pop.count[, -grep('NAG', colnames(pop.count))]


pop.pca = prcomp(t(pop.count), scale=T)
pop.pca.scores = as.data.frame(pop.pca$x)[,1:5]
pop.pca.scores$Name<-rownames(pop.pca.scores)

plot<-merge(pop.pca.scores, Meta, by='Name')
plot$Out<-gsub('OutLiers', 'Outliers',plot$Out)

p1 <- ggplot(plot, aes(x= average_quality_of_mapped_bases, y=PC1, color=Pop, shape=Out))+
geom_point()+theme_classic()+labs(y='PC1',x='Average quality of mapped bases')+scale_colour_manual(values= MyColour)+guides(color=F, shape=F)#+geom_smooth(method = "lm", se = FALSE,size=0.5)

p2<-ggplot(plot, aes(x= SUBMISSION.DATE, y=plot$PC1, color=Pop, shape=Out))+
geom_point()+theme_classic()+labs(y='PC1',x='Submission date')+scale_colour_manual(values= MyColour)+guides(color=F, shape=F)

                     
p3<-ggplot(plot, aes(x= mean_insert_size, y=plot$PC1, color=Pop, shape=Out))+
geom_point()+theme_classic()+labs(y='PC1',x='Mean insert size')+scale_colour_manual(values= MyColour)+guides(color=F, shape=F)

p4<-ggplot(plot, aes(x=X._of_mismatched_bases, y=plot$PC1, color=Pop, shape=Out))+
geom_point()+theme_classic()+labs(x='Number of mismatched bases',y='PC1')+scale_x_reverse()+scale_colour_manual(values= MyColour)+guides(shape=guide_legend(title=""),color=guide_legend(title="Population"))

pf<-plot_grid(p1,p2,p3,p4, labels=c('A','B','C','D'), nrow=1,rel_widths=c(1,1,1,1.2))

ggsave('~/Documents/QualityPaper/Figures/PC1_Correlation.jpg',pf,height=6,width=14)
ggsave('~/Documents/QualityPaper/Figures/PC1_Correlation.tiff',pf,height=6,width=14)

pt<-plot[,c('Name','PC1','X._total_reads','X._mapped_reads_properly_paired','mean_insert_size','insert_size_median_absolute_deviation','STUDY_ID','ANALYSIS_GROUP','Transitions','Hom_ref','X._total_bases','X._mapped_reads','X._of_mismatched_bases','X._mapped_reads_paired_in_sequencing','median_insert_size','X._duplicate_bases')]
pt<-unique(melt(pt, id=c('Name','PC1')))

ggplot(pt, aes(x=PC1, y=value))+facet_wrap(.~variable)+
geom_point()+geom_smooth(method='lm',formula=y~x)+
theme_classic()+labs(title=paste("PC1 vs Quality for JPT"), subtitle=t1, x='PC1',y='')+
theme(plot.title = element_text(hjust= 0.5),plot.subtitle = element_text(hjust= 0.5))

ggsave('~/Documents/QualityPaper/PC1_Correlation.jpg',height=4,width=8)