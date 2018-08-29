library(ggplot2)
library(ggthemes)
library(cowplot)

Pops<-read.table("~/Dropbox/LukeTemp/NamePop.txt")
Meta<-read.table("~/Dropbox/LukeTemp/1000GenomesMetaData.txt")
Meta$SUBMISSION.DATE<-as.Date(Meta$SUBMISSION.DATE, "%Y-%m-%d")
phase<-read.table('~/genomes/genomes/PHASE1_igsr_samples.tsv', header=T, sep='\t')
phase$Phase=1
Meta<-merge(Meta, phase[c('Sample.name','Phase')],by.x='Name',by.y='Sample.name', all=T)
Meta[which(is.na(Meta$Phase)),]$Phase=3
colors<-read.table("~/Dropbox/LukeTemp/NamePopBigPop_color.html.txt", comment.char='>',col.names=c('Pop','hex','BigPop'))
MyColour <- as.character(colors $hex)
names(MyColour) <- colors $Pop

pop.batch<-unique(
Meta[,c('Name','SUBMISSION.DATE','average_quality_of_mapped_bases','Pop','BigPop','CENTER_NAME','Phase','INSTRUMENT_MODEL')])
pop.batch$CENTER_NAME<-gsub('Illumina','ILLUMINA',pop.batch$CENTER_NAME)
pop.batch$INSTRUMENT_MODEL <-gsub('Illumina','',pop.batch$INSTRUMENT_MODEL)
pop.batch$INSTRUMENT_MODEL <-gsub('Genome Analyzer','GA',pop.batch$INSTRUMENT_MODEL)
pop.batch$INSTRUMENT_MODEL <-gsub('HiSeq 2000','HS 2k',pop.batch$INSTRUMENT_MODEL)

pop.batch<-pop.batch[complete.cases(pop.batch),]
pop.batch<-pop.batch[!duplicated(pop.batch$Name),]

Var<-aggregate(pop.batch$average_quality_of_mapped_bases,list(pop.batch$Pop),var)
Var<-merge(Var, pop.batch, by.x='Group.1', by.y='Pop')
plt <- Var[with(Var, order(x, average_quality_of_mapped_bases)),]

Day<-aggregate(plt$SUBMISSION.DATE, list(plt$Group.1), mean)
Day<-merge(Day, pop.batch, by.x='Group.1', by.y='Pop')
plt <- Day[with(Day, order(x,Group.1, average_quality_of_mapped_bases)),]

plt$Name <- factor(plt$Name, levels = unique(plt$Name))

a=ggplot(plt, aes(x=Name, y= average_quality_of_mapped_bases, color=Group.1, shape=as.factor(Phase)))+geom_point()+
scale_colour_manual(values= MyColour, breaks=unique(plt$Group.1))+
scale_x_discrete(expand=c(0.1,0))+scale_shape_manual(values=c(1,3))+
guides(color=guide_legend(title="Population",
override.aes = list(shape = 15, size=3),ncol=1),shape=F)+
labs(title="Data quality over time",x="Populations ranked by sequencing date", y="Average quality of mapped bases")+
theme(plot.title = element_text(hjust= 0.5), 
axis.text.x=element_blank(), axis.ticks.x=element_blank())+geom_hline(yintercept=30, linetype=2, color='grey')

b=ggplot(pop.batch, aes(
x=SUBMISSION.DATE, y= average_quality_of_mapped_bases, color=CENTER_NAME, shape=as.factor(Phase)))+geom_point()+scale_shape_manual(values=c(1,3),name='Phase')+
guides(color=guide_legend(title="Sequencing\n    Center",
override.aes = list(shape = 15, size=5)))+
labs(y="Average quality of mapped bases")+geom_hline(yintercept=30, linetype=2, color='grey')+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.line.x=element_blank())

c=ggplot(plt, aes(x= SUBMISSION.DATE, y= INSTRUMENT_MODEL,shape= INSTRUMENT_MODEL))+geom_point()+theme(axis.line.y=element_blank())+labs(x="Individuals ranked by sequencing date", y='Sequencer')+ guides(shape=guide_legend(title="",label=F,override.aes = list(alpha = 0))) +scale_x_date(date_labels = "%Y",date_breaks = "1 year")+theme(axis.text.y=element_text(size=10,hjust=0))

d=plot_grid(a,b,c, nrow=3, labels=c('A','B',''),rel_heights=c(4,4,1),align='v')
ggsave("~/Documents/QualityPaper/MapQualOverTime.jpg",d,height=14,width=9)
