library(ggplot2)
#read colors per pop (used in 1kG Nature)
colors<-read.table("~/Dropbox/LukeTemp/NamePopBigPop_color.html.txt", comment.char='>',col.names=c('Pop','hex','BigPop'))
MyColour <- as.character(colors $hex)
names(MyColour) <- colors $Pop


pops<-list.files(path="/Users/luke/Documents/SigSnps",pattern=".frq$", full.names=T)

Freq<-data.frame()
for(f in seq(1,length(pops))){
pop<-unique(read.table(pops[f], row.names=NULL,header=T))
names(pop)=c('CHROM', 'POS', 'N_ALLELES','N_CHR', 'Y.ALLELE.FREQ','X.ALLELE.FREQ')
pop$Pop=substr(pops[f],31,33)
Freq<-rbind(Freq, pop)
}

Freq$X.ALLELE.FREQ<-gsub('T:','',Freq$X.ALLELE.FREQ)
Freq$X.ALLELE.FREQ<-gsub('C:','',Freq$X.ALLELE.FREQ)
Freq$X.ALLELE.FREQ<-gsub('G:','',Freq$X.ALLELE.FREQ)
Freq$X.ALLELE.FREQ<-gsub('A:','',Freq$X.ALLELE.FREQ)

perPop<-read.table('/Users/luke/Documents/SigSnps/PerPop.txt', col.names=c('CHROM','POS','Pop'))
perPop$Sig<-'Yes'

Freq<-merge(Freq,perPop, by=c('CHROM', 'POS','Pop'),all=T)

ordr<-c('FIN', 'GBR', 'CEU', 'IBS', 'TSI', 'CHS', 'CDX', 'CHB', 'JPT', 'KHV', 'GIH', 'STU', 'PJL', 'ITU', 'BEB', 'PEL', 'MXL', 'CLM', 'PUR', 'ASW', 'ACB', 'GWD', 'YRI', 'LWK', 'ESN', 'MSL')

Freq $ordered <- factor(Freq$Pop, levels =ordr)
#Freq[with(Freq, order(ordered),]


plt<-Freq[,c('CHROM','POS','X.ALLELE.FREQ','ordered','Sig')]
dodge <- position_dodge(width=2)  

ggplot(plt, aes(x=ordered, y= as.numeric(X.ALLELE.FREQ), color=ordered, group=1))+geom_hline(yintercept=0, color='grey70', linetype=1)+geom_hline(yintercept=0.5, color='grey70', linetype=3)+geom_line()+geom_point(shape=15)+scale_color_manual(values=MyColour)+scale_y_continuous(breaks=c(0,0.5))+facet_wrap(.~ POS,ncol=1)+theme_classic()+theme(strip.background=element_blank(), strip.text=element_blank(), axis.title=element_blank(), axis.ticks.y=element_blank(),axis.line=element_blank(), axis.text.x=element_text(angle=90, size=5, vjust=-0.05))+guides(color=F)+geom_point(data=plt[complete.cases(plt),], shape=1, size=4, color='black')

ggsave('~/Documents/QualityPaper/FreqPerSnp.jpg', height=10, width=3)