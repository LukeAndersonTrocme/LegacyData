##JPT Manhattan

dir='/Volumes/gravel-1/luke_projects/1000Genomes/Regression/'
out='/Volumes/gravel-1/luke_projects/1000Genomes/MeanDev/'


jptNames = list.files(path=dir,pattern='JPT', full.names = T)
jReg = do.call(rbind, lapply(jptNames, function(x) read.table(x, header=T, sep=' ')))
jReg<-jReg[which(jReg$Pos != 'Pos'),]
jReg$Chr<-as.numeric(as.character(jReg$Chr))
jReg$Pos<-as.numeric(as.character(jReg$Pos))
jReg$dev<-as.numeric(as.character(jReg$dev))

jReg$plog10 <- - pchisq(jReg$dev, 1, lower.tail=F, log.p=T)/log(10)

jReg$p <- pchisq(jReg$dev, 1, lower.tail=F)

jsig.6<-jReg[which(jReg$plog10 >= 6),]
jNotSig<-jReg[which(jReg$plog10 < 6),]

ggplot(jNotSig, aes(x=Pos, y = plog10, color=as.factor(Chr)))+
facet_grid(~Chr, scales='free_x', space='free_x', switch='x')+
scale_fill_manual (values=getPalette(colourCount))+
scale_y_continuous(expand=c(0,0))+
geom_point(alpha=0.3, size=1)+theme_classic()+
geom_point(data=jsig.6, aes(x=Pos, y = plog10,),color='black',shape=3)+
geom_hline(yintercept = 8, color='red')+
geom_hline(yintercept = 6, color='blue')+
labs(y='-log10(p)', x='Chromosome')+
theme(plot.title = element_text(hjust = 0.5), 
plot.subtitle = element_text(hjust = 0.5), 
axis.text.x=element_blank(), 
axis.ticks.x=element_blank(),
axis.line.x=element_blank(), 
strip.background = element_blank(),
strip.text.x = element_text(size = 6))+
guides(color=F)

ggsave('~/Documents/Regression/JPT_Manhattan.jpg', height=7, width=12)
