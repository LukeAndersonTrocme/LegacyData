##Manhattan
fileNames = list.files(path='/Volumes/gravel/luke_projects/1000Genomes/Regression/', pattern="*.Regression_*.csv", full.names = T)
Reg = do.call(rbind, lapply(fileNames, function(x) read.table(x, header=T)))

Reg<-Reg[which(Reg$Pos != 'Pos'),]
Reg$Chr<-as.numeric(as.character(Reg$Chr))
Reg$Pos<-as.numeric(as.character(Reg$Pos))
Reg$dev<-as.numeric(as.character(Reg$dev))
Reg$Pop=NULL
Reg$p=NULL

Reg <- Reg %>% group_by(Chr,Pos) %>% summarize(SumDev = sum(dev), MeanDev = mean(dev), Count = n())

Reg$plog10 <- - pchisq(Reg$dev, 1, lower.tail=F, log.p=T)/log(10)
sig.6<-Reg[which(Reg$plog10 > 6),]
NotSig<-Reg[which(Reg$plog10 < 6),]

ggplot(NotSig, aes(x=Pos, y = plog10, color=as.factor(Chr)))+
facet_grid(~Chr, scales='free_x', space='free_x', switch='x')+
scale_fill_manual (values=getPalette(colourCount))+
scale_y_continuous(expand=c(0,0))+
geom_point(alpha=0.3, size=1)+theme_classic()+
geom_point(data=sig.6, aes(x=Pos, y = plog10,),color='black',shape=3)+
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
