##rejected from Manhattan


sig.6<-Reg[which((Reg$plog10 >= 6)&(Reg$plog10 < 20)&(Reg$Count > 1)),]
sig.20<-Reg[which((Reg$plog10 >= 20)&(Reg$Count > 1)),]
sig.20$plog10=20
NotSig<-Reg[which((Reg$plog10 < 6)&(Reg$Count > 1)),]

write.table(Reg[which((Reg$plog10 >= 6)&(Reg$Count > 1)),], '~/Documents/QualityPaper/Significant6SNPs.txt', quote=F, row.names=F)

ggplot(NotSig, aes(x=Pos, y = plog10, color=as.factor(Chr)))+
facet_grid(~Chr, scales='free_x', space='free_x', switch='x')+
scale_fill_manual (values=getPalette(colourCount))+
scale_y_continuous(expand=c(0,0))+
geom_point(alpha=0.3, size=1)+theme_classic()+
geom_point(data=sig.6, aes(x=Pos, y = plog10,),color='black',shape=3)+
geom_point(data=sig.20, aes(x=Pos, y = plog10,),color='black',shape=1)+
labs(y='-log10(p)', x='Chromosome')+
theme(plot.title = element_text(hjust = 0.5), 
plot.subtitle = element_text(hjust = 0.5), 
axis.text.x=element_blank(), 
axis.ticks.x=element_blank(),
axis.line.x=element_blank(), 
strip.background = element_blank(),
strip.text.x = element_text(size = 6))+
guides(color=F)+expand_limits(y=c(0,21))
ggsave('~/Documents/QualityPaper/ManhattanPlot.jpg', height=5, width=10)

