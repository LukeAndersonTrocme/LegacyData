library(data.table)
library(ggplot2)
library(cowplot)

Pops <- fread('~/genomes/genomes/PopNames/NamePop.txt', header=F, nrows=26)

Qual <- fread('~/genomes/genomes/misc/1000GenomesQual_INT.pheno')

# GET EQUATION AND R-SQUARED AS STRING
# SOURCE: http://goo.gl/K4yh

lm_eqn <- function(df){
	names(df)=c('x','y')
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~pv,
         list(a = format(coef(m)[1], digits = 2), 
              b = format(coef(m)[2], digits = 2), 
             r2 = format(summary(m)$r.squared, digits = 3),
             pv = format(summary(m)$coefficients[2,4], digits = 3)))
    as.character(as.expression(eq));                 
}

p1 <- p + geom_text(x = 25, y = 300, label = lm_eqn(df), parse = TRUE)


for(i in seq(1,26)){
	pop = as.character(Pops$V1[i])
	PC = fread(paste('~/Documents/PCAperPop/',pop,'.eigenvec',sep=''))

	plt = merge(Qual, PC, by=c('V1','V2'))
	p1 = ggplot(plt, aes(x=V3.y,y=V3.x))+geom_point()+
	geom_smooth(method='lm',formula=y~x, se=F)+
	theme_classic()+theme(plot.title = element_text(hjust = 0.5))+
	labs(y='Normalized Quality',x='PC1',title=pop)+
	annotate("text", x = -0.25, y = 2.5, 
		label = lm_eqn(plt[,c('V3.y','V3.x')]), size = 2, parse=TRUE)
	
	p2 = ggplot(plt, aes(x=V4,y=V3.x))+geom_point()+
	geom_smooth(method='lm',formula=y~x, se=F)+
	theme_classic()+
	labs(y='Normalized Quality',x='PC2')+
	annotate("text", x = -0.25, y = 2.5, 
		label = lm_eqn(plt[,c('V4','V3.x')]), size = 2, parse=TRUE)
		
	p3 = plot_grid(p1,p2,ncol = 1)
	
	ggsave(paste('~/Documents/PCAperPop/PC_Qual',pop,'.jpg',sep = ''), p3, height=7, width=4)		
}

#annotate("text",x=0.75,y=0.25, label=deparse(lm_eqn(plt[,c('V3.x','V3.y')])), parse=TRUE)



for(i in seq(1,26)){
	pop = as.character(Pops$V1[i])
	PC = fread(paste('~/Documents/PCAperPop/',pop,'.eigenvec',sep=''))
	ggplot(PC, aes(V3,V4))+geom_point()+theme_classic()
	ggsave(paste('~/Documents/PCAperPop/PC1vPC2',pop,'.jpg',sep = ''), height=4, width=4)
	}