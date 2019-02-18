##Compare GWAS p-values
library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)
colors <- read.table(
	"~/Dropbox/LukeTemp/NamePopBigPop_color.html.txt", comment.char='>', 
	col.names=c('Pop','hex','BigPop'))
MyColour <- as.character(colors $hex)
names(MyColour) <- colors $Pop

rsID <- fread(
	'/Users/luke/genomes/genomes/1kGP_4bedFiltered/1kGP_GenomeWide_Chr.Pos.rsID.txt',
	col.names=c('Chr','Pos','rsID'))

gbr <- fread(
	'/Users/luke/Documents/QualityPaper/sig/_Format_sig.csv',
	 col.names=c('Chr','Pos','ID','Genotype','Pop','Quality', 'Count'))

pubs <- fread(
	'~/Dropbox/LukeTemp/Significant6SNPs_adjusted_published.txt',
	col.names=c('row', 'Author','Journal','link','rsID','AF','pubPval'))	

pubs$pubPval <- NULL

plt<-unique(left_join(gbr, rsID))

plt<-unique(left_join(plt, pubs))
plt<-plt[,complete.cases(plt$Genotype)]
ggplot(plt, aes(y=Genotype, x=Quality))+facet_wrap(Journal~Pos)+geom_point()+geom_smooth(method = "glm",     method.args = list(family = "binomial"), se = FALSE,size=0.5)