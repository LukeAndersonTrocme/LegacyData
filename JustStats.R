##Just the Stats
library(data.table)
args = commandArgs(trailingOnly=T)

#read input file
#Input <- fread('/Users/luke/Documents/Regression/CHR22.Format_JPT.csv',sep=',')
Input <- fread(paste("gunzip -c ",args[1]),sep=',', header=T)
Pop <-args[3]
Positions<-unique(Input $Pos)

chunks=split(Positions, ceiling(seq_along(Positions)/10000))
index <- 1
for(chunk in chunks){
print(paste(Pop,' Progress :',round(index/length(chunks)*100, digits=2),' %'))
subChunk<-Input[which(Input $Pos %in% chunk),]	
stat <- with(subChunk, by(subChunk,Pos, function(x) 
		glm(value ~ average_quality_of_mapped_bases, data=x, family=binomial)))
#extract deviance
less = unique(subChunk[,c('Chr','Pos','Pop')])
less$dev = as.numeric(as.character(sapply(stat,function(x) anova(x, test='Chi')[2,2])))

#write results to file
if(index == 1){
	write.table(less,file=args[2], quote=F, row.names=F, append=FALSE)}
if(index > 1){
    write.table(less,file=args[2], quote=F, row.names=F, append=TRUE)}
index <- index + 1
}

##code dump
# for f in `seq 1 22`; do echo $f; java -jar ~/bin/snpEff_latest_core/snpEff/SnpSift.jar extractFields /Users/luke/genomes/genomes/1kGP_4bedFiltered/1kGP_Chr${f}.4bed_filtered.vcf.gz CHROM POS "GEN[*].GT"   | sed 's/\//|/g ; s/0|0/0/g ; s/0|1/1/g ; s/1|0/1/g ; s/1|1/2/g ; /|/d'| tail -n +2 | gzip > /Users/luke/Documents/Regression/CHR${f}.Genotypes.txt.gz;  python /Users/luke/Documents/QualityPaper/GT_wideLong.py -i  /Users/luke/Documents/Regression/CHR${f}.Genotypes.txt.gz -o /Users/luke/Documents/Regression/CHR${f}.Format; done
