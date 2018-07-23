##Just the Stats
library(data.table)
args = commandArgs(trailingOnly=T)

#read input file #/Users/luke/Documents/Regression/CHR22.Format_JPT.csv
Input <- fread(args[1],sep=',')

Positions<-unique(Input $Pos)

chunks=split(Positions, ceiling(seq_along(Positions)/10000))
index <- 1
for(chunk in chunks){
print(paste('Progress :',index/length(chunks)*100,' %'))
subChunk<-Input[which(Input $Pos %in% chunk),]	
stat <- with(subChunk, by(subChunk,Pos, function(x) 
		glm(value ~ average_quality_of_mapped_bases, data=x, family=binomial)))
#extract P value
less = unique(subChunk[,c('Chr','Pos','Pop')])
less$p = -log10(as.numeric(as.character(sapply(stat,function(x) summary(x)$coef[2,4]))))
#write results to file
if(index == 1){
	write.table(less,file=args[2], quote=F, row.names=F, append=FALSE)}
if(index > 1){
    write.table(less,file=args[2], quote=F, row.names=F, append=TRUE)}
index <- index + 1
}

##code dump
# for f in `seq 14 21`; do echo $f; java -jar ~/bin/snpEff_latest_core/snpEff/SnpSift.jar extractFields /Users/luke/genomes/genomes/1kGP_4bedFiltered/1kGP_Chr${f}.4bed_filtered.vcf.gz CHROM POS "GEN[*].GT"   | sed 's/\//|/g ; s/0|0/0/g ; s/0|1/1/g ; s/1|0/1/g ; s/1|1/2/g ; /|/d'| tail -n +2 | gzip > /Users/luke/Documents/Regression/CHR${f}.Genotypes.txt.gz;  python /Users/luke/Documents/QualityPaper/GT_wideLong.py -i  /Users/luke/Documents/Regression/CHR${f}.Genotypes.txt.gz -o /Users/luke/Documents/Regression/CHR${f}.Format; done
