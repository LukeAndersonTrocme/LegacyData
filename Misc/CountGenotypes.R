library(data.table)
args = commandArgs(trailingOnly=T)

ColNames = read.table('/Users/luke/genomes/genomes/hg19/Genotypes/ColNames.txt', header=F)
fileName = paste("gzcat ",args[1], sep='')
print(fileName)
genotypes <- fread(fileName, sep='\t', header=F)
names(genotypes)= as.character(unlist(ColNames))

chunks = split(genotypes, cumsum((1:nrow(genotypes)-1)%%10000==0))
i=1
for(chunk in chunks){

chunk $Count = rowSums(chunk[,-(1:3)])

	#write results to file
	if(i == 1){
		write.table(chunk[,c('CHROM','POS','ID','Count')],file=args[2], 
		quote=F, row.names=F, col.names=T, append=FALSE)}
	if(i > 1){
    	write.table(chunk[,c('CHROM','POS','ID','Count')],file=args[2], 
    	quote=F, row.names=F, col.names=F, append=TRUE)}
    
    i=i+1
    	}