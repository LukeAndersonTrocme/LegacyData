pops<-fread('~/Dropbox/LukeTemp/Misc/Name.BigPop.Pop.txt', col.names=c('Name','BigPop','Pop'), header=F)

for(con in unique(pops$BigPop)){
	p=pops[which(pops$BigPop==con),]
	print(as.vector(unique(p$Pop)))
	#pp<-cbind(p$Name, p$Name)
	#write.table(pp,paste('~/genomes/genomes/PopNames/',con,'_names.txt', sep=''), quote=F, col.names=F, row.names=F)
}

a<-unique(pops[,c("Pop","BigPop")])
write.table(a, '~/Dropbox/LukeTemp/Misc/BigPop.Pop.txt', quote=F, row.names=F)

as.character(a[which(a$Pop=='JPT'),"BigPop"])