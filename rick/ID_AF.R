#Taking Rick's Code and putting it in an Rscript

#load libaries
library(Matrix)
library(irlba)
library(threejs)
library(tidyverse)
library(lfa)
library(furrr)
library(glm2)
library(cowplot)
library(drake)

#For running in terminal
args = commandArgs(trailingOnly=T)
Chr <- args[1]
#Chr <- '22'
#set directories and paths
GenoPath='/Users/luke/genomes/genomes/hg19/phase3/'
Dir = '/Users/luke/Documents/QualityPaper/rick/'
setwd(Dir)

#Get ID names
if(!file.exists("ids.txt")){
ids <- readLines(pipe(paste("gzcat ",GenoPath,"ALL.chr",Chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz  |
                      sed -n /^#CHROM/p | 
                      tr '\t' '\n' | 
                      tail -n +10",sep='')))

writeLines(ids, "ids.txt")}

if(file.exists("ids.txt")){
  ids <- unlist(read.table("ids.txt"))
  head(ids)
}

#Get Quality Scores
read_basfiles <- function(bas_filepath) {
  
  bas_file <- read_tsv(bas_filepath,
                  col_types = cols_only('sample' = 'c',        
                                       'average_quality_of_mapped_bases' = col_guess())) 
  return( bas_file)
}

bas_files <- list.files(path = "/Users/luke/genomes/bas_file/1000GenomesBAS/",
                        pattern="*bam.bas", full.names = TRUE)

bas_dt <- 
  map_dfr(bas_files,read_basfiles ) %>%
  rename(avg_qual="average_quality_of_mapped_bases") %>%
  group_by(sample) %>% summarize(avg_qual = mean(avg_qual))

#keep samples in 1kGP VCF files
bas_dt <- bas_dt[which((bas_dt$sample %in% ids)==T),]

#get populations of each ID
ped <- 
  read_tsv(url("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped"))%>% 
  rename(sample="Individual ID")

bas_dt <- inner_join( ped  %>% 
                       select(sample, Population),
                       bas_dt,
                     by="sample")
                     
#Load Genotype Data
#Load VCF file into a sparse matrix                     
p <- pipe(paste("gzcat ",GenoPath,"ALL.chr",Chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz  |
          sed /^#/d  |
          cut  -f '10-' |
          ./a.out | 
          cut -f '1-2'",sep=''))
x <- read.table(p,
                colClasses=c("integer","integer"), 
                fill=TRUE,
                row.names=NULL)

# Convert to a sparse matrix of people (rows) x variant (columns)
GT <- sparseMatrix(i=x[,2], j=x[,1], x=1.0)

# Inspect the dimensions of this matrix
print(dim(GT))                     

# PCA compute the first three principal component vectors
cm <- colMeans(GT)
p <- irlba(GT,
           nv=3,
           nu=3,
           tol=0.1,
           center=cm)
           
#Add to dataset
colnames(p$u) <- paste0("pc",1:3)
bas_dt <- bind_cols(bas_dt, 
                    as_tibble(p$u))
                     
# Logistic factor analysis
GT_mtx <- t(as.matrix(GT[,1:50000]))
LF <- lfa(GT_mtx, 4)
print(dim(LF))
#Add to dataset
bas_dt <- bind_cols(bas_dt,
                    tibble(lf1=LF[,1], lf2=LF[,2],lf3=LF[,3]))

#compute Individual Specific Allele Frequency
sub_ial <- af(GT_mtx, LF)
                    
lfa <- ggplot(bas_dt,
       aes(lf1, lf2, color=Population)) +
  geom_point() +
  theme_bw()

pca <- ggplot(bas_dt,
       aes(pc1, pc3, color=Population)) +
  geom_point() +
  theme_bw()
     
plot_grid(pca, lfa)
 
to_datalist <- function(i, data_genotype, data_iaf, qual){
  
    data <- tibble(genotype=data_genotype[i,],
               iaf=data_iaf[i,],
               qual=qual)
               
    return(data)       
}

fit_model_1 <- function(data){
	
	fit <- glm2(genotype ~ bas_dt$avg_qual + offset(iaf), 
				family ="binomial", data=data)
	
	return(fit)
}

plot_data <- function(data){
  
      p <- ggplot(data,
       aes(qual, genotype, color=iaf)) +
      geom_point() +
      theme_bw()
  return(p)
  
}

loci <- seq(1,50000)
names(loci) <- loci

data_list <- map(loci, 
                 to_datalist,
                 data_genotype=GT_mtx,
                 data_iaf=sub_ial)

#save(data_list, file = '~/Documents/QualityPaper/rick/data_list.temp')
#load('~/Documents/QualityPaper/rick/data_list.temp')  
               
future::plan(sequential)

dt <- head(data_list,100)

fits <- future_map(dt, fit_model_1)

coefLF = lapply(fits, function(x) coef(x)[[2]])  
                                                              